#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <vector>
#include <sstream>
#include <algorithm>
#include <fstream> // for writing results to a file
#include <numeric>

#include "../include/math.hpp"
#include "../include/mc_types.hpp"
#include "../include/option.hpp"
#include "../include/stepper.hpp"
#include "../include/mc_engine.hpp"
#include "../include/convergence.hpp"

// Small stable hash to tag RNG streams by (method, optType, AV, CV, Npaths, Nsteps)
inline uint64_t rng_tag(const char* method, const std::string& optType,
                        bool av, bool cv, std::size_t Npaths, std::size_t Nsteps)
{
    uint64_t h = 1469598103934665603ULL; // FNV offset basis
    auto mix = [&](uint64_t x){ h ^= x; h *= 1099511628211ULL; };
    for (const char* p = method; *p; ++p) mix((unsigned char)*p);
    for (char c : optType) mix((unsigned char)c);
    mix(av); mix(cv); mix(Npaths); mix(Nsteps);
    return h;
}


// Generic XY CSV writer: header "x_name,y_name" then pairs
static void write_xy_csv(const std::string& filename,
                         const std::string& x_name,
                         const std::string& y_name,
                         const std::vector<double>& xs,
                         const std::vector<double>& ys)
{
    if (xs.size() != ys.size())
        throw std::runtime_error("write_xy_csv: xs.size() != ys.size()");
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("write_xy_csv: cannot open file: " + filename);
    out << x_name << "," << y_name << "\n";
    for (std::size_t i = 0; i < xs.size(); ++i)
        out << xs[i] << "," << ys[i] << "\n";
}

// Export 10 (or user-chosen) sequences; each sequence runs N=1..Nmax and records the MC price.
// We write 4 CSVs: <base>_none.csv, <base>_av.csv, <base>_cv.csv, <base>_both.csv
// Each CSV: header "N,seq,price" with seq in [0..numSeq-1].
#include <fstream>
#include <memory>

// Export numSeq sequences; each sequence runs N = Nmin..Nmax (step = stride) and records the MC price
// using GBM *terminal* (exact) draws. Writes 4 CSVs:
//   <base>_none.csv, <base>_av.csv, <base>_cv.csv, <base>_both.csv
// Each CSV: "N,seq,price"
// NOTE: In terminal mode, CV is inert (no discretization error to correlate with a control).
inline std::unique_ptr<std::ofstream> open_csv_(const std::string& path) {
    auto f = std::make_unique<std::ofstream>(path);
    if (!*f) throw std::runtime_error("Cannot open CSV: " + path);
    *f << "N,seq,price\n";
    return f;
}

static void export_price_sequences_terminal_gbm(const std::string& base_tag,
                                                EuropeanOption<double>& opt,
                                                const std::string& optionType,
                                                std::size_t Nmax,
                                                int numSeq,
                                                std::size_t Nmin = 2,
                                                std::size_t stride = 1)
{
    if (Nmax < 2) Nmax = 2;
    if (Nmin < 2) Nmin = 2;
    if (stride == 0) stride = 1;

    auto f_none = open_csv_(base_tag + "_none.csv");
    auto f_av   = open_csv_(base_tag + "_av.csv");
    auto f_cv   = open_csv_(base_tag + "_cv.csv");
    auto f_both = open_csv_(base_tag + "_both.csv");

    GBM<double> gbm;

    for (int seq = 0; seq < numSeq; ++seq) {
        // If your MC engine supports seeding, set a different seed here for independence
        // uint64_t seed = 1234567ULL + 7919ULL * seq;
        rng_reseed_thread(static_cast<uint64_t>(seq));
        // Before each variant, reset to the same (seq,N)-tagged stream:
        auto tag = [](int seq, std::size_t N, int var){ 
            // small hash / tag composer; any stable mapping works
            return (uint64_t(seq) << 32) ^ (uint64_t(N) << 2) ^ uint64_t(var);
        };

        for (std::size_t N = Nmin; N <= Nmax; N += stride) {
            // NONE
            rng_reseed_thread(tag(seq, N, 0));
            {
                MonteCarlo<double, EuropeanOption<double>, GBM<double>> mc(opt, gbm);
                mc.useAntitheticVariates(false);
                mc.useControlVariates(false); // inert here
                auto r = mc.priceTerminal(optionType, N);
                *f_none << N << "," << seq << "," << r.price << "\n";
            }
            // AV
            rng_reseed_thread(tag(seq, N, 1));
            {
                MonteCarlo<double, EuropeanOption<double>, GBM<double>> mc(opt, gbm);
                mc.useAntitheticVariates(true);
                mc.useControlVariates(false); // inert
                auto r = mc.priceTerminal(optionType, N);
                *f_av << N << "," << seq << "," << r.price << "\n";
            }
            // CV (inert with terminal draw; kept for completeness)
            rng_reseed_thread(tag(seq, N, 2));
            {
                MonteCarlo<double, EuropeanOption<double>, GBM<double>> mc(opt, gbm);
                mc.useAntitheticVariates(false);
                mc.useControlVariates(true);
                auto r = mc.priceTerminal(optionType, N);
                *f_cv << N << "," << seq << "," << r.price << "\n";
            }
            // BOTH (AV effective; CV inert)
            rng_reseed_thread(tag(seq, N, 3));
            {
                MonteCarlo<double, EuropeanOption<double>, GBM<double>> mc(opt, gbm);
                mc.useAntitheticVariates(true);
                mc.useControlVariates(true);
                auto r = mc.priceTerminal(optionType, N);
                *f_both << N << "," << seq << "," << r.price << "\n";
            }
        }
    }
}



// Simple XY CSV writer: header + pairs.
// static void write_xy_csv(const std::string& filename,
//                          const std::string& x_name,
//                          const std::string& y_name,
//                          const std::vector<double>& xs,
//                          const std::vector<double>& ys)
// {
//     if (xs.size() != ys.size())
//         throw std::runtime_error("write_xy_csv: xs.size() != ys.size()");
//     std::ofstream out(filename);
//     if (!out) throw std::runtime_error("write_xy_csv: cannot open " + filename);
//     out << x_name << "," << y_name << "\n";
//     for (std::size_t i = 0; i < xs.size(); ++i) out << xs[i] << "," << ys[i] << "\n";
// }

// Write flat (path-major) vector as a tidy CSV: path_id, step, t, S
static void write_paths_csv(const std::string& filename,
                            const std::vector<double>& flat,
                            std::size_t Npaths,
                            std::size_t Nsteps,
                            double T)
{
    if (Npaths == 0 || Nsteps == 0 || flat.size() != Npaths * Nsteps) {
        throw std::runtime_error("write_paths_csv: size mismatch (Npaths*Nsteps != flat.size()).");
    }

    std::ofstream out(filename);
    if (!out) throw std::runtime_error("write_paths_csv: cannot open file: " + filename);

    out << "path_id,step,t,S\n";

    const double dt = T / static_cast<double>(Nsteps);
    for (std::size_t p = 0; p < Npaths; ++p) {
        for (std::size_t k = 0; k < Nsteps; ++k) {
            const std::size_t idx = p * Nsteps + k; // path-major
            const double t = (k + 1) * dt;
            out << p << ',' << k << ',' << t << ',' << flat[idx] << '\n';
        }
    }
}

// Convenience wrapper: run MC full-path pricing to fill paths and dump them.
// Now Npaths and Nsteps are explicit user arguments.
template <typename Step>
static void dump_mc_paths_csv(const std::string& filename,
                              EuropeanOption<double>& opt,
                              const Step& stepper,
                              const std::string& optionType,
                              std::size_t Npaths,
                              std::size_t Nsteps,
                              bool useAV = false,
                              bool useCV = false)
{
    MonteCarlo<double, EuropeanOption<double>, Step> mc(opt, stepper);
    mc.useAntitheticVariates(useAV);
    mc.useControlVariates(useCV);

    // Run with user’s Nsteps
    (void) mc.pricePath(optionType, Npaths, Nsteps);  // <- ensure your pricePath overload accepts Nsteps

    const std::vector<double>& flat = mc.getPricePathsFlat();

    if (flat.size() != Npaths * Nsteps) {
        if (flat.size() == Npaths * (Nsteps + 1)) {
            std::vector<double> trimmed; trimmed.reserve(Npaths * Nsteps);
            for (std::size_t p = 0; p < Npaths; ++p) {
                const std::size_t base = p * (Nsteps + 1);
                for (std::size_t k = 1; k <= Nsteps; ++k)
                    trimmed.push_back(flat[base + k]);
            }
            write_paths_csv(filename, trimmed, Npaths, Nsteps, opt.getMaturityTime());
            return;
        }
        throw std::runtime_error("Unexpected flat paths size in MonteCarlo; adjust indexing.");
    }

    write_paths_csv(filename, flat, Npaths, Nsteps, opt.getMaturityTime());
}

// === STREAMING STRONG CONVERGENCE ===
// Calls ca.strongConvergence on growing prefixes of Ms and prints the newly-available row each time.
template <typename Analyzer>
static void stream_strong_conv(const char* tag, Analyzer& ca,
                               std::size_t numPaths, const std::vector<std::size_t>& Ms)
{
    using std::cout; using std::left; using std::setw; using std::fixed; using std::setprecision;

    // Header once
    cout << "\n== STRONG convergence (" << tag << ") ==\n";
    cout << std::string(70, '-') << "\n";
    cout << left
         << setw(12) << "Steps (M)"
         << setw(12) << "dt"
         << setw(22) << "E[|ΔS_T|]"
         << setw(22) << "RMS(ΔS_T)" << std::endl;
    cout << std::string(70, '-') << std::endl;

    StrongConvergenceResult Rfull;
    for (std::size_t i = 0; i < Ms.size(); ++i) {
        std::vector<std::size_t> Ms_prefix(Ms.begin(), Ms.begin() + i + 1); // this is the prefix, meaning
        auto R = ca.strongConvergence(numPaths, Ms_prefix);

        // Print ONLY the last row (new result)
        const std::size_t k = R.M.size() - 1;
        cout << left
             << setw(12) << R.M[k]
             << setw(12) << fixed << setprecision(6) << R.dt[k]
             << setw(22) << fixed << setprecision(6) << R.meanAbs[k]
             << setw(22) << fixed << setprecision(6) << R.rms[k]
             << std::endl;

        // keep last full result for slope reporting at the end
        if (i + 1 == Ms.size()) Rfull = std::move(R);
    }

    cout << std::string(70, '-') << std::endl;
    cout << "Estimated order on E[|ΔS_T|]: " << fixed << setprecision(6) << Rfull.slope_meanAbs
         << "   |   on RMS: " << fixed << setprecision(6) << Rfull.slope_rms << std::endl;

        // Export for plotting: x = dt, y = E[|ΔS_T|]  (you can also export RMS similarly)
    try {
        write_xy_csv("strong_meanabs.csv", "dt", "E_abs_dS",
                     Rfull.dt, Rfull.meanAbs);
        // Optional: RMS too
        // write_xy_csv("strong_rms.csv", "dt", "RMS_dS", Rfull.dt, Rfull.rms);
        std::cout << "[csv] strong_meanabs.csv written.\n";
    } catch (const std::exception& ex) {
        std::cerr << "[csv] strong export error: " << ex.what() << "\n";
    }

}

// === STREAMING WEAK CONVERGENCE ===
template <typename Analyzer>
static void stream_weak_conv(const char* tag, Analyzer& ca, const std::string& optionType,
                             std::size_t numPaths, const std::vector<std::size_t>& Ms)
{
    using std::cout; using std::left; using std::setw; using std::fixed; using std::scientific; using std::setprecision;

    cout << "\n== WEAK convergence (" << tag << ") ==\n";
    cout << std::string(76, '-') << "\n";
    cout << left
         << setw(12) << "Steps (M)"
         << setw(12) << "dt"
         << setw(26) << "Bias"
         << setw(26) << "SE(Bias)" << std::endl;
    cout << std::string(76, '-') << std::endl;

    WeakConvergenceResult Rfull;
    for (std::size_t i = 0; i < Ms.size(); ++i) {
        std::vector<std::size_t> Ms_prefix(Ms.begin(), Ms.begin() + i + 1);
        auto R = ca.weakConvergence(optionType, numPaths, Ms_prefix);

        const std::size_t k = R.M.size() - 1;
        cout << left
             << setw(12) << R.M[k]
             << setw(12) << fixed      << setprecision(6) << R.dt[k]
             << setw(26) << scientific << setprecision(3) << R.bias[k]
             << setw(26) << scientific << setprecision(3) << R.seBias[k]
             << std::endl;

        if (i + 1 == Ms.size()) Rfull = std::move(R);
    }

    cout << std::string(76, '-') << std::endl;
    cout << "Estimated WEAK order (slope of log|bias| vs log dt): "
         << fixed << setprecision(6) << Rfull.slope_bias << std::endl;
    
    // Export for plotting: x = dt, y = |bias|
    std::vector<double> absBias = Rfull.bias;
    for (auto& b : absBias) b = std::abs(b);
    try {
        write_xy_csv("weak_bias.csv", "dt", "abs_bias", Rfull.dt, absBias);
        std::cout << "[csv] weak_bias.csv written.\n";
    } catch (const std::exception& ex) {
        std::cerr << "[csv] weak export error: " << ex.what() << "\n";
    }

}

// === STREAMING MC STATISTICAL CONVERGENCE ===
template <typename Analyzer>
static void stream_mc_stat_conv(const char* tag, Analyzer& ca, const std::string& optionType,
                                const std::vector<std::size_t>& Ns, bool useAV, bool useCV)
{
    using std::cout; using std::left; using std::setw; using std::fixed; using std::setprecision;

    cout << "\n== MC statistical convergence (" << tag << ") ==  [expect slope ≈ -0.5]\n";
    cout << std::string(76, '-') << "\n";
    cout << left
         << setw(12) << "Paths (N)"
         << setw(16) << "Price"
         << setw(16) << "StdErr"
         << setw(16) << "CI95 width" << std::endl;
    cout << std::string(76, '-') << std::endl;

    MCStatConvergence Rfull;
    for (std::size_t i = 0; i < Ns.size(); ++i) {
        std::vector<std::size_t> Ns_prefix(Ns.begin(), Ns.begin() + i + 1);
        auto R = ca.mcConvergence(optionType, Ns_prefix, useAV, useCV);

        const std::size_t k = R.Ns.size() - 1;
        const double ci95 = 1.96 * 2.0 * R.stdErrs[k];
        cout << left
             << setw(12) << R.Ns[k]
             << setw(16) << fixed << setprecision(6) << R.prices[k]
             << setw(16) << fixed << setprecision(6) << R.stdErrs[k]
             << setw(16) << fixed << setprecision(6) << ci95
             << std::endl;

        if (i + 1 == Ns.size()) Rfull = std::move(R);
    }

    cout << std::string(76, '-') << std::endl;
    cout << "Estimated order (slope of log SE vs log N): "
         << fixed << setprecision(4) << Rfull.slope_loglog << std::endl;

        // Export for plotting: x = N, y = StdErr
    // (Convert Ns to double for CSV)
    std::vector<double> Ns_d(Rfull.Ns.begin(), Rfull.Ns.end());
    try {
        write_xy_csv("mc_stat_stderr.csv", "N", "StdErr", Ns_d, Rfull.stdErrs);
        std::cout << "[csv] mc_stat_stderr.csv written.\n";
    } catch (const std::exception& ex) {
        std::cerr << "[csv] mc_stat export error: " << ex.what() << "\n";
    }

}


// ---------- Lightweight table reporter (unchanged) ----------
struct Reporter {
    bool header_printed = false;
    void print_header_once() {
        if (header_printed) return; header_printed = true;
        using std::setw; using std::left; using std::right;
        std::cout << "\n";
        std::cout << std::string(118, '=') << "\n";
        std::cout << left
                  << setw(10)  << "Method"
                  << setw(6)   << "Type"
                  << setw(6)   << "AV"
                  << setw(6)   << "CV"
                  << right
                  << setw(12)  << "Paths"
                  << setw(10)  << "Steps"
                  << setw(16)  << "Price"
                  << setw(16)  << "95% CI"
                  << setw(18)  << "Time (ms)" << "\n";
        std::cout << std::string(118, '-') << "\n";
    }
    void print_row(const char* method, const std::string& optType, bool av, bool cv,
                   std::size_t Npaths, std::size_t Nsteps, double price, double ci95, double ms)
    {
        using namespace std;
        if (!header_printed) print_header_once();
        cout << fixed << setprecision(6)
             << left
             << setw(10) << method
             << setw(6)  << optType
             << setw(6)  << (av ? "on" : "off")
             << setw(6)  << (cv ? "on" : "off")
             << right
             << setw(12) << Npaths
             << setw(10) << Nsteps
             << setw(16) << price
             << setw(16) << ci95
             << setw(18) << ms << "\n";
    }
    void finish() { std::cout << std::string(118, '=') << "\n"; }
};

inline void print_scenario_banner(const std::string& tag,
                                  double S0, double K, double r, double sigma, double T,
                                  double bs_call, double bs_put)
{
    std::ostringstream lineParams, lineTitle, lineBS;
    lineParams << std::fixed
               << "S0=" << std::setprecision(6) << S0 << "  "
               << "K=" << std::setprecision(6) << K << "  "
               << "S0/K=" << std::setprecision(6) << (S0 / K) << "  "
               << "r=" << std::setprecision(6) << r << "  "
               << "sigma=" << std::setprecision(6) << sigma << "  "
               << "T=" << std::setprecision(6) << T;
    lineTitle << " SCENARIO: " << tag << "  ";
    lineBS << " BS Call=" << std::fixed << std::setprecision(6) << bs_call
           << "    BS Put=" << bs_put;

    const std::size_t w = std::max({ lineTitle.str().size(), lineParams.str().size(), lineBS.str().size() });
    const std::string eq(w, '='), dash(w, '-');
    std::cout << "\n" << eq << "\n" << lineTitle.str() << "\n" << dash << "\n"
              << lineParams.str() << "\n" << dash << "\n"
              << lineBS.str() << "\n" << eq << "\n";
}

// Run one pricing row (unchanged core)
template <typename Step>
void run_fullpath_row(const char* methodName, EuropeanOption<double>& opt, const Step& stepper,
                      const std::string& optType, bool useAV, bool useCV,
                      std::size_t Npaths, std::size_t Nsteps, Reporter& rep)
{
    MonteCarlo<double, EuropeanOption<double>, Step> mc(opt, stepper);
    mc.useAntitheticVariates(useAV);
    mc.useControlVariates(useCV);
    const auto t0 = std::chrono::steady_clock::now();
    auto r = mc.pricePath(optType, Npaths);
    const auto t1 = std::chrono::steady_clock::now();
    const double ms  = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double ci95= 1.96 * r.stdDev;
    rep.print_row(methodName, optType, useAV, useCV, Npaths, Nsteps, r.price, ci95, ms);
}

template <typename Step>
void run_terminal_row(const char* methodName, EuropeanOption<double>& opt, const Step& stepper,
                      const std::string& optType, bool useAV, bool useCV, std::size_t Npaths, Reporter& rep)
{
    MonteCarlo<double, EuropeanOption<double>, Step> mc(opt, stepper);
    mc.useAntitheticVariates(useAV);
    mc.useControlVariates(useCV);
    const auto t0 = std::chrono::steady_clock::now();
    auto r = mc.priceTerminal(optType, Npaths);
    const auto t1 = std::chrono::steady_clock::now();
    const double ms  = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double ci95= 1.96 * r.stdDev;
    rep.print_row(methodName, optType, useAV, useCV, Npaths, /*Steps*/1, r.price, ci95, ms);
}

template <typename Step>
void run_fullpath_groups_for_type(const char* methodName, EuropeanOption<double>& opt, const Step& stepper,
                                  const std::string& optType, std::size_t Npaths, Reporter& rep)
{
    const std::size_t Nsteps = std::max<std::size_t>(1, static_cast<std::size_t>(TRADING_DAYS_PER_YEAR * opt.getMaturityTime()));

    rng_reseed_thread(rng_tag(methodName, optType, false, false, Npaths, Nsteps));
    run_fullpath_row(methodName, opt, stepper, optType, false, false, Npaths, Nsteps, rep);

    rng_reseed_thread(rng_tag(methodName, optType, true,  false, Npaths, Nsteps));
    run_fullpath_row(methodName, opt, stepper, optType, true,  false, Npaths, Nsteps, rep);

    rng_reseed_thread(rng_tag(methodName, optType, false, true,  Npaths, Nsteps));
    run_fullpath_row(methodName, opt, stepper, optType, false, true,  Npaths, Nsteps, rep);

    rng_reseed_thread(rng_tag(methodName, optType, true,  true,  Npaths, Nsteps));
    run_fullpath_row(methodName, opt, stepper, optType, true,  true,  Npaths, Nsteps, rep);
}


template <typename Step>
void run_terminal_groups_for_type(const char* methodName, EuropeanOption<double>& opt, const Step& stepper,
                                  const std::string& optType, std::size_t Npaths, Reporter& rep)
{
    rng_reseed_thread(rng_tag(methodName, optType, false, false, Npaths, 1));
    run_terminal_row(methodName, opt, stepper, optType, false, false, Npaths, rep);

    rng_reseed_thread(rng_tag(methodName, optType, true,  false, Npaths, 1));
    run_terminal_row(methodName, opt, stepper, optType, true,  false, Npaths, rep);

    rng_reseed_thread(rng_tag(methodName, optType, false, true,  Npaths, 1));
    run_terminal_row(methodName, opt, stepper, optType, false, true,  Npaths, rep);

    rng_reseed_thread(rng_tag(methodName, optType, true,  true,  Npaths, 1));
    run_terminal_row(methodName, opt, stepper, optType, true,  true, Npaths, rep);
}


// ---------- MAIN ----------
int main()
{
    rng_set_seed(42);

    struct {
        double S0 = 100.0;
        double r  = 0.05;
        double v  = 0.20;
        double T  = 1.0;
        std::size_t Npaths = 1000000;
    } cfg;

    EuropeanOption<double> euro(cfg.S0, /*K=*/cfg.S0, cfg.r, cfg.v, cfg.T);
    GBM<double>      gbm;
    Euler<double>    euler;
    Milstein<double> milstein;

    const double bs_call = euro.closedForm("call");
    const double bs_put  = euro.closedForm("put");
    print_scenario_banner("ATM",
                          euro.getSpotPrice(), euro.getStrikePrice(),
                          euro.getRiskFreeRate(), euro.getVolatility(), euro.getMaturityTime(),
                          bs_call, bs_put);

    // ================= Full-path (discretized) =================
    std::cout << "\n********** FULL-PATH (discretized) MONTE CARLO **********\n";
    {
        Reporter rep;
        rep.print_header_once();
        run_fullpath_groups_for_type("GBM",      euro, gbm,      "call", cfg.Npaths, rep);
        run_fullpath_groups_for_type("Euler",    euro, euler,    "call", cfg.Npaths, rep);
        run_fullpath_groups_for_type("Milstein", euro, milstein, "call", cfg.Npaths, rep);
        run_fullpath_groups_for_type("GBM",      euro, gbm,      "put",  cfg.Npaths, rep);
        run_fullpath_groups_for_type("Euler",    euro, euler,    "put",  cfg.Npaths, rep);
        run_fullpath_groups_for_type("Milstein", euro, milstein, "put",  cfg.Npaths, rep);
        rep.finish();
    }

    // ================= Terminal (single-step, exact GBM) =================
    std::cout << "\n********** TERMINAL (single-step, exact GBM) **********\n";
    {
        Reporter rep; rep.print_header_once();
        run_terminal_groups_for_type("Terminal", euro, gbm, "call", cfg.Npaths, rep);
        run_terminal_groups_for_type("Terminal", euro, gbm, "put",  cfg.Npaths, rep);
        rep.finish();
    }

    // ================= Convergence Analyzer (STREAMING) =================
    {
        const std::vector<std::size_t> Ms = { 4, 8, 16, 32, 64, 128, 256};
        const std::size_t NpathsConv = 1000000;

        ConvergenceAnalyzer<double, EuropeanOption<double>, Euler<double>>    ca_eu(euro, euler);
        ConvergenceAnalyzer<double, EuropeanOption<double>, Milstein<double>> ca_mi(euro, milstein);
        ConvergenceAnalyzer<double, EuropeanOption<double>, GBM<double>>      ca_gb(euro, gbm);

        // STRONG / WEAK — Euler
        stream_strong_conv("Euler", ca_eu, NpathsConv, Ms);
        stream_weak_conv("Euler / call", ca_eu, "call", NpathsConv, Ms);

        // STRONG / WEAK — Milstein
        stream_strong_conv("Milstein", ca_mi, NpathsConv, Ms);
        stream_weak_conv("Milstein / call", ca_mi, "call", NpathsConv, Ms);

        // MC statistical (any stepper)
        std::vector<std::size_t> Ns = { 500, 1000, 2000, 5000, 10000, 20000, 50000 };
        stream_mc_stat_conv("GBM (full-path)",      ca_gb, "call", Ns, /*AV*/false, /*CV*/false);
        stream_mc_stat_conv("Euler (full-path)",    ca_eu, "call", Ns, /*AV*/false, /*CV*/false);
        stream_mc_stat_conv("Milstein (full-path)", ca_mi, "call", Ns, /*AV*/false, /*CV*/false);
    }

        // ===== PATH DUMPS FOR PYTHON PLOTS =====
    // Choose a manageable Npaths for file size when plotting (e.g., 200–2000).
    // ====== EXPORT SAMPLE PATHS ======
    try {
        int Npaths = 100;
        int Nsteps = 252; // e.g., 1 year of daily steps
        // JUST TO SHOW HOW MC WORKS
        dump_mc_paths_csv("paths_gbm_call.csv", euro, gbm, "call", Npaths, Nsteps);
        dump_mc_paths_csv("paths_euler_call.csv", euro, euler, "call", Npaths, Nsteps);
        dump_mc_paths_csv("paths_milstein_call.csv", euro, milstein, "call", Npaths, Nsteps);

        std::cout << "\n[paths] Exported " << Npaths << " paths with "
                  << Nsteps << " steps each to CSV files.\n";
    } catch (const std::exception& ex) {
        std::cerr << "[paths] ERROR: " << ex.what() << "\n";
    }

        // ====== EXPORT mean(S_T) vs N FOR 4 VR SETTINGS ======
        // ====== EXPORT price vs N for 10 sequences and 4 VR configurations ======
    try {
    const std::size_t Nmax = 5000;   // x-axis upper bound
    const int numSeq       = 10;      // number of independent sequences
    const std::size_t Nmin = 2;       // avoid N=1 (your engine requires >=2)
    const std::size_t stride = 2;     // increase to 2/5/10 to speed up further

    export_price_sequences_terminal_gbm("PriceSeq_GBMterm_call", euro, "call",
                                        Nmax, numSeq, Nmin, stride);

    std::cout << "[csv] Wrote PriceSeq_GBMterm_call_{none,av,cv,both}.csv\n";
} catch (const std::exception& ex) {
    std::cerr << "[csv] export error (GBM terminal): " << ex.what() << "\n";
}





    return 0;
}
