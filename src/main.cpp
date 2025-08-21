#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <vector>
#include <sstream>
#include <algorithm>
#include "../include/utilities.hpp"

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
        std::vector<std::size_t> Ms_prefix(Ms.begin(), Ms.begin() + i + 1);
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
    run_fullpath_row(methodName, opt, stepper, optType, false, false, Npaths, Nsteps, rep);
    run_fullpath_row(methodName, opt, stepper, optType, true,  false, Npaths, Nsteps, rep);
    run_fullpath_row(methodName, opt, stepper, optType, false, true,  Npaths, Nsteps, rep);
    run_fullpath_row(methodName, opt, stepper, optType, true,  true,  Npaths, Nsteps, rep);
}

template <typename Step>
void run_terminal_groups_for_type(const char* methodName, EuropeanOption<double>& opt, const Step& stepper,
                                  const std::string& optType, std::size_t Npaths, Reporter& rep)
{
    run_terminal_row(methodName, opt, stepper, optType, false, false, Npaths, rep);
    run_terminal_row(methodName, opt, stepper, optType, true,  false, Npaths, rep);
    run_terminal_row(methodName, opt, stepper, optType, false, true,  Npaths, rep);
    run_terminal_row(methodName, opt, stepper, optType, true,  true,  Npaths, rep);
}

// ---------- MAIN ----------
int main()
{
    struct {
        double S0 = 100.0;
        double r  = 0.05;
        double v  = 0.20;
        double T  = 1.0;
        std::size_t Npaths = 100000;
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
        const std::size_t NpathsConv = 500000;

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


    return 0;
}


