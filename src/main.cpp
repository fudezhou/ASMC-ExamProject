#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <vector>
#include <sstream>
#include <algorithm>
#include "../include/utilities.hpp"
#include "../include/pde.hpp"

// === NEW: pretty-printer for statistical convergence
static void print_stat_conv(const char* tag, const MCStatConvergence& R)
{
    std::cout << "\n-- MC statistical convergence ("
              << tag << ")  [expect slope ≈ -0.5] --\n";
    std::cout << std::left << std::setw(12) << "Paths"
              << std::setw(18) << "Price"
              << std::setw(18) << "StdErr"
              << std::setw(18) << "CI95 width" << "\n";

    for (std::size_t i=0;i<R.Ns.size();++i) {
        const double ci95 = 2.0 * 1.96 * R.stdErrs[i];
        std::cout << std::left << std::setw(12) << R.Ns[i]
                  << std::setw(18) << std::fixed << std::setprecision(6) << R.prices[i]
                  << std::setw(18) << R.stdErrs[i]
                  << std::setw(18) << ci95 << "\n";
    }
    std::cout << "Estimated order (slope log SE vs log N) ≈ "
              << std::fixed << std::setprecision(4) << R.slope_loglog << "\n";
}
// === NEW: run the statistical convergence for ANY stepper using pricePath
template <typename Step>
void run_statistical_convergence_any_stepper(const char* name,
                                             EuropeanOption<double>& opt,
                                             const Step& stepper,
                                             const std::string& optType,
                                             bool AV, bool CV)
{
    MonteCarlo<double, EuropeanOption<double>, Step> mc(opt, stepper);
    std::vector<std::size_t> Ns = { 500, 1000, 2000, 5000, 10000, 20000, 50000 };

    auto R = mc.statisticalConvergenceByPaths(optType, Ns, AV, CV);

    std::ostringstream tag;
    tag << name << " / " << optType << " / AV=" << (AV?"on":"off") << " / CV=" << (CV?"on":"off");
    print_stat_conv(tag.str().c_str(), R);
}


static double slope_loglog(const std::vector<double>& h,
                           const std::vector<double>& e)
{
    const int n = static_cast<int>(h.size());
    double Sx=0, Sy=0, Sxx=0, Sxy=0;
    for (int i=0;i<n;++i) {
        const double x = std::log(h[i]);
        const double y = std::log(e[i]);
        Sx+=x; Sy+=y; Sxx+=x*x; Sxy+=x*y;
    }
    return (n*Sxy - Sx*Sy) / (n*Sxx - Sx*Sx); // estimated order
}

template <typename Step>
void run_convergence(const char* name,
                     EuropeanOption<double>& opt,
                     const Step& stepper,
                     const std::string& optType,
                     std::size_t Npaths)
{
    MonteCarlo<double, EuropeanOption<double>, Step> mc(opt, stepper);
    mc.useAntitheticVariates(false);
    mc.useControlVariates(false);

    std::vector<std::size_t> tSteps = { 4, 8, 16, 32, 64, 128, 256 , 512 , 1024 };  // time steps
    std::vector<double> h, estrL1, eweak;

    std::cout << "\n-- Convergence (NO AV/CV): " << name << " / " << optType << " --\n";
    std::cout << std::left << std::setw(8)  << "M"
              << std::setw(12) << "dt"
              << std::setw(18) << "Strong L1"
              << std::setw(18) << "Weak bias"
              << std::setw(18) << "SE(bias)" << "\n";

    for (auto step : tSteps) {
        const double dt = opt.getMaturityTime() / double(step);
        auto s = mc.strongErrorOnState(Npaths, step);
        auto w = mc.weakErrorOnPayoff(optType, Npaths, step);

        std::cout << std::left << std::setw(8)  << step
                  << std::setw(12) << std::fixed << std::setprecision(6) << dt
                  << std::setw(18) << s.meanAbs
                  << std::setw(18) << w.bias
                  << std::setw(18) << w.seBias << "\n";

        h.push_back(dt);
        estrL1.push_back(s.meanAbs);
        eweak.push_back(w.bias);
    }

    const double p_strong = slope_loglog(h, estrL1);
    const double p_weak   = slope_loglog(h, eweak);
    std::cout << "Estimated STRONG order ≈ " << p_strong
              << "   |   Estimated WEAK order ≈ " << p_weak << "\n";
}


// ---------- Lightweight table reporter ----------
struct Reporter {
    bool header_printed = false;

    void print_header_once() {
        if (header_printed) return;
        header_printed = true;

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
                  << setw(18)  << "Time (ms)"
                  << "\n";
        std::cout << std::string(118, '-') << "\n";
    }

    void print_row(const char* method,
                   const std::string& optType,
                   bool av, bool cv,
                   std::size_t Npaths, std::size_t Nsteps,
                   double price, double ci95,
                   double ms)
    {
        using std::setw; using std::left; using std::right; using std::fixed; using std::setprecision;
        if (!header_printed) print_header_once();

        std::cout << fixed << setprecision(6)
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
                  << setw(18) << ms
                  << "\n";
    }

    void finish() {
        std::cout << std::string(118, '=') << "\n";
    }
};

// ---------- Scenario banner ----------
inline void print_scenario_banner(const std::string& tag,
                                  double S0, double K, double r, double sigma, double T,
                                  double bs_call, double bs_put)
{
    std::ostringstream lineParams;
    lineParams << std::fixed
               << "S0=" << std::setprecision(6) << S0 << "  "
               << "K=" << std::setprecision(6) << K << "  "
               << "S0/K=" << std::setprecision(6) << (S0 / K) << "  "
               << "r=" << std::setprecision(6) << r << "  "
               << "sigma=" << std::setprecision(6) << sigma << "  "
               << "T=" << std::setprecision(6) << T;

    std::ostringstream lineTitle;
    lineTitle << " SCENARIO: " << tag << "  ";

    std::ostringstream lineBS;
    lineBS << " BS Call=" << std::fixed << std::setprecision(6) << bs_call
           << "    BS Put=" << bs_put;

    const std::size_t w = std::max({ lineTitle.str().size(), lineParams.str().size(), lineBS.str().size() });
    const std::string eq(w, '=');
    const std::string dash(w, '-');

    std::cout << "\n" << eq << "\n";
    std::cout << lineTitle.str() << "\n";
    std::cout << dash << "\n";
    std::cout << lineParams.str() << "\n";
    std::cout << dash << "\n";
    std::cout << lineBS.str() << "\n";
    std::cout << eq << "\n";
}

// ---------- Helpers to run a single pricing row ----------
template <typename Step>
void run_fullpath_row(const char* methodName,
                      EuropeanOption<double>& opt,
                      const Step& stepper,
                      const std::string& optType,
                      bool useAV, bool useCV,
                      std::size_t Npaths, std::size_t Nsteps,
                      Reporter& rep)
{
    MonteCarlo<double, EuropeanOption<double>, Step> mc(opt, stepper);
    mc.useAntitheticVariates(useAV);
    mc.useControlVariates(useCV);

    const auto t0 = std::chrono::steady_clock::now();
    auto r = mc.pricePath(optType, Npaths);
    const auto t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double ci95 = 1.96 * r.stdDev;

    rep.print_row(methodName, optType, useAV, useCV, Npaths, Nsteps, r.price, ci95, ms);
}

template <typename Step>
void run_terminal_row(const char* methodName,
                      EuropeanOption<double>& opt,
                      const Step& stepper,
                      const std::string& optType,
                      bool useAV, bool useCV,
                      std::size_t Npaths,
                      Reporter& rep)
{
    MonteCarlo<double, EuropeanOption<double>, Step> mc(opt, stepper);
    mc.useAntitheticVariates(useAV);
    mc.useControlVariates(useCV);

    const auto t0 = std::chrono::steady_clock::now();
    auto r = mc.priceTerminal(optType, Npaths);
    const auto t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double ci95 = 1.96 * r.stdDev;

    rep.print_row(methodName, optType, useAV, useCV, Npaths, /*Steps*/ 1, r.price, ci95, ms);
}

// ---------- Print groups in the requested order ----------
template <typename Step>
void run_fullpath_groups_for_type(const char* methodName,
                                  EuropeanOption<double>& opt,
                                  const Step& stepper,
                                  const std::string& optType,
                                  std::size_t Npaths,
                                  Reporter& rep)
{
    const std::size_t Nsteps = std::max<std::size_t>(
        1, static_cast<std::size_t>(TRADING_DAYS_PER_YEAR * opt.getMaturityTime())
    );

    // 1) No AV/CV
    run_fullpath_row(methodName, opt, stepper, optType, /*AV*/false, /*CV*/false, Npaths, Nsteps, rep);
    // 2) AV only
    run_fullpath_row(methodName, opt, stepper, optType, /*AV*/true,  /*CV*/false, Npaths, Nsteps, rep);
    // 3) CV only
    run_fullpath_row(methodName, opt, stepper, optType, /*AV*/false, /*CV*/true,  Npaths, Nsteps, rep);
    // 4) AV + CV
    run_fullpath_row(methodName, opt, stepper, optType, /*AV*/true,  /*CV*/true,  Npaths, Nsteps, rep);
}

template <typename Step>
void run_terminal_groups_for_type(const char* methodName,
                                  EuropeanOption<double>& opt,
                                  const Step& stepper,
                                  const std::string& optType,
                                  std::size_t Npaths,
                                  Reporter& rep)
{
    // 1) No AV/CV
    run_terminal_row(methodName, opt, stepper, optType, /*AV*/false, /*CV*/false, Npaths, rep);
    // 2) AV only
    run_terminal_row(methodName, opt, stepper, optType, /*AV*/true,  /*CV*/false, Npaths, rep);
    // 3) CV only
    run_terminal_row(methodName, opt, stepper, optType, /*AV*/false, /*CV*/true,  Npaths, rep);
    // 4) AV + CV
    run_terminal_row(methodName, opt, stepper, optType, /*AV*/true,  /*CV*/true,  Npaths, rep);
}

// ---------- MAIN (ATM only; ordered output) ----------
int main()
{
    struct {
        double S0 = 100.0;
        double r  = 0.05;
        double v  = 0.20;
        double T  = 1.0;
        std::size_t Npaths = 10000; // tweak freely
    } cfg;

    // ATM: K = S0
    EuropeanOption<double> euro(cfg.S0, /*K=*/cfg.S0, cfg.r, cfg.v, cfg.T);

    GBM<double>      gbm;
    Euler<double>    euler;
    Milstein<double> milstein;

    // Black–Scholes references for banner
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

        // ---- All CALL rows grouped by VR setting (No/AV/CV/AV+CV) across all methods
        rep.print_header_once();
        run_fullpath_groups_for_type("GBM",      euro, gbm,      "call", cfg.Npaths, rep);
        run_fullpath_groups_for_type("Euler",    euro, euler,    "call", cfg.Npaths, rep);
        run_fullpath_groups_for_type("Milstein", euro, milstein, "call", cfg.Npaths, rep);

        // ---- All PUT rows grouped by VR setting (No/AV/CV/AV+CV) across all methods
        run_fullpath_groups_for_type("GBM",      euro, gbm,      "put", cfg.Npaths, rep);
        run_fullpath_groups_for_type("Euler",    euro, euler,    "put", cfg.Npaths, rep);
        run_fullpath_groups_for_type("Milstein", euro, milstein, "put", cfg.Npaths, rep);

        rep.finish();
    }

    // ================= Terminal (single-step, exact GBM) =================
    std::cout << "\n********** TERMINAL (single-step, exact GBM) **********\n";
    {
        Reporter rep;
        rep.print_header_once();

        // Only GBM stepper is meaningful here
        run_terminal_groups_for_type("Terminal", euro, gbm, "call", cfg.Npaths, rep);
        run_terminal_groups_for_type("Terminal", euro, gbm, "put",  cfg.Npaths, rep);

        rep.finish();
    }

    {
    const std::size_t NpathsConv = 100000; // tune for stability/speed
    run_convergence("Euler",    euro, euler,    "call", NpathsConv);
    run_convergence("Milstein", euro, milstein, "call", NpathsConv);
    run_convergence("GBM", euro, gbm, "call", NpathsConv);
    }

    {
    // Demonstrate on all steppers; AV/CV off to show the raw -1/2 slope.
    run_statistical_convergence_any_stepper("GBM (full-path)",      euro, gbm,      "call", /*AV*/false, /*CV*/false);
    run_statistical_convergence_any_stepper("Euler (full-path)",    euro, euler,    "call", /*AV*/false, /*CV*/false);
    run_statistical_convergence_any_stepper("Milstein (full-path)", euro, milstein, "call", /*AV*/false, /*CV*/false);

    // Optional: show that VR changes the constant, not the slope
    run_statistical_convergence_any_stepper("Euler (full-path)",    euro, euler,    "call", /*AV*/true,  /*CV*/false);
    run_statistical_convergence_any_stepper("Euler (full-path)",    euro, euler,    "call", /*AV*/false, /*CV*/true);
    }

    //     // ================= PDE (Crank–Nicolson) =================
    // std::cout << "\n********** PDE (Crank–Nicolson) **********\n";
    // {
    //     auto t0 = std::chrono::steady_clock::now();
    //     PDECNPricer<double> pde(euro);
    //     const double Smax = 5.0 * euro.getSpotPrice(); // good default; tweak freely
    //     auto cn_call = pde.price("call", /*M*/600, /*N*/600, Smax, /*rannacher*/true);
    //     auto cn_put  = pde.price("put",  /*M*/600, /*N*/600, Smax, /*rannacher*/true);
    //     auto t1 = std::chrono::steady_clock::now();
    //     const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    //     std::cout << std::fixed << std::setprecision(6);
    //     std::cout << "CN Call: " << cn_call.price
    //               << "   (Δ=" << cn_call.delta
    //               << ", Γ=" << cn_call.gamma
    //               << ", Θ=" << cn_call.theta
    //               << ")   [M=" << cn_call.M << ", N=" << cn_call.N
    //               << ", Smax=" << cn_call.Smax << ", " << ms << " ms]\n";

    //     std::cout << "CN  Put: " << cn_put.price
    //               << "   (Δ=" << cn_put.delta
    //               << ", Γ=" << cn_put.gamma
    //               << ", Θ=" << cn_put.theta
    //               << ")   [M=" << cn_put.M << ", N=" << cn_put.N
    //               << ", Smax=" << cn_put.Smax << "]\n";
    // }


    return 0;
}


