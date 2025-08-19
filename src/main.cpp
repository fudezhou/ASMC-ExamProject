// #include <iostream>
// // #include <iomanip> // for std::setprecision
// #include "../include/utilities.hpp"

// int main() {

//     EuropeanOption<double> euroOption(100.0, 100.0, 0.05, 0.2, 1.0);
//     euroOption.info();
//     euroOption.closedForm("call");
//     euroOption.closedForm("put");

//     std::cout << "\nEuropean Option Call Price: " << euroOption.closedForm("call") << std::endl;
//     std::cout << "European Option Put Price: " << euroOption.closedForm("put") << std::endl;

//     // testing payoff function
//     euroOption.setType("call");
//     auto callPayoff = euroOption.payOff();
//     euroOption.setType("put");
//     auto putPayoff = euroOption.payOff();

//     // Montecarlo full path with GBM
//     GBM<double> gbm;
//     MonteCarlo<double, EuropeanOption<double>, GBM<double>> mc(euroOption, gbm);

//     auto result = mc.pricePath("call", 10000);
//     std::cout << "Monte Carlo with " << gbm.getName() << " for Call Option Price: " << result.price << " ± " << 1.96 * result.stdDev << " (95% CI)\n";

//     result = mc.pricePath("put", 10000);
//     std::cout << "Monte Carlo with " << gbm.getName() << " for Put Option Price: " << result.price << " ± " << 1.96 * result.stdDev << " (95% CI)\n";

//     // Monte Carlo with Euler
//     Euler<double> euler;
//     MonteCarlo<double, EuropeanOption<double>, Euler<double>> mcEuler(euroOption, euler);

//     result = mcEuler.pricePath("call", 10000);
//     std::cout << "Monte Carlo with " << euler.getName() << " for Call Option Price: " << result.price << " ± " << 1.96 * result.stdDev << " (95% CI)\n";

//     result = mcEuler.pricePath("put", 10000);
//     std::cout << "Monte Carlo with " << euler.getName() << " for Put Option Price: " << result.price << " ± " << 1.96 * result.stdDev << " (95% CI)\n";

//     // Monte Carlo with Milstein
//     Milstein<double> milstein;
//     MonteCarlo<double, EuropeanOption<double>, Milstein<double>> mcMilstein(euroOption, milstein);

//     result = mcMilstein.pricePath("call", 10000);
//     std::cout << "Monte Carlo with " << milstein.getName() << " for Call Option Price: " << result.price << " ± " << 1.96 * result.stdDev << " (95% CI)\n";

//     result = mcMilstein.pricePath("put", 10000);
//     std::cout << "Monte Carlo with " << milstein.getName() << " for Put Option Price: " << result.price << " ± " << 1.96 * result.stdDev << " (95% CI)\n";

//     return 0;
// }

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <vector>
#include <sstream>
#include <algorithm>
#include "../include/utilities.hpp"

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

// ---------- Scenario banner (exact style you requested) ----------
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

    // Choose width = max(len of title line, params line, BS line), then pad borders to that width
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

// ---------- Full-path suite (discretized) ----------
template <typename Step>
void run_suite(const char* methodName,
               EuropeanOption<double>& opt,
               const Step& stepper,
               std::size_t Npaths,
               Reporter& rep)
{
    const std::size_t Nsteps = std::max<std::size_t>(
        1, static_cast<std::size_t>(TRADING_DAYS_PER_YEAR * opt.getMaturityTime())
    );

    auto run = [&](const std::string& optType, bool useAV, bool useCV) {
        MonteCarlo<double, EuropeanOption<double>, Step> mc(opt, stepper);
        mc.useAntitheticVariates(useAV);
        mc.useControlVariates(useCV);

        const auto t0 = std::chrono::steady_clock::now();
        auto r = mc.pricePath(optType, Npaths);
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        const double ci95 = 1.96 * r.stdDev;
        rep.print_row(methodName, optType, useAV, useCV, Npaths, Nsteps, r.price, ci95, ms);
    };

    run("call", false, false);
    run("call", true,  false);
    run("call", false, true);
    run("call", true,  true);

    run("put",  false, false);
    run("put",  true,  false);
    run("put",  false, true);
    run("put",  true,  true);
}

// ---------- Terminal (single-step, exact GBM) ----------
template <typename Step>
void run_terminal_suite(const char* methodName,
                        EuropeanOption<double>& opt,
                        const Step& stepper,
                        std::size_t Npaths,
                        Reporter& rep)
{
    auto run = [&](const std::string& optType, bool useAV, bool useCV) {
        MonteCarlo<double, EuropeanOption<double>, Step> mc(opt, stepper);
        mc.useAntitheticVariates(useAV);
        mc.useControlVariates(useCV);

        const auto t0 = std::chrono::steady_clock::now();
        auto r = mc.priceTerminal(optType, Npaths);
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        const double ci95 = 1.96 * r.stdDev;
        rep.print_row(methodName, optType, useAV, useCV, Npaths, /*Steps*/ 1, r.price, ci95, ms);
    };

    run("call", false, false);
    run("call", true,  false);
    run("call", false, true);
    run("call", true,  true);

    run("put",  false, false);
    run("put",  true,  false);
    run("put",  false, true);
    run("put",  true,  true);
}

// ---------- MAIN ----------
int main()
{
    struct {
        double S0 = 100.0;
        double r  = 0.05;
        double v  = 0.20;
        double T  = 1.0;
        std::size_t Npaths = 10000; // tweak freely
    } cfg;

    EuropeanOption<double> euro(cfg.S0, /*K=*/cfg.S0, cfg.r, cfg.v, cfg.T);
    GBM<double>      gbm;
    Euler<double>    euler;
    Milstein<double> milstein;

    struct Scenario { const char* tag; double Kmult; };
    const std::vector<Scenario> scenarios = {
        {"ITM (for calls) / OTM (for puts)", 0.90}, // K = 0.9 * S0
        {"ATM",                               1.00}, // K = 1.0 * S0
        {"OTM (for calls) / ITM (for puts)", 1.10}  // K = 1.1 * S0
    };

    for (const auto& sc : scenarios) {
        euro.setStrikePrice(sc.Kmult * cfg.S0);

        // BS refs (just for the banner)
        const double bs_call = euro.closedForm("call");
        const double bs_put  = euro.closedForm("put");

        // Scenario banner (exact layout)
        print_scenario_banner(sc.tag,
                              euro.getSpotPrice(), euro.getStrikePrice(),
                              euro.getRiskFreeRate(), euro.getVolatility(), euro.getMaturityTime(),
                              bs_call, bs_put);

        // FULL-PATH
        std::cout << "\n********** FULL-PATH (discretized) MONTE CARLO **********\n";
        Reporter rep_path;
        rep_path.print_header_once();
        run_suite("GBM",      euro, gbm,      cfg.Npaths, rep_path);
        run_suite("Euler",    euro, euler,    cfg.Npaths, rep_path);
        run_suite("Milstein", euro, milstein, cfg.Npaths, rep_path);
        rep_path.finish();

        // TERMINAL
        std::cout << "\n********** TERMINAL (single-step, exact GBM) **********\n";
        Reporter rep_term;
        rep_term.print_header_once();
        run_terminal_suite("Terminal", euro, gbm, cfg.Npaths, rep_term);
        rep_term.finish();
    }

    return 0;
}




