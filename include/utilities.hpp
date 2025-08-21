#ifndef UTILITIES_HPP
#define UTILITIES_HPP

/*
    In this .hpp file only declarations are made.
    The implementation is done in the corresponding .cpp file.
*/

#include <iostream>
#include <chrono>   // using C++ dates from <chrono> STL library
#include <cmath>    // for CDF: Cumulative Distribution Function
#include <numbers>  // for mathematical constants
#include <random>
#include <functional>
#include <algorithm>
#include <stdexcept>


inline constexpr double INV_SQRT_2PI = 0.39894228040143267793994605993438; // 1 / sqrt(2 * pi)
inline constexpr double INV_SQRT_2   = 0.70710678118654752440084436210485; // 1 / sqrt(2)
inline constexpr double TRADING_DAYS_PER_YEAR = 252.0;

// Standard normal PDF
inline double norm_pdf(double x) 
{
    return INV_SQRT_2PI * std::exp(-0.5 * x * x);
}

// Standard normal CDF
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x * INV_SQRT_2);
}

inline std::mt19937_64& global_rng() {
    static std::mt19937_64 rng(42); // fixed seed for reproducibility
    return rng;
}

inline double standard_normal_sample() {
    static thread_local std::normal_distribution<double> dist(0.0, 1.0);
    return dist(global_rng());
}

template <typename T>
struct MCReturn
{
    T price  = 0.0; // option price
    T var    = 0.0; // variance of the price
    T stdDev = 0.0; // standard deviation of the price
};

// template <typename T>
// inline T gbm_ST (T S0, T r, T sigma, T maturity, T Z)
// {
//     return S0 * std::exp((r - 0.5 * sigma * sigma) * maturity + sigma * Z * std::sqrt(maturity));
// }

// inline double standard_normal_sample() {
//     static thread_local std::mt19937_64 rng{std::random_device{}()};  // engine (64-bit Mersenne Twister)
//     static thread_local std::normal_distribution<double> dist(0.0, 1.0); // mean=0, stddev=1
//     return dist(rng);
// }

// === NEW: small helper for log–log slope (kept header-only) ===
inline double slope_loglog_xy(const std::vector<double>& x,
                              const std::vector<double>& y)
{
    const int n = static_cast<int>(x.size());
    double Sx=0, Sy=0, Sxx=0, Sxy=0;
    for (int i=0;i<n;++i) {
        const double lx = std::log(x[i]);
        const double ly = std::log(y[i]);
        Sx+=lx; Sy+=ly; Sxx+=lx*lx; Sxy+=lx*ly;
    }
    return (n*Sxy - Sx*Sy) / (n*Sxx - Sx*Sx);
}

template <typename T>
class Option { // base class for all options
public: 
    Option(const std::string& name)
        : _name(name) {}
    Option(const std::string& name,
           T spotPrice,
           T strikePrice,
           T riskFreeRate,
           T volatility,
           T maturityTime)
        : _name(name)
        , _spotPrice(spotPrice)
        , _strikePrice(strikePrice)
        , _riskFreeRate(riskFreeRate)
        , _volatility(volatility)
        , _maturityTime(maturityTime) {}

    // setters declarations
    void setName(const std::string& name);
    void setType(const std::string& optionType);
    void setSpotPrice(T spotPrice);
    void setStrikePrice(T strikePrice);
    void setRiskFreeRate(T riskFreeRate);
    void setVolatility(T volatility);
    void setMaturityTime(T maturityTime);

    // getters declarations
    const std::string& getName() const;
    const std::string& getType() const;
    const T& getSpotPrice() const;
    const T& getStrikePrice() const;
    const T& getRiskFreeRate() const;
    const T& getVolatility() const;
    const T& getMaturityTime() const;

    // info
    virtual void info() const;

    // declaring a pure abstract method for payoff calculation
    // this method must be implemented by all derived option classes
    // this is polymorphic behavior
    // I want payOff to return a function type
    virtual const std::function<T(T)> payOff() const = 0;

    virtual ~Option() = default;

private:
    std::string _name; // name of the option (e.g., "European" , "Asian", "American")
    std::string _optionType; // type of the option (e.g., "Call", "Put")    
    T           _spotPrice;
    T           _strikePrice;
    T           _riskFreeRate;
    T           _volatility;
    T           _maturityTime; // time to maturity in years
};

// inheritance of European class
template <typename T>
class EuropeanOption : public Option<T> {
public:
    EuropeanOption() : Option<T>("European") {}
    EuropeanOption(T spotPrice,
                   T strikePrice,
                   T riskFreeRate,
                   T volatility,
                   T maturityTime)
        : Option<T>("European", spotPrice, strikePrice, riskFreeRate, volatility, maturityTime) {}

    // Closed form solution Black-Scholes
    T closedForm(const std::string& optionType);

    // Override info function
    void info() const override;

    virtual const std::function<T(T)> payOff() const override; 
//     unsigned long _numPaths = 100000; // old First and Second naive implementations
//     std::vector<std::vector<T>> _pricePaths;
// private:
//     std::size_t _stride; 
//     std::vector<T> _pricePaths_flat; // THIS IS A FAR BETTER DESIGN CHOICE!!!

};

// Stepper (Abstract Base Class)
// A stepper or path generator for simulating price paths
template <typename T>
class Stepper {
public:
    Stepper() = default;
    virtual ~Stepper() = default;

    // Advance the process by one time step
    virtual T advance(T S_t, T r, T sigma, T dt) const = 0;

    virtual T advanceWithZ(T S_t, T r, T sigma, T dt, T /*Z*/) const
    {
        return advance(S_t, r, sigma, dt); // default implementation
    }

    // Generate a terminal random draw (Only for GBM)
    virtual T terminalDraw(T S0, T r, T sigma, T maturity) const
    {
        throw std::logic_error("terminalDraw() not available for this stepper");
    }

    virtual T terminalDrawWithZ(T S0, T r, T sigma, T maturity, T Z) const
    {
        throw std::logic_error("terminalDrawWithZ() not available for this stepper");
    }

    // info
    virtual const char* getName() const = 0;
};

template <typename T>
class GBM : public Stepper<T> {
public:
    GBM() = default;
    ~GBM() override = default;

    T advance(T S_t, T r, T sigma, T dt) const override;
    T advanceWithZ(T S_t, T r, T sigma, T dt, T Z) const override;
    T terminalDraw(T S0, T r, T sigma, T maturity) const override;
    T terminalDrawWithZ(T S0, T r, T sigma, T maturity, T Z) const override;

    virtual const char *getName() const override { return "GBM"; };
};

template <typename T>
class Euler : public Stepper<T> {
public:
    Euler() = default;
    ~Euler() override = default;

    T advance(T S_t, T r, T sigma, T dt) const override;
    T advanceWithZ(T S_t, T r, T sigma, T dt, T Z) const override;

    virtual const char *getName() const override { return "Euler"; }
};

template <typename T>
class Milstein : public Stepper<T> {
public:
    Milstein() = default;
    ~Milstein() override = default;

    T advance(T S_t, T r, T sigma, T dt) const override;
    T advanceWithZ(T S_t, T r, T sigma, T dt, T Z) const override;

    virtual const char *getName() const override { return "Milstein"; }
};

// ---------- Convergence statistics ----------
struct StrongStats {
    double meanAbs;  // E[ |S*_T - S^h_T| ]
    double rms;      // sqrt( E[ (S*_T - S^h_T)^2 ] )
    double seAbs;    // std error of meanAbs
};

struct WeakStats {
    double bias;       // | E[g(S^h_T)] - E[g(S*_T)] |
    double seBias;     // std error of the mean of D = g(S^h_T) - g(S*_T)
    double meanApprox; // E[g(S^h_T)] estimate
    double meanExact;  // E[g(S*_T)] estimate
};

// === NEW: MC statistical convergence result ===
struct MCStatConvergence {
    std::vector<std::size_t> Ns;   // path counts
    std::vector<double> prices;    // MC prices (per N)
    std::vector<double> stdErrs;   // reported standard errors (per N)
    double slope_loglog = 0.0;     // slope of log(SE) vs log(N), expect ≈ -0.5
};

// ========= Convergence result containers =========
struct StrongConvergenceResult {
    std::vector<std::size_t> M;   // time steps per year
    std::vector<double> dt;       // timestep size
    std::vector<double> meanAbs;  // E|S_T^exact - S_T^num|
    std::vector<double> rms;      // sqrt(E[(S_T^exact - S_T^num)^2])
    double slope_meanAbs = 0.0;   // slope of log(meanAbs) vs log(dt)
    double slope_rms     = 0.0;   // slope of log(rms) vs log(dt)
};

struct WeakConvergenceResult {
    std::vector<std::size_t> M;
    std::vector<double> dt;
    std::vector<double> bias;     // |E[g(S_T^num)] - E[g(S_T^exact)]|
    std::vector<double> seBias;   // SE of bias estimator (for CI)
    double slope_bias = 0.0;      // slope of log(bias) vs log(dt)
};

// ========= Convergence Analyzer =========
template <typename T, typename Opt, typename Step>
class ConvergenceAnalyzer {
public:
    ConvergenceAnalyzer(const Opt& opt, const Step& stepper)
    : _opt(opt), _stepper(stepper) {}

    // MC statistical: uses MonteCarlo::statisticalConvergenceByPaths
    MCStatConvergence mcConvergence(const std::string& optionType,
                                    const std::vector<std::size_t>& Ns,
                                    bool useAV=false, bool useCV=false) const;

    // Strong: uses MonteCarlo::strongErrorOnState over an M-grid
    StrongConvergenceResult strongConvergence(std::size_t numPaths,
                                              const std::vector<std::size_t>& Ms) const;

    // Weak: uses MonteCarlo::weakErrorOnPayoff over an M-grid
    WeakConvergenceResult weakConvergence(const std::string& optionType,
                                          std::size_t numPaths,
                                          const std::vector<std::size_t>& Ms) const;

private:
    Opt       _opt;
    Step      _stepper;
};

// Utility: slope (simple OLS) of log(y) vs log(x)
double slope_loglog_xy(const std::vector<double>& x,
                       const std::vector<double>& y);


template <typename T, typename Opt, typename Step>
class MonteCarlo {
public:
    MonteCarlo(Opt& option, const Step& stepper)
        : _option(option), _stepper(stepper) {}

    ~MonteCarlo() = default;

    void setNumPaths(std::size_t numPaths);
    void setNumSteps(std::size_t numSteps);
    void setStride(std::size_t stride);
    void useAntitheticVariates(bool useAV);
    void useControlVariates(bool useCV);

    // NEW SIGNATURE: terminal (single-step) Monte Carlo using exact GBM terminal draw.
    MCReturn<T> priceTerminal(const std::string& optionType, std::size_t numPaths);
    // NEW SIGNATURE: full path Monte Carlo using stepper
    MCReturn<T> pricePath(const std::string& optionType, std::size_t numPaths);

    StrongStats strongErrorOnState(std::size_t numPaths, std::size_t numSteps);

    WeakStats weakErrorOnPayoff(const std::string& optionType, std::size_t numPaths, std::size_t numSteps);

    MCStatConvergence statisticalConvergenceByPaths(const std::string& optionType,
                                                const std::vector<std::size_t>& Ns,
                                                bool useAV=false,
                                                bool useCV=false);

private: 
    Opt& _option; // reference to the option
    const Step& _stepper; // reference to the stepper
    std::size_t _numPaths = 100000; // number of paths to generate
    std::size_t _numSteps = 252;
    bool _useAV = false;
    bool _useCV = false;
    std::size_t _stride;
    std::vector<T> _pricePaths_flat;
};



#endif // UTILITIES_HPP