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

template <typename T>
inline T gbm_ST (T S0, T r, T sigma, T maturity, T Z)
{
    return S0 * std::exp((r - 0.5 * sigma * sigma) * maturity + sigma * Z * std::sqrt(maturity));
}

// inline double standard_normal_sample() {
//     static thread_local std::mt19937_64 rng{std::random_device{}()};  // engine (64-bit Mersenne Twister)
//     static thread_local std::normal_distribution<double> dist(0.0, 1.0); // mean=0, stddev=1
//     return dist(rng);
// }

template <typename T>
class Option { // base class for all options
public: 
    Option() = default;
    Option(T spotPrice, 
           T strikePrice, 
           T riskFreeRate, 
           T volatility, 
           T maturityTime)
        : _spotPrice(spotPrice)
        , _strikePrice(strikePrice)
        , _riskFreeRate(riskFreeRate)
        , _volatility(volatility)
        , _maturityTime(maturityTime) {}

    // setters declarations
    void setType(const std::string& optionType);
    void setSpotPrice(T spotPrice);
    void setStrikePrice(T strikePrice);
    void setRiskFreeRate(T riskFreeRate);
    void setVolatility(T volatility);
    void setMaturityTime(T maturityTime);

    // getters declarations
    const std::string& getType() const;
    const T& getSpotPrice() const;
    const T& getStrikePrice() const;
    const T& getRiskFreeRate() const;
    const T& getVolatility() const;
    const T& getMaturityTime() const;

    // info
    virtual void info() const;

    virtual ~Option() = default;

private:
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
    EuropeanOption() = default;
    EuropeanOption(T spotPrice, T strikePrice, T riskFreeRate, T volatility, T maturityTime)
        : Option<T>(spotPrice, strikePrice, riskFreeRate, volatility, maturityTime) {}

    // setter for _numPaths
    void setNumPaths(const unsigned long& numPaths);

    // getter for _numPaths
    unsigned long getNumPaths() const;

    // Closed form solution
    T closedForm(const std::string& optionType);

    // Montecarlo full path with GBM
    const MCReturn<T> naiveMonteCarloGBM(const std::string& optionType, int numPaths);

    const MCReturn<T> fastLowVarMCGBM(const std::string& optionType, int numPaths, bool antithetic = true, bool controlVariates = true);

    // Override info function
    void info() const override;

private: 
//     unsigned long _numPaths = 100000; // old First and Second naive implementations
//     std::vector<std::vector<T>> _pricePaths;
    std::size_t _stride; 
    std::vector<T> _pricePaths_flat; // THIS IS A FAR BETTER DESIGN CHOICE!!!

};


#endif // UTILITIES_HPP