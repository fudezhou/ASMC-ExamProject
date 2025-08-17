#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>
#include <chrono>   // using C++ dates from <chrono> STL library
#include <cmath>    // for CDF: Cumulative Distribution Function
#include <numbers>  // for mathematical constants


inline constexpr double INV_SQRT_2PI = 0.39894228040143267793994605993438; // 1 / sqrt(2 * pi)
inline constexpr double INV_SQRT_2   = 0.70710678118654752440084436210485; // 1 / sqrt(2)

// Standard normal PDF
inline double norm_pdf(double x) 
{
    return INV_SQRT_2PI * std::exp(-0.5 * x * x);
}

// Standard normal CDF
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x * INV_SQRT_2);
}

/*
    In this .hpp file only declarations are made.
    The implementation is done in the corresponding .cpp file.
*/

// Forward declaration of the BlackAndScholes class template
// (Since 1973)
template <typename T>
class BlackAndScholes {
public:
    // constructor
    BlackAndScholes(const T&S, 
                    const T& K, 
                    const T& r, 
                    const T& sigma,    
                    const T& t_T)
        : _spotPrice(S)
        , _strikePrice(K)
        , _riskFreeRate(r)
        , _volatility(sigma) 
        , _t_T(t_T)
    {
        std::cout << "===============================================" << std::endl;
        std::cout << "Black & Scholes object created with parameters:" << std::endl;
        info();
        std::cout << "===============================================" << std::endl;
    }

    BlackAndScholes() = default;

    // setters
    void setSpotPrice(const T& S);
    void setStrikePrice(const T& K);
    void setRiskFreeRate(const T& r);
    void setTimeToExpiration(const T& t_T);
    void setVolatility(const T& sigma);

    // getters
    const T&    getSpotPrice() const;
    const T&    getStrikePrice() const;
    const T&    getRiskFreeRate() const;
    const T&    getTimeToExpiration() const;
    const T&    getVolatility() const;

    virtual void info() const;

    // the following price() method takes a string argument
    // that specifies the type of option either "call" or "put"
    virtual double price(const std::string& optionType) const;

    // destructor
    virtual ~BlackAndScholes() = default;

private:
    // member variables
    T     _spotPrice;
    T     _strikePrice;
    T     _riskFreeRate;
    T     _volatility;
    T     _t_T;
};

#endif // UTILITIES_HPP