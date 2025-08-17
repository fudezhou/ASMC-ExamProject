#include <../include/utilities.hpp>

template <typename T>
void BlackAndScholes<T>::setSpotPrice(const T& S) 
{
    _spotPrice = S;
}

template <typename T>
void BlackAndScholes<T>::setStrikePrice(const T& K) 
{
    _strikePrice = K;
}

template <typename T>
void BlackAndScholes<T>::setRiskFreeRate(const T& r) 
{
    _riskFreeRate = r;
}

template <typename T>
void BlackAndScholes<T>::setTimeToExpiration(const T& t_T) 
{
    _t_T = t_T;
}

template <typename T>
void BlackAndScholes<T>::setVolatility(const T& sigma) 
{
    _volatility = sigma;
}

template <typename T>
const T& BlackAndScholes<T>::getSpotPrice() const 
{
    return _spotPrice;
}

template <typename T>
const T& BlackAndScholes<T>::getStrikePrice() const 
{
    return _strikePrice;
}


template <typename T>
const T& BlackAndScholes<T>::getRiskFreeRate() const 
{
    return _riskFreeRate;
}

template <typename T>
const T& BlackAndScholes<T>::getTimeToExpiration() const 
{
    return _t_T;
}

template <typename T>
const T& BlackAndScholes<T>::getVolatility() const 
{
    return _volatility;
}

template <typename T>
void BlackAndScholes<T>::info() const 
{
    std::cout << "Spot Price         : " << _spotPrice << std::endl;
    std::cout << "Strike Price       : " << _strikePrice << std::endl;
    std::cout << "Risk-Free Rate     : " << _riskFreeRate << std::endl;
    std::cout << "Volatility         : " << _volatility << std::endl;
    std::cout << "Time to Expiration : " << _t_T << std::endl;
}

// template <typename T>
// double BlackAndScholes<T>::price(const std::string& optionType) const 
// {
//     // check if the string is a valid one
//     if (optionType != "call" && optionType != "put") {
//         throw std::invalid_argument("Invalid option type");
//     }
//     // Implement the Black-Scholes formula for call and put options
//     double d1 = (std::log(_spotPrice / _strikePrice) + (_riskFreeRate + 0.5 * _volatility * _volatility) * _t_T) / (_volatility * std::sqrt(_t_T));
//     double d2 = d1 - _volatility * std::sqrt(_t_T);

//     if (optionType == "call") {
//         return _spotPrice * norm_cdf(d1) - _strikePrice * std::exp(-_riskFreeRate * _t_T) * norm_cdf(d2);
//     } else {
//         return _strikePrice * std::exp(-_riskFreeRate * _t_T) * norm_cdf(-d2) - _spotPrice * norm_cdf(-d1);
//     }
// }

template <typename T>
double BlackAndScholes<T>::price(const std::string& optionType) const 
{
    // check if the string is a valid one
    if (optionType != "call" && optionType != "put") {
        throw std::invalid_argument("Invalid option type");
    }

    // map to payoff sign φ: call=+1, put=-1
    const double phi = (optionType == "call") ? +1.0 : -1.0;

    return phi * (
           _spotPrice * norm_cdf(
               phi * (
                   (std::log(_spotPrice / _strikePrice) +
                   (_riskFreeRate + 0.5 * _volatility * _volatility) * _t_T)
                   / (_volatility * std::sqrt(_t_T))
               )
           )
           - _strikePrice * std::exp(-_riskFreeRate * _t_T) *
           norm_cdf(
               phi * (
                   (std::log(_spotPrice / _strikePrice) +
                   (_riskFreeRate - 0.5 * _volatility * _volatility) * _t_T)
                   / (_volatility * std::sqrt(_t_T))
               )
           )
    );
};

template class BlackAndScholes<double>;