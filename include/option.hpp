#ifndef OPTION_HPP
#define OPTION_HPP

#include <string>
#include <functional>
#include <iostream>
#include <cmath>
#include "math.hpp"

template <typename T>
class Option {
public:
    Option() = default;
    Option(T S0, T K, T r, T sigma, T Tmat)
    : _spotPrice(S0), _strikePrice(K), _riskFreeRate(r), _volatility(sigma), _maturityTime(Tmat) {}

    // Setters
    void setName(const std::string& name)        { _name = name; }
    void setType(const std::string& optionType)  { _optionType = optionType; }
    void setSpotPrice(T spotPrice)               { _spotPrice = spotPrice; }
    void setStrikePrice(T strikePrice)           { _strikePrice = strikePrice; }
    void setRiskFreeRate(T riskFreeRate)         { _riskFreeRate = riskFreeRate; }
    void setVolatility(T volatility)             { _volatility = volatility; }
    void setMaturityTime(T maturityTime)         { _maturityTime = maturityTime; }

    // Getters
    const std::string& getName()        const { return _name; }
    const std::string& getType()        const { return _optionType; }
    const T&           getSpotPrice()   const { return _spotPrice; }
    const T&           getStrikePrice() const { return _strikePrice; }
    const T&           getRiskFreeRate()const { return _riskFreeRate; }
    const T&           getVolatility()  const { return _volatility; }
    const T&           getMaturityTime()const { return _maturityTime; }

    void info() const {
        std::cout << "======================================\n"
                  << "Spot Price     : " << _spotPrice << "\n"
                  << "Strike Price   : " << _strikePrice << "\n"
                  << "Risk-Free Rate : " << _riskFreeRate << "\n"
                  << "Volatility     : " << _volatility << "\n"
                  << "Maturity Time  : " << _maturityTime << " years\n"
                  << "======================================\n";
    }

protected:
    std::string _name = "European";
    std::string _optionType = "call";
    T _spotPrice = 100;
    T _strikePrice = 100;
    T _riskFreeRate = 0.0;
    T _volatility = 0.2;
    T _maturityTime = 1.0;
};

template <typename T>
class EuropeanOption : public Option<T> {
public:
    using Option<T>::Option; // inherit ctors

    const std::function<T(T)> payOff() const {
        const auto phi = (this->getType() == "call") ? 1.0 : -1.0;
        return [phi, K = this->getStrikePrice()](T ST) {
            return std::fmax(phi * (ST - K), 0.0);
        };
    }

    T closedForm(const std::string& optionType)
    {
        if (optionType != "call" && optionType != "put")
            throw std::invalid_argument("Invalid option type");
        this->setType(optionType);

        const double phi = (optionType == "call") ? +1.0 : -1.0;
        const double S0 = this->getSpotPrice();
        const double K  = this->getStrikePrice();
        const double r  = this->getRiskFreeRate();
        const double v  = this->getVolatility();
        const double Tm = this->getMaturityTime();

        const double d1 = (std::log(S0 / K) + (r + 0.5*v*v)*Tm) / (v*std::sqrt(Tm));
        const double d2 = d1 - v*std::sqrt(Tm);

        return phi * ( S0 * norm_cdf(phi * d1) - K * std::exp(-r*Tm) * norm_cdf(phi * d2) );
    }

    T computeDelta() {
        const double phi = (this->getType() == "call") ? +1.0 : -1.0;
        const double S0 = this->getSpotPrice();
        const double K  = this->getStrikePrice();
        const double r  = this->getRiskFreeRate();
        const double v  = this->getVolatility();
        const double Tm = this->getMaturityTime();

        const double d1 = (std::log(S0 / K) + (r + 0.5*v*v)*Tm) / (v*std::sqrt(Tm));
        return phi * norm_cdf(phi * d1);
    }

    T computeVega() {
        const double S0 = this->getSpotPrice();
        const double K  = this->getStrikePrice();
        const double r  = this->getRiskFreeRate();
        const double v  = this->getVolatility();
        const double Tm = this->getMaturityTime();

        const double d1 = (std::log(S0 / K) + (r + 0.5*v*v)*Tm) / (v*std::sqrt(Tm));
        return S0 * std::sqrt(Tm) * norm_pdf(d1);
    }

    void info() const {
        std::cout << "======================================\n"
                  << this->getName() << " Option Info\n"
                  << "Spot Price     : " << this->getSpotPrice() << "\n"
                  << "Strike Price   : " << this->getStrikePrice() << "\n"
                  << "Risk-Free Rate : " << this->getRiskFreeRate() << "\n"
                  << "Volatility     : " << this->getVolatility() << "\n"
                  << "Maturity Time  : " << this->getMaturityTime() << " years\n"
                  << "======================================\n";
    }
};

#endif // OPTION_HPP

