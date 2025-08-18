#include <../include/utilities.hpp>

template <typename T>
void Option<T>::setType(const std::string& optionType)
{
    _optionType = optionType;
}

template <typename T>
void Option<T>::setSpotPrice(T spotPrice)
{
    _spotPrice = spotPrice;
}

template <typename T>
void Option<T>::setStrikePrice(T strikePrice)
{
    _strikePrice = strikePrice;
}

template <typename T>
void Option<T>::setRiskFreeRate(T riskFreeRate)
{
    _riskFreeRate = riskFreeRate;
}

template <typename T>
void Option<T>::setVolatility(T volatility)
{
    _volatility = volatility;
}

template <typename T>
void Option<T>::setMaturityTime(T maturityTime)
{
    _maturityTime = maturityTime;
}

template <typename T>
const std::string& Option<T>::getType() const
{
    return _optionType;
}

template <typename T>
const T& Option<T>::getSpotPrice() const
{
    return _spotPrice;
}

template <typename T>
const T& Option<T>::getStrikePrice() const
{
    return _strikePrice;
}

template <typename T>
const T& Option<T>::getRiskFreeRate() const
{
    return _riskFreeRate;
}

template <typename T>
const T& Option<T>::getVolatility() const
{
    return _volatility;
}

template <typename T>
const T& Option<T>::getMaturityTime() const
{
    return _maturityTime;
}

template <typename T>
void Option<T>::info() const
{
    std::cout << "======================================\n"
              << "Spot Price     : " << _spotPrice << "\n"
              << "Strike Price   : " << _strikePrice << "\n"
              << "Risk-Free Rate : " << _riskFreeRate << "\n"
              << "Volatility     : " << _volatility << "\n"
              << "Maturity Time  : " << _maturityTime << " years\n";
    std::cout << "======================================\n";
}

// template <typename T>
// void EuropeanOption<T>::setNumPaths(const unsigned long& numPaths)
// {
//     _numPaths = numPaths;
// }

// template <typename T>
// unsigned long EuropeanOption<T>::getNumPaths() const
// {
//     return _numPaths;
// }

template <typename T>
T EuropeanOption<T>::closedForm(const std::string& optionType)
{
    // check if the string is a valid one
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");

    this->setType(optionType);

    // map to pay off phi: call = +1, put = -1
    const double phi = (optionType == "call") ? 1.0 : -1.0;

    return phi * (
           this->getSpotPrice() * norm_cdf(
               phi * (
                   (std::log(this->getSpotPrice() / this->getStrikePrice()) +
                   (this->getRiskFreeRate() + 0.5 * this->getVolatility() * this->getVolatility()) * this->getMaturityTime()) /
                   (this->getVolatility() * std::sqrt(this->getMaturityTime()))
               )
           )
           - this->getStrikePrice() * std::exp(-this->getRiskFreeRate() * this->getMaturityTime()) *
           norm_cdf(
               phi * (
                   (std::log(this->getSpotPrice() / this->getStrikePrice()) +
                   (this->getRiskFreeRate() - 0.5 * this->getVolatility() * this->getVolatility()) * this->getMaturityTime()) /
                   (this->getVolatility() * std::sqrt(this->getMaturityTime()))
               )
           )
    );

}

/*
    First Naive implementation
*/
// template <typename T>
// T EuropeanOption<T>::naiveMonteCarloGBM(const std::string& optionType, int numPaths)
// {
//     if (optionType != "call" && optionType != "put")
//         throw std::invalid_argument("Invalid option type");

//     this->setType(optionType);

//     if (numPaths <= 0)
//         throw std::invalid_argument("Number of paths must be positive");

//     setNumPaths(numPaths);

//     // Get option parameters
//     const double S0       = this->getSpotPrice();
//     const double K        = this->getStrikePrice();
//     const double r        = this->getRiskFreeRate();
//     const double sigma    = this->getVolatility();
//     const double maturity = this->getMaturityTime();

//     // Set the number of time steps
//     const std::size_t numSteps = std::max<std::size_t>(1, static_cast<std::size_t>(TRADING_DAYS_PER_YEAR * maturity));
//     // Calculate the time step size
//     const double dt = maturity / static_cast<double>(numSteps);
//     const double sqrtDt = std::sqrt(dt);

//     _pricePaths.clear(); // Clear any existing price paths
//     _pricePaths.resize(_numPaths, std::vector<T>(numSteps, 0.0)); // Resize to accommodate the number of paths

//     // Generate price paths
//     // update rule is: 
//     // S_t 
//     for (std::size_t i = 0; i < _numPaths; ++i) // for each simulation path...
//     {
//         double Wt = 0.0; // initiall each path start with W_t = 0, also W_0 = 0
//         double t  = 0.0;

//         for (std::size_t j = 0; j < numSteps; ++j) // ..for each time step of a given simulation path
//         {
//             Wt += standard_normal_sample() * sqrtDt; // Brownian motion increment 
//             t = (j + 1) * dt; // (j+1) because in the first column we want the next time step update

//             _pricePaths[i][j] = S0 * std::exp((r - 0.5 * sigma * sigma) * t + sigma * Wt);
//         }
//     }

//     double cumPayoff = 0.0;

//     for (std::size_t i = 0; i < _numPaths; ++i) // for each simulation path
//     {
//         const double ST = _pricePaths[i].back(); // get the final price of the path, same as _pricePaths[i][numSteps - 1]
//         const double payoff = (optionType == "call") 
//             ? std::max(0.0, ST - K) // Call option payoff
//             : std::max(0.0, K - ST); // Put option payoff
//         cumPayoff += payoff; // accumulate the payoff
//     }

//     // return discounted average payoff
//     return std::exp(-r * maturity) * (cumPayoff / static_cast<double>(_numPaths));
// }

/*
    Second Naive implementation
*/
// template <typename T>
// T EuropeanOption<T>::naiveMonteCarloGBM(const std::string& optionType, int numPaths)
// {
//     if (optionType != "call" && optionType != "put")
//         throw std::invalid_argument("Invalid option type");

//     this->setType(optionType);

//     if (numPaths <= 0)
//         throw std::invalid_argument("Number of paths must be positive");

//     setNumPaths(numPaths);

//     // Get option parameters
//     const double S0       = this->getSpotPrice();
//     const double K        = this->getStrikePrice();
//     const double r        = this->getRiskFreeRate();
//     const double sigma    = this->getVolatility();
//     const double maturity = this->getMaturityTime();

//     // Set the number of time steps
//     const std::size_t numSteps = std::max<std::size_t>(1, static_cast<std::size_t>(TRADING_DAYS_PER_YEAR * maturity));
//     // Calculate the time step size
//     const double dt = maturity / static_cast<double>(numSteps);
//     const double sqrtDt = std::sqrt(dt);

//     _pricePaths.clear(); // Clear any existing price paths
//     _pricePaths.resize(_numPaths, std::vector<T>(numSteps + 1)); // Resize to accommodate the number of paths

//     // Generate price paths
//     double cumPayoff = 0.0;

//     for (std::size_t i = 0; i < _numPaths; ++i) // for each simulation path...
//     {
//         double S = S0; // initial price
//         _pricePaths[i][0] = S; // store the initial price

//         for (std::size_t j = 1; j <= numSteps; ++j) // ..for each time step of a given simulation path
//         {
//             const double Z = standard_normal_sample();
//             S *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrtDt * Z);
//             _pricePaths[i][j] = S; // store the price at the next time step
//         }

//         const double ST = S; // final price of the last step in the last path
//         cumPayoff += (optionType == "call") 
//                     ? std::max(0.0, ST - K) // Call option payoff
//                     : std::max(0.0, K - ST); // Put option payoff
//     }

//     // return discounted average payoff
//     return std::exp(-r * maturity) * (cumPayoff / static_cast<double>(_numPaths));
// }

/*
    Third Naive implementation
*/
template <typename T>
const MCReturn<T> EuropeanOption<T>::naiveMonteCarloGBM(const std::string& optionType, int numPaths)
{
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");

    this->setType(optionType);

    const double phi = (optionType == "call") ? 1.0 : -1.0;

    if (numPaths <= 1)
        throw std::invalid_argument("Number of paths must be positive and greater or equal to 2");

    // Get option parameters
    const double S0       = this->getSpotPrice();
    const double K        = this->getStrikePrice();
    const double r        = this->getRiskFreeRate();
    const double sigma    = this->getVolatility();
    const double maturity = this->getMaturityTime();

    // Set the number of time steps
    const std::size_t numSteps = std::max<std::size_t>(1, static_cast<std::size_t>(TRADING_DAYS_PER_YEAR * maturity));

    std::cout << "\nRunning Naive Monte Carlo GBM for " << numPaths << " paths and " << numSteps << " time steps...\n";

    // Calculate the time step size
    const double dt = maturity / static_cast<double>(numSteps);
    const double sqrtDt = std::sqrt(dt);

    _stride = numSteps + 1; // stride is the number of time steps + 1 for the initial price
    _pricePaths_flat.resize(numPaths * _stride, 0.0); // with .assign() you are initializing the vector with zeros, so you are clearing all previous data as well

    double sum = 0.0;
    double sum2 = 0.0;

    for (std::size_t i = 0; i < numPaths; ++i)
    {
        const std::size_t base = i * _stride;
        double S = S0; // initial price
        _pricePaths_flat[base] = S0;

        for (std::size_t j = 1; j <= numSteps; ++j)
        {
            S *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrtDt * standard_normal_sample());
            _pricePaths_flat[base + j] = S;
        }

        const double discountedPayoff = std::exp(-r * maturity) * std::fmax(0.0, phi * (S - K));
        sum += discountedPayoff;
        sum2 += discountedPayoff * discountedPayoff;
    }

    const double mean = sum / static_cast<double>(numPaths);
    const double variance = (sum2  - static_cast<double>(numPaths) * (mean * mean)) / (static_cast<double>(numPaths) - 1.0); // unbiased
    const double stdDev = std::sqrt(variance )/ std::sqrt(static_cast<double>(numPaths));

    // return discounted average payoff
    return { mean, variance, stdDev };
}

template <typename T>
const MCReturn<T> EuropeanOption<T>::fastLowVarMCGBM(const std::string& optionType, 
                                                     int numPaths, 
                                                     bool antithetic, 
                                                     bool controlVariates)
{
    /*
       Fast Low Variance Control Variate Monte Carlo for GBM
       Exploits:
        - equivalence in law --> no need to simulate the entire path
        - antithetic variates on the final payoff (variance reduction)
        - control variates on the initial price (variance reduction)
    */
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");

    this->setType(optionType);

    const double phi = (optionType == "call") ? 1.0 : -1.0;

    if (numPaths <= 1)
        throw std::invalid_argument("Number of paths must be positive and greater or equal to 2");
    
    std::cout << "\nRunning Fast Low Variance Monte Carlo GBM for " 
              << numPaths << " paths with antithetic = " << std::boolalpha << antithetic 
              << " and control variates = " << controlVariates << "...\n";

    // Get option parameters
    const double S0       = this->getSpotPrice();
    const double K        = this->getStrikePrice();
    const double r        = this->getRiskFreeRate();
    const double sigma    = this->getVolatility();
    const double maturity = this->getMaturityTime();

    // new variables
    const bool useAntithetic = antithetic;
    const bool useControlVariates = controlVariates;

    // stats accumulator
    double sumA=0.0, sumA2=0.0; // cumulative sum and sum of squares for antithetic variates
    double sumC=0.0, sumC2=0.0; // cumulative sum and sum of squares for control variates
    double sumAC=0.0; // cumulative sum for antithetic and control variates (covariance)

    int used = 0;
    if (useAntithetic)
    {
        const int half = numPaths / 2;
        used = 2 * half;

        for (std::size_t i = 0; i < half; ++i)
        {
            const double Z = standard_normal_sample();

            // Calculate the final price for the original path
            const double ST = gbm_ST(S0, r, sigma, maturity, Z);
            // Calculate the final price for the antithetic path
            const double ST_antithetic = gbm_ST(S0, r, sigma, maturity, -Z);

            // Calculate the discounted payoff for both paths
            const double discountedPayoff = std::exp(-r * maturity) * std::fmax(0.0, phi * (ST - K));
            const double discountedPayoff_antithetic = std::exp(-r * maturity) * std::fmax(0.0, phi * (ST_antithetic - K));

            // Calculate control variates
            const double controlVariate = std::exp(-r * maturity) * ST;
            const double controlVariate_antithetic = std::exp(-r * maturity) * ST_antithetic;

            sumA += discountedPayoff + discountedPayoff_antithetic;
            sumA2 += discountedPayoff * discountedPayoff + discountedPayoff_antithetic * discountedPayoff_antithetic;
            sumC += controlVariate + controlVariate_antithetic;
            sumC2 += controlVariate * controlVariate + controlVariate_antithetic * controlVariate_antithetic;
            sumAC += discountedPayoff * controlVariate + discountedPayoff_antithetic * controlVariate_antithetic;
        }
    } else {
        used = numPaths;

        for (std::size_t i = 0; i < numPaths; ++i)
        {
            const double Z = standard_normal_sample();

            // Calculate the final price
            const double ST = gbm_ST(S0, r, sigma, maturity, Z);

            // Calculate the discounted payoff
            const double discountedPayoff = std::exp(-r * maturity) * std::fmax(0.0, phi * (ST - K));

            // Calculate control variates
            const double controlVariate = std::exp(-r * maturity) * ST;

            sumA += discountedPayoff;
            sumA2 += discountedPayoff * discountedPayoff;
            sumC += controlVariate;
            sumC2 += controlVariate * controlVariate;
            sumAC += discountedPayoff * controlVariate;
        }
    }

    const double N = static_cast<double>(used);

    const double E_A = sumA / N; // computed price for antithetic variates (mean)
    const double E_C = sumC / N; // computed price for control variates (mean)

    // unbiased stats
    const double Var_A = (sumA2 - N * E_A * E_A) / (N - 1); // variance of the price for antithetic variates
    const double Var_C = (sumC2 - N * E_C * E_C) / (N - 1); // variance of the price for control variates
    const double Cov_AC = (sumAC - N * E_A * E_C) / (N - 1); // covariance between antithetic and control variates

    double price = E_A;
    double var = Var_A;

    if (useControlVariates && Var_C > 0.0)
    {
        // Adjust the price using control variates
        const double beta = Cov_AC / Var_C; // optimal beta
        const double adjustedPrice = E_A - beta * (E_C - S0);
        const double adjustedVar = Var_A - Cov_AC * Cov_AC / Var_C; // adjusted variance
        price = adjustedPrice;
        var = adjustedVar;
    }

    return { price, var, std::sqrt(var/N) }; // return the price, variance, and standard deviation
}

template <typename T>
void EuropeanOption<T>::info() const
{
    std::cout << "======================================\n"
              << "European Option Info:\n"
              << "Spot Price     : " << this->getSpotPrice() << "\n"
              << "Strike Price   : " << this->getStrikePrice() << "\n"
              << "Risk-Free Rate : " << this->getRiskFreeRate() << "\n"
              << "Volatility     : " << this->getVolatility() << "\n"
              << "Maturity Time  : " << this->getMaturityTime() << " years\n";
    std::cout << "======================================\n";
}

template class Option<double>;
template class EuropeanOption<double>;
