#include <../include/utilities.hpp>

template <typename T>
void Option<T>::setName(const std::string& name)
{
    _name = name;
}

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
const std::string& Option<T>::getName() const
{
    return _name;
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

template <typename T>
const std::function<T(T)> EuropeanOption<T>::payOff() const
{
    const auto phi = (this -> getType() == "call") ? 1.0 : -1.0;
    return [phi, strikePrice = this->getStrikePrice()](T finalPrice) {
        return std::fmax(phi * (finalPrice - strikePrice), 0.0);
    };
}
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

template <typename T>
void EuropeanOption<T>::info() const
{
    std::cout << "======================================\n"
              << this->getName() << " Option Info\n"
              << "Spot Price     : " << this->getSpotPrice() << "\n"
              << "Strike Price   : " << this->getStrikePrice() << "\n"
              << "Risk-Free Rate : " << this->getRiskFreeRate() << "\n"
              << "Volatility     : " << this->getVolatility() << "\n"
              << "Maturity Time  : " << this->getMaturityTime() << " years\n";
    std::cout << "======================================\n";
}

template <typename T>
T GBM<T>::advance(T S_t, T r, T sigma, T dt) const // one step GBM
{
    const T Z = static_cast<T>(standard_normal_sample());
    return S_t * std::exp((r - 0.5 * sigma * sigma) * dt + sigma * Z * std::sqrt(dt));
}

template <typename T>
T GBM<T>::advanceWithZ(T S_t, T r, T sigma, T dt, T Z) const
{
    return S_t * std::exp((r - 0.5 * sigma * sigma) * dt + sigma * Z * std::sqrt(dt));
}

template <typename T>
T GBM<T>::terminalDraw(T S0, T r, T sigma, T maturity) const
{
    // Implement the GBM terminal draw
    const T Z = static_cast<T>(standard_normal_sample());
    return S0 * std::exp((r - 0.5 * sigma * sigma) * maturity + sigma * Z * std::sqrt(maturity));
}

template <typename T>
T GBM<T>::terminalDrawWithZ(T S0, T r, T sigma, T maturity, T Z) const
{
    return S0 * std::exp((r - 0.5 * sigma * sigma) * maturity + sigma * Z * std::sqrt(maturity));
}

template <typename T>
T Euler<T>::advance(T S_t, T r, T sigma, T dt) const
{
    const T Z = static_cast<T>(standard_normal_sample());
    return S_t * (1.0 + r * dt + sigma * Z * std::sqrt(dt));
}

template <typename T>
T Euler<T>::advanceWithZ(T S_t, T r, T sigma, T dt, T Z) const
{
    return S_t * (1.0 + r * dt + sigma * Z * std::sqrt(dt));
}

template <typename T>
T Milstein<T>::advance(T S_t, T r, T sigma, T dt) const
{
    // Implement the Milstein advance step
    const T Z = static_cast<T>(standard_normal_sample());
    return  S_t * (1.0 + r * dt + sigma * Z * std::sqrt(dt) +
            0.5 * sigma * sigma * (Z * Z - 1.0) * dt);
}

template <typename T>
T Milstein<T>::advanceWithZ(T S_t, T r, T sigma, T dt, T Z) const
{
    return S_t * (1.0 + r * dt + sigma * Z * std::sqrt(dt) +
            0.5 * sigma * sigma * (Z * Z - 1.0) * dt);
}

template <typename T, typename Opt, typename Step>
void MonteCarlo<T, Opt, Step>::setNumPaths(std::size_t numPaths)
{
    _numPaths = numPaths;
}

template <typename T, typename Opt, typename Step>
void MonteCarlo<T, Opt, Step>::useAntitheticVariates(bool useAV)
{
    _useAV = useAV;
}

template <typename T, typename Opt, typename Step>
void MonteCarlo<T, Opt, Step>::useControlVariates(bool useCV)
{
    _useCV = useCV;
}

/*
// ---------- Full-path (discretized) Monte Carlo ----------
template <typename T, typename Opt, typename Step>
MCReturn<T> MonteCarlo<T, Opt, Step>::pricePath(const std::string& optionType, std::size_t numPaths)
{
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");
    if (numPaths <= 1)
        throw std::invalid_argument("Number of paths must be >= 2");

    _option.setType(optionType);
    setNumPaths(numPaths);

    auto payOffFun = _option.payOff();

    const T S0       = _option.getSpotPrice();
    const T r        = _option.getRiskFreeRate();
    const T sigma    = _option.getVolatility();
    const T maturity = _option.getMaturityTime();
    const T disc     = std::exp(-r * maturity);

    // time grid
    const std::size_t numSteps = std::max<std::size_t>(1, static_cast<std::size_t>(TRADING_DAYS_PER_YEAR * maturity));
    const T dt = maturity / static_cast<T>(numSteps);

    _stride = numSteps + 1;
    _pricePaths_flat.resize(numPaths * _stride, T(0));

    // No AV: identical to your current code (legs are i.i.d.)
    if (!_useAV) {
        long N = 0;
        T sumY = 0, sumY2 = 0;
        T sumX = 0, sumX2 = 0, sumXY = 0;

        for (std::size_t i = 0; i < numPaths; ++i) {
            T S = S0;
            const std::size_t base = i * _stride;
            _pricePaths_flat[base + 0] = S0;

            for (std::size_t j = 1; j <= numSteps; ++j) {
                S = _stepper.advance(S, r, sigma, dt);
                _pricePaths_flat[base + j] = S;
            }

            const T Y = disc * payOffFun(S);
            sumY  += Y;
            sumY2 += Y*Y;

            if (_useCV) {
                const T X = disc * S;
                sumX  += X;
                sumX2 += X*X;
                sumXY += Y*X;
            }
        }
        N = static_cast<long>(numPaths);

        const T Nf = static_cast<T>(N);
        const T EY   = sumY / Nf;
        const T varY = (sumY2 - Nf * EY * EY) / (Nf - T(1));

        if (!_useCV) return { EY, varY, std::sqrt(varY / Nf) };

        const T EX    = sumX / Nf;
        const T varX  = (sumX2 - Nf * EX * EX) / (Nf - T(1));
        const T covXY = (sumXY - Nf * EY * EX) / (Nf - T(1));

        const T beta    = (varX > T(0)) ? (covXY / varX) : T(0);
        const T EY_adj  = EY - beta * (EX - S0);
        const T var_adj = varY - ((varX > T(0)) ? (covXY * covXY / varX) : T(0));

        return { EY_adj, var_adj, std::sqrt(var_adj / Nf) };
    }

    // ===== AV branch: work with PAIRS =====
    const std::size_t pairs = (numPaths / 2);
    long Npairs = static_cast<long>(pairs);

    T sumA = 0, sumA2 = 0;               // A = (Y1+Y2)/2
    T sumXb = 0, sumXb2 = 0, sumAXb = 0; // Xb = (X1+X2)/2

    for (std::size_t i = 0; i < pairs; ++i) {
        // Path 1
        T S1 = S0;
        const std::size_t base1 = (2*i) * _stride;
        _pricePaths_flat[base1 + 0] = S1;

        // Path 2 (antithetic)
        T S2 = S0;
        const std::size_t base2 = (2*i + 1) * _stride;
        _pricePaths_flat[base2 + 0] = S2;

        for (std::size_t j = 1; j <= numSteps; ++j) {
            const T Z = static_cast<T>(standard_normal_sample());
            S1 = _stepper.advanceWithZ(S1, r, sigma, dt, +Z);
            S2 = _stepper.advanceWithZ(S2, r, sigma, dt, -Z);
            _pricePaths_flat[base1 + j] = S1;
            _pricePaths_flat[base2 + j] = S2;
        }

        const T Y1 = disc * payOffFun(S1);
        const T Y2 = disc * payOffFun(S2);
        const T A  = (Y1 + Y2) * T(0.5);   // pair average

        sumA  += A;
        sumA2 += A*A;

        if (_useCV) {
            const T X1  = disc * S1;      // control variate per leg
            const T X2  = disc * S2;
            const T Xb  = (X1 + X2) * T(0.5); // pair-averaged CV

            sumXb  += Xb;
            sumXb2 += Xb*Xb;
            sumAXb += A*Xb;
        }
    }

    const T Nf = static_cast<T>(Npairs);
    const T EA   = sumA / Nf;
    const T varA = (sumA2 - Nf * EA * EA) / (Nf - T(1));

    if (!_useCV) {
        return { EA, varA, std::sqrt(varA / Nf) };
    }

    const T EXb   = sumXb / Nf;
    const T varXb = (sumXb2 - Nf * EXb * EXb) / (Nf - T(1));
    const T covAX = (sumAXb - Nf * EA * EXb) / (Nf - T(1));

    const T beta    = (varXb > T(0)) ? (covAX / varXb) : T(0);
    const T EA_adj  = EA - beta * (EXb - S0);                 // E[X]=S0
    const T var_adj = varA - ((varXb > T(0)) ? (covAX * covAX / varXb) : T(0));

    return { EA_adj, var_adj, std::sqrt(var_adj / Nf) };
}
*/

// ---------- Strong error on STATE vs exact GBM (CRN) ----------
template <typename T, typename Opt, typename Step>
StrongStats MonteCarlo<T, Opt, Step>::strongErrorOnState(std::size_t numPaths,
                                                         std::size_t numSteps)
{
    if (numPaths <= 1) throw std::invalid_argument("numPaths >= 2");
    if (numSteps == 0) throw std::invalid_argument("numSteps >= 1");

    const T S0   = _option.getSpotPrice();
    const T r    = _option.getRiskFreeRate();
    const T sig  = _option.getVolatility();
    const T Tmat = _option.getMaturityTime();
    const T dt   = Tmat / static_cast<T>(numSteps);

    // No AV/CV here by design
    T sumAbs = 0, sumAbs2 = 0, sumSq = 0;

    for (std::size_t i = 0; i < numPaths; ++i) {
        T S  = S0;  // approximate path terminal
        T WT = 0;   // to build exact terminal from SAME shocks (A shock is
                    // a standard normal sample, Z)

        for (std::size_t j = 0; j < numSteps; ++j) {
            const T Z = static_cast<T>(standard_normal_sample());
            S  = _stepper.advanceWithZ(S, r, sig, dt, Z);
            WT += std::sqrt(dt) * Z;
        }
        const T Sexact = S0 * std::exp((r - T(0.5)*sig*sig)*Tmat + sig*WT);

        const T diff  = Sexact - S;
        const T adiff = std::fabs(diff);

        sumAbs  += adiff; // absolute error
        sumAbs2 += adiff*adiff; // squared absolute error
        sumSq   += diff*diff; // squared error
    }

    const T Nf     = static_cast<T>(numPaths); // number of paths
    const T meanAbs= sumAbs / Nf; // mean absolute error
    const T varAbs = (sumAbs2 - Nf * meanAbs * meanAbs) / (Nf - T(1)); // variance of absolute error
    const T seAbs  = std::sqrt( varAbs / Nf ); // standard error of absolute error
    const T rms    = std::sqrt( sumSq / Nf ); // root mean square error

    return { static_cast<double>(meanAbs),
             static_cast<double>(rms),
             static_cast<double>(seAbs) };
}

// ---------- Weak error (payoff bias) vs exact GBM (CRN) ----------
template <typename T, typename Opt, typename Step>
WeakStats MonteCarlo<T, Opt, Step>::weakErrorOnPayoff(const std::string& optionType,
                                                      std::size_t numPaths,
                                                      std::size_t numSteps)
{
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");
    if (numPaths <= 1) throw std::invalid_argument("numPaths >= 2");
    if (numSteps == 0) throw std::invalid_argument("numSteps >= 1");

    _option.setType(optionType);

    const auto g  = _option.payOff();
    const T S0    = _option.getSpotPrice();
    const T r     = _option.getRiskFreeRate();
    const T sig   = _option.getVolatility();
    const T Tmat  = _option.getMaturityTime();
    const T dt    = Tmat / static_cast<T>(numSteps);
    const T disc  = std::exp(-r * Tmat);

    // No AV/CV here by design
    T sumYh = 0, sumYh2 = 0;
    T sumYs = 0, sumYs2 = 0;
    T sumD  = 0, sumD2  = 0;

    for (std::size_t i = 0; i < numPaths; ++i) {
        T S  = S0;  // approximate terminal
        T WT = 0;   // exact via same shocks

        for (std::size_t j = 0; j < numSteps; ++j) {
            const T Z = static_cast<T>(standard_normal_sample());
            S  = _stepper.advanceWithZ(S, r, sig, dt, Z);
            WT += std::sqrt(dt) * Z;
        }
        const T Sexact = S0 * std::exp((r - T(0.5)*sig*sig)*Tmat + sig*WT);

        const T Yh = disc * g(S);
        const T Ys = disc * g(Sexact);
        const T D  = Yh - Ys;

        sumYh += Yh; sumYh2 += Yh*Yh;
        sumYs += Ys; sumYs2 += Ys*Ys;
        sumD  += D;  sumD2  += D*D;
    }

    const T Nf    = static_cast<T>(numPaths);
    const T EYh   = sumYh / Nf;
    const T EYs   = sumYs / Nf;
    const T ED    = sumD  / Nf;
    const T varD  = (sumD2 - Nf * ED * ED) / (Nf - T(1));
    const T seD   = std::sqrt( varD / Nf );

    WeakStats out;
    out.bias       = std::fabs(static_cast<double>(ED));
    out.seBias     = static_cast<double>(seD);
    out.meanApprox = static_cast<double>(EYh);
    out.meanExact  = static_cast<double>(EYs);
    return out;
}

// === NEW: Stepper-agnostic statistical convergence via pricePath ===
template <typename T, typename Opt, typename Step>
MCStatConvergence
MonteCarlo<T, Opt, Step>::statisticalConvergenceByPaths(const std::string& optionType,
                                                        const std::vector<std::size_t>& Ns,
                                                        bool useAV,
                                                        bool useCV)
{
    // Configure VR flags locally for the study
    this->useAntitheticVariates(useAV);
    this->useControlVariates(useCV);

    MCStatConvergence out;
    out.Ns = Ns;
    out.prices.resize(Ns.size());
    out.stdErrs.resize(Ns.size());

    std::vector<double> xN; xN.reserve(Ns.size());
    std::vector<double> se; se.reserve(Ns.size());

    for (std::size_t i = 0; i < Ns.size(); ++i) {
        const std::size_t N = Ns[i];

        // pricePath works for ANY stepper (GBM/Euler/Milstein); it internally sets a daily grid.
        auto r = this->pricePath(optionType, N);

        out.prices[i] = static_cast<double>(r.price);
        out.stdErrs[i]= static_cast<double>(r.stdDev);

        xN.push_back(static_cast<double>(N));
        se.push_back(out.stdErrs[i]);
    }

    // Expect ~ -0.5
    out.slope_loglog = slope_loglog_xy(xN, se);
    return out;
}

// ---------- Full-path (discretized) Monte Carlo ----------
template <typename T, typename Opt, typename Step>
MCReturn<T> MonteCarlo<T, Opt, Step>::pricePath(const std::string& optionType, std::size_t numPaths)
{
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");
    if (numPaths <= 1)
        throw std::invalid_argument("Number of paths must be >= 2");

    _option.setType(optionType);
    setNumPaths(numPaths);

    const auto payOffFun = _option.payOff();

    const T S0       = _option.getSpotPrice();
    const T r        = _option.getRiskFreeRate();
    const T sigma    = _option.getVolatility();
    const T Tmat     = _option.getMaturityTime();
    const T disc     = std::exp(-r * Tmat);

    // Time grid
    const std::size_t numSteps = std::max<std::size_t>(1, static_cast<std::size_t>(TRADING_DAYS_PER_YEAR * Tmat));
    const T dt = Tmat / static_cast<T>(numSteps);

    _stride = numSteps + 1;
    _pricePaths_flat.resize(numPaths * _stride, T(0));

    if (!_useAV) {
        // ----- No AV: leg-wise i.i.d. sampling -----
        long N = 0;
        T sumY = 0, sumY2 = 0;
        T sumX = 0, sumX2 = 0, sumXY = 0;

        for (std::size_t i = 0; i < numPaths; ++i) {
            T S = S0;
            T WT = 0; // accumulate Brownian motion using the SAME Zs we use to step
            const std::size_t base = i * _stride;
            _pricePaths_flat[base + 0] = S0;

            for (std::size_t j = 1; j <= numSteps; ++j) {
                const T Z = static_cast<T>(standard_normal_sample());
                S  = _stepper.advanceWithZ(S, r, sigma, dt, Z);
                WT += std::sqrt(dt) * Z;
                _pricePaths_flat[base + j] = S;
            }

            const T Y = disc * payOffFun(S);
            sumY  += Y;
            sumY2 += Y*Y;

            if (_useCV) {
                // Perfect CV from the SAME shocks
                const T ST_exact = S0 * std::exp((r - T(0.5)*sigma*sigma)*Tmat + sigma*WT); // I know exact ST under GBM
                const T X = disc * ST_exact; // E[X] = S0 exactly
                sumX  += X;  sumX2 += X*X;  sumXY += Y*X;
            }
        }
        N = static_cast<long>(numPaths);

        const T Nf = static_cast<T>(N);
        const T EY   = sumY / Nf;
        const T varY = (sumY2 - Nf * EY * EY) / (Nf - T(1));

        if (!_useCV) return { EY, varY, std::sqrt(varY / Nf) };

        const T EX    = sumX / Nf;
        const T varX  = (sumX2 - Nf * EX * EX) / (Nf - T(1));
        const T covXY = (sumXY - Nf * EY * EX) / (Nf - T(1));
        const T beta  = (varX > T(0)) ? (covXY / varX) : T(0);

        const T EY_adj  = EY - beta * (EX - S0);
        const T var_adj = varY - ((varX > T(0)) ? (covXY * covXY / varX) : T(0));
        return { EY_adj, var_adj, std::sqrt(var_adj / Nf) };
    }

    // ----- AV branch: PAIR-BASED statistics over pair means -----
    const std::size_t pairs = (numPaths / 2); // if odd, last leg is dropped intentionally
    if (pairs == 0) throw std::invalid_argument("Need at least 2 paths when AV is on.");
    long Npairs = static_cast<long>(pairs);

    T sumA = 0, sumA2 = 0;               // A = (Y1+Y2)/2
    T sumXb = 0, sumXb2 = 0, sumAXb = 0; // Xb = (X1+X2)/2

    for (std::size_t i = 0; i < pairs; ++i) {
        // Leg 1
        T S1 = S0, WT1 = 0;
        const std::size_t base1 = (2*i) * _stride;
        _pricePaths_flat[base1 + 0] = S1;

        // Leg 2 (antithetic)
        T S2 = S0, WT2 = 0;
        const std::size_t base2 = (2*i + 1) * _stride;
        _pricePaths_flat[base2 + 0] = S2;

        for (std::size_t j = 1; j <= numSteps; ++j) {
            const T Z = static_cast<T>(standard_normal_sample());
            S1  = _stepper.advanceWithZ(S1, r, sigma, dt, +Z);
            S2  = _stepper.advanceWithZ(S2, r, sigma, dt, -Z);
            WT1 += std::sqrt(dt) * (+Z);
            WT2 += std::sqrt(dt) * (-Z);
            _pricePaths_flat[base1 + j] = S1;
            _pricePaths_flat[base2 + j] = S2;
        }

        const T Y1 = disc * payOffFun(S1);
        const T Y2 = disc * payOffFun(S2);
        const T A  = T(0.5) * (Y1 + Y2); // pair mean payoff

        sumA  += A;
        sumA2 += A*A;

        if (_useCV) {
            // Perfect pair-averaged control variate Xb
            const T ST1_exact = S0 * std::exp((r - T(0.5)*sigma*sigma)*Tmat + sigma*WT1);
            const T ST2_exact = S0 * std::exp((r - T(0.5)*sigma*sigma)*Tmat + sigma*WT2);
            const T X1 = disc * ST1_exact;
            const T X2 = disc * ST2_exact;
            const T Xb = T(0.5) * (X1 + X2);

            sumXb  += Xb;
            sumXb2 += Xb*Xb;
            sumAXb += A*Xb;
        }
    }

    const T Nf = static_cast<T>(Npairs);
    const T EA   = sumA / Nf;
    const T varA = (sumA2 - Nf * EA * EA) / (Nf - T(1));

    if (!_useCV) return { EA, varA, std::sqrt(varA / Nf) };

    const T EXb   = sumXb / Nf;
    const T varXb = (sumXb2 - Nf * EXb * EXb) / (Nf - T(1));
    const T covAX = (sumAXb - Nf * EA * EXb) / (Nf - T(1));
    const T beta  = (varXb > T(0)) ? (covAX / varXb) : T(0);

    const T EA_adj  = EA - beta * (EXb - S0);                 // E[X]=S0
    const T var_adj = varA - ((varXb > T(0)) ? (covAX * covAX / varXb) : T(0));
    return { EA_adj, var_adj, std::sqrt(var_adj / Nf) };
}

/*
// ---------- Terminal (single-step, exact GBM) ----------
template <typename T, typename Opt, typename Step>
MCReturn<T> MonteCarlo<T, Opt, Step>::priceTerminal(const std::string& optionType, std::size_t numPaths)
{
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");
    if (numPaths <= 1)
        throw std::invalid_argument("Number of paths must be >= 2 for priceTerminal().");

    _option.setType(optionType);
    setNumPaths(numPaths);

    const auto payOffFun = _option.payOff();
    const T S0       = _option.getSpotPrice();
    const T r        = _option.getRiskFreeRate();
    const T sigma    = _option.getVolatility();
    const T maturity = _option.getMaturityTime();
    const T disc     = std::exp(-r * maturity);

    if (!_useAV) {
        long N = 0;
        T sumY = 0, sumY2 = 0;
        T sumX = 0, sumX2 = 0, sumXY = 0;

        for (std::size_t i = 0; i < _numPaths; ++i) {
            const T Z  = static_cast<T>(standard_normal_sample());
            const T ST = _stepper.terminalDrawWithZ(S0, r, sigma, maturity, Z);

            const T Y = disc * payOffFun(ST);
            sumY  += Y;     sumY2 += Y*Y;

            if (_useCV) {
                const T X = disc * ST;
                sumX  += X;  sumX2 += X*X;  sumXY += Y*X;
            }
        }
        N = static_cast<long>(_numPaths);

        const T Nf = static_cast<T>(N);
        const T EY   = sumY / Nf;
        const T varY = (sumY2 - Nf * EY * EY) / (Nf - T(1));

        if (!_useCV) return { EY, varY, std::sqrt(varY / Nf) };

        const T EX    = sumX / Nf;
        const T varX  = (sumX2 - Nf * EX * EX) / (Nf - T(1));
        const T covXY = (sumXY - Nf * EY * EX) / (Nf - T(1));

        const T beta    = (varX > T(0)) ? (covXY / varX) : T(0);
        const T EY_adj  = EY - beta * (EX - S0);
        const T var_adj = varY - ((varX > T(0)) ? (covXY * covXY / varX) : T(0));

        return { EY_adj, var_adj, std::sqrt(var_adj / Nf) };
    }

    // ===== AV branch: pairs =====
    const std::size_t pairs = (_numPaths / 2);
    long Npairs = static_cast<long>(pairs);

    T sumA = 0, sumA2 = 0;
    T sumXb = 0, sumXb2 = 0, sumAXb = 0;

    for (std::size_t i = 0; i < pairs; ++i) {
        const T Z   = static_cast<T>(standard_normal_sample());
        const T ST1 = _stepper.terminalDrawWithZ(S0, r, sigma, maturity, +Z);
        const T ST2 = _stepper.terminalDrawWithZ(S0, r, sigma, maturity, -Z);

        const T Y1 = disc * payOffFun(ST1);
        const T Y2 = disc * payOffFun(ST2);
        const T A  = (Y1 + Y2) * T(0.5);

        sumA  += A;
        sumA2 += A*A;

        if (_useCV) {
            const T X1 = disc * ST1;
            const T X2 = disc * ST2;
            const T Xb = (X1 + X2) * T(0.5);

            sumXb  += Xb;
            sumXb2 += Xb*Xb;
            sumAXb += A*Xb;
        }
    }

    const T Nf = static_cast<T>(Npairs);
    const T EA   = sumA / Nf;
    const T varA = (sumA2 - Nf * EA * EA) / (Nf - T(1));

    if (!_useCV) return { EA, varA, std::sqrt(varA / Nf) };

    const T EXb   = sumXb / Nf;
    const T varXb = (sumXb2 - Nf * EXb * EXb) / (Nf - T(1));
    const T covAX = (sumAXb - Nf * EA * EXb) / (Nf - T(1));

    const T beta    = (varXb > T(0)) ? (covAX / varXb) : T(0);
    const T EA_adj  = EA - beta * (EXb - S0);
    const T var_adj = varA - ((varXb > T(0)) ? (covAX * covAX / varXb) : T(0));

    return { EA_adj, var_adj, std::sqrt(var_adj / Nf) };
}
*/
// ---------- Terminal (single-step, exact GBM) ----------
template <typename T, typename Opt, typename Step>
MCReturn<T> MonteCarlo<T, Opt, Step>::priceTerminal(const std::string& optionType, std::size_t numPaths)
{
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");
    if (numPaths <= 1)
        throw std::invalid_argument("Number of paths must be >= 2 for priceTerminal().");

    _option.setType(optionType);
    setNumPaths(numPaths);

    const auto payOffFun = _option.payOff();
    const T S0       = _option.getSpotPrice();
    const T r        = _option.getRiskFreeRate();
    const T sigma    = _option.getVolatility();
    const T Tmat     = _option.getMaturityTime();
    const T disc     = std::exp(-r * Tmat);

    if (!_useAV) {
        long N = 0;
        T sumY = 0, sumY2 = 0;
        T sumX = 0, sumX2 = 0, sumXY = 0;

        for (std::size_t i = 0; i < _numPaths; ++i) {
            const T Z  = static_cast<T>(standard_normal_sample());
            const T ST = _stepper.terminalDrawWithZ(S0, r, sigma, Tmat, Z);

            const T Y = disc * payOffFun(ST);
            sumY  += Y;     sumY2 += Y*Y;

            if (_useCV) {
                const T X = disc * ST;  // exact terminal sampling -> E[X]=S0
                sumX  += X;  sumX2 += X*X;  sumXY += Y*X;
            }
        }
        N = static_cast<long>(_numPaths);

        const T Nf = static_cast<T>(N);
        const T EY   = sumY / Nf;
        const T varY = (sumY2 - Nf * EY * EY) / (Nf - T(1));

        if (!_useCV) return { EY, varY, std::sqrt(varY / Nf) };

        const T EX    = sumX / Nf;
        const T varX  = (sumX2 - Nf * EX * EX) / (Nf - T(1));
        const T covXY = (sumXY - Nf * EY * EX) / (Nf - T(1));
        const T beta  = (varX > T(0)) ? (covXY / varX) : T(0);

        const T EY_adj  = EY - beta * (EX - S0);
        const T var_adj = varY - ((varX > T(0)) ? (covXY * covXY / varX) : T(0));
        return { EY_adj, var_adj, std::sqrt(var_adj / Nf) };
    }

    // AV branch: pair means
    const std::size_t pairs = (_numPaths / 2); // if odd, last leg is dropped intentionally
    if (pairs == 0) throw std::invalid_argument("Need at least 2 paths when AV is on.");
    long Npairs = static_cast<long>(pairs);

    T sumA = 0, sumA2 = 0;
    T sumXb = 0, sumXb2 = 0, sumAXb = 0;

    for (std::size_t i = 0; i < pairs; ++i) {
        const T Z   = static_cast<T>(standard_normal_sample());
        const T ST1 = _stepper.terminalDrawWithZ(S0, r, sigma, Tmat, +Z);
        const T ST2 = _stepper.terminalDrawWithZ(S0, r, sigma, Tmat, -Z);

        const T Y1 = disc * payOffFun(ST1);
        const T Y2 = disc * payOffFun(ST2);
        const T A  = T(0.5) * (Y1 + Y2);

        sumA  += A;
        sumA2 += A*A;

        if (_useCV) {
            const T X1 = disc * ST1;
            const T X2 = disc * ST2;
            const T Xb = T(0.5) * (X1 + X2);

            sumXb  += Xb;
            sumXb2 += Xb*Xb;
            sumAXb += A*Xb;
        }
    }

    const T Nf = static_cast<T>(Npairs);
    const T EA   = sumA / Nf;
    const T varA = (sumA2 - Nf * EA * EA) / (Nf - T(1));

    if (!_useCV) return { EA, varA, std::sqrt(varA / Nf) };

    const T EXb   = sumXb / Nf;
    const T varXb = (sumXb2 - Nf * EXb * EXb) / (Nf - T(1));  // sample var of Xb
    const T covAX = (sumAXb - Nf * EA * EXb) / (Nf - T(1));  // sample cov(A, Xb)
    const T beta  = (varXb > T(0)) ? (covAX / varXb) : T(0);

    const T EA_adj  = EA - beta * (EXb - S0);
    const T var_adj = varA - ((varXb > T(0)) ? (covAX * covAX / varXb) : T(0));
    return { EA_adj, var_adj, std::sqrt(var_adj / Nf) };
}



// // Pricing with direct terminal sampling (single step, exact GBM)
// template <typename T, typename Opt, typename Step>
// MCReturn<T> MonteCarlo<T, Opt, Step>::priceTerminal(const std::string& optionType, std::size_t numPaths)
// {
//     // Validate inputs
//     if (optionType != "call" && optionType != "put")
//         throw std::invalid_argument("Invalid option type");
//     if (numPaths <= 1)
//         throw std::invalid_argument("Number of paths must be >= 2 for priceTerminal().");

//     // Apply inputs to engine state
//     _option.setType(optionType);
//     setNumPaths(numPaths);

//     // Get payoff and contract params
//     const auto payOffFun = _option.payOff();
//     const T S0       = _option.getSpotPrice();
//     const T r        = _option.getRiskFreeRate();
//     const T sigma    = _option.getVolatility();
//     const T maturity = _option.getMaturityTime();
//     const T disc     = std::exp(-r * maturity);

//     // Accumulators: Y = disc * payoff(S_T), control variate X = disc * S_T (E[X]=S0)
//     long N = 0;
//     T sumY = T(0), sumY2 = T(0);
//     T sumX = T(0), sumX2 = T(0), sumXY = T(0);

//     if (_useAV) {
//         const std::size_t half = _numPaths / 2;
//         for (std::size_t i = 0; i < half; ++i) {
//             const T Z   = static_cast<T>(standard_normal_sample());
//             const T ST1 = gbm_ST(S0, r, sigma, maturity, +Z);
//             const T ST2 = gbm_ST(S0, r, sigma, maturity, -Z);

//             const T Y1 = disc * payOffFun(ST1);
//             const T Y2 = disc * payOffFun(ST2);

//             sumY  += (Y1 + Y2);
//             sumY2 += (Y1*Y1 + Y2*Y2);

//             if (_useCV) {
//                 const T X1 = disc * ST1;   // E[X]=S0 under RN
//                 const T X2 = disc * ST2;
//                 sumX  += (X1 + X2);
//                 sumX2 += (X1*X1 + X2*X2);
//                 sumXY += (Y1*X1 + Y2*X2);
//             }
//         }
//         N = static_cast<long>(2 * half); // drop last if odd for pairing
//     } else {
//         for (std::size_t i = 0; i < _numPaths; ++i) {
//             const T Z  = static_cast<T>(standard_normal_sample());
//             const T ST = gbm_ST(S0, r, sigma, maturity, Z);

//             const T Y = disc * payOffFun(ST);
//             sumY  += Y;
//             sumY2 += Y*Y;

//             if (_useCV) {
//                 const T X = disc * ST;
//                 sumX  += X;
//                 sumX2 += X*X;
//                 sumXY += Y*X;
//             }
//         }
//         N = static_cast<long>(_numPaths);
//     }

//     const T Nf   = static_cast<T>(N);
//     const T EY   = sumY / Nf;
//     const T varY = (sumY2 - Nf * EY * EY) / (Nf - T(1));

//     if (!_useCV) {
//         return { EY, varY, std::sqrt(varY / Nf) };
//     }

//     // Control variate adjustment with X = disc * S_T and E[X]=S0
//     const T EX    = sumX / Nf;
//     const T varX  = (sumX2 - Nf * EX * EX) / (Nf - T(1));
//     const T covXY = (sumXY - Nf * EY * EX) / (Nf - T(1));

//     const T beta    = (varX > T(0)) ? (covXY / varX) : T(0);
//     const T EY_adj  = EY - beta * (EX - S0);
//     const T var_adj = varY - ((varX > T(0)) ? (covXY * covXY / varX) : T(0));

//     return { EY_adj, var_adj, std::sqrt(var_adj / Nf) };
// }

// // ------ Terminal (exact) pricing: now uses Step::terminalDraw ------
// template <typename T, typename Opt, typename Step>
// MCReturn<T> MonteCarlo<T, Opt, Step>::priceTerminal(const std::string& optionType, std::size_t numPaths)
// {
//     if (optionType != "call" && optionType != "put")
//         throw std::invalid_argument("Invalid option type");
//     if (numPaths <= 1)
//         throw std::invalid_argument("Number of paths must be >= 2 for priceTerminal().");

//     _option.setType(optionType);
//     setNumPaths(numPaths);

//     const auto payOffFun = _option.payOff();
//     const T S0       = _option.getSpotPrice();
//     const T r        = _option.getRiskFreeRate();
//     const T sigma    = _option.getVolatility();
//     const T maturity = _option.getMaturityTime();
//     const T disc     = std::exp(-r * maturity);

//     long N = 0;
//     T sumY = T(0), sumY2 = T(0);
//     T sumX = T(0), sumX2 = T(0), sumXY = T(0);

//     if (_useAV) {
//         // Antithetic flag keeps pair-wise accumulation. We still draw via stepper twice.
//         const std::size_t half = _numPaths / 2;
//         for (std::size_t i = 0; i < half; ++i) {
//             const T Z = static_cast<T>(standard_normal_sample());
//             const T ST1 = _stepper.terminalDrawWithZ(S0, r, sigma, maturity, Z);
//             const T ST2 = _stepper.terminalDrawWithZ(S0, r, sigma, maturity, -Z);

//             const T Y1 = disc * payOffFun(ST1);
//             const T Y2 = disc * payOffFun(ST2);

//             sumY  += (Y1 + Y2);
//             sumY2 += (Y1*Y1 + Y2*Y2);

//             if (_useCV) {
//                 const T X1 = disc * ST1;
//                 const T X2 = disc * ST2;
//                 sumX  += (X1 + X2);
//                 sumX2 += (X1*X1 + X2*X2);
//                 sumXY += (Y1*X1 + Y2*X2);
//             }
//         }
//         N = static_cast<long>(2 * half);
//     } else {
//         for (std::size_t i = 0; i < _numPaths; ++i) {
//             const T ST = _stepper.terminalDraw(S0, r, sigma, maturity);

//             const T Y = disc * payOffFun(ST);
//             sumY  += Y;
//             sumY2 += Y*Y;

//             if (_useCV) {
//                 const T X = disc * ST;
//                 sumX  += X;
//                 sumX2 += X*X;
//                 sumXY += Y*X;
//             }
//         }
//         N = static_cast<long>(_numPaths);
//     }

//     const T Nf   = static_cast<T>(N);
//     const T EY   = sumY / Nf;
//     const T varY = (sumY2 - Nf * EY * EY) / (Nf - T(1));

//     if (!_useCV) {
//         return { EY, varY, std::sqrt(varY / Nf) };
//     }

//     const T EX    = sumX / Nf;
//     const T varX  = (sumX2 - Nf * EX * EX) / (Nf - T(1));
//     const T covXY = (sumXY - Nf * EY * EX) / (Nf - T(1));

//     const T beta    = (varX > T(0)) ? (covXY / varX) : T(0);
//     const T EY_adj  = EY - beta * (EX - S0); // E[X] = S0
//     const T var_adj = varY - ((varX > T(0)) ? (covXY * covXY / varX) : T(0));

//     return { EY_adj, var_adj, std::sqrt(var_adj / Nf) };
// }


template class Option<double>;
template class EuropeanOption<double>;
template class GBM<double>;
template class MonteCarlo<double, EuropeanOption<double>, GBM<double>>;
template class Euler<double>;
template class MonteCarlo<double, EuropeanOption<double>, Euler<double>>;
template class Milstein<double>;
template class MonteCarlo<double, EuropeanOption<double>, Milstein<double>>;


