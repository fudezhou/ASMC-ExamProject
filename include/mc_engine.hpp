#ifndef MC_ENGINE_HPP
#define MC_ENGINE_HPP

#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <span>          // <-- added

#if defined(USE_OMP)
  #include <omp.h>
#endif

#include "math.hpp"
#include "mc_types.hpp"
#include "option.hpp"
#include "stepper.hpp"

template <typename T, typename Opt, typename Step>
class MonteCarlo {
public:
    MonteCarlo(Opt& option, const Step& stepper) : _option(option), _stepper(stepper) {}

    void setNumPaths(std::size_t numPaths)   { _numPaths = numPaths; }
    void useAntitheticVariates(bool useAV)   { _useAV = useAV; }
    void useControlVariates(bool useCV)      { _useCV = useCV; }

    // Full-path discretized pricing (any stepper)
    MCReturn<T> pricePath(const std::string& optionType, 
                          std::size_t numPaths, 
                          std::size_t numSteps = TRADING_DAYS_PER_YEAR,
                          bool storePaths = false);

    // ===== Terminal (single-step exact GBM) APIs =====

    // (A) Legacy overload — generates normals internally (backward-compatible)
    MCReturn<T> priceTerminal(const std::string& optionType, std::size_t numPaths);

    // (B) New overload — uses pregenerated normals via span (no copies)
    MCReturn<T> priceTerminal(const std::string& optionType,
                              std::size_t numPaths,
                              std::span<const T> preGeneratedZ);

    // (C) Convenience overload — accept std::vector and forward to span
    MCReturn<T> priceTerminal(const std::string& optionType,
                              std::size_t numPaths,
                              const std::vector<T>& preGeneratedZ)
    {
        return priceTerminal(optionType, numPaths,
                             std::span<const T>(preGeneratedZ.data(), preGeneratedZ.size()));
    }

    // Strong/Weak errors (CRN vs exact GBM)
    StrongStats strongErrorOnState(std::size_t numPaths, std::size_t numSteps);

    WeakStats   weakErrorOnPayoff(const std::string& optionType,
                                  std::size_t numPaths, std::size_t numSteps);

    // Empirical MC statistical convergence across N
    MCStatConvergence statisticalConvergenceByPaths(const std::string& optionType,
                                                    const std::vector<std::size_t>& Ns,
                                                    bool useAV, bool useCV);

    // (optional) flattened path storage (per path, stride = steps+1)
    const std::vector<T>& pricePathsFlat() const { return _pricePaths_flat; }
    const std::vector<T>& getPricePathsFlat() const { return _pricePaths_flat; }

    MCgreeksFD computeGreeksFD(const std::string& optionType,
                               std::size_t numPaths,
                               std::size_t numSteps,
                               T h);

    std::size_t stride() const { return _stride; }
    void enablePathStorage(bool storePaths) { _storePaths = storePaths; }

private:
    // Shared implementation for terminal pricing (consumes a Z provider)
    template <typename ZProvider>
    MCReturn<T> priceTerminal_impl_(const std::string& optionType,
                                    std::size_t numPaths,
                                    ZProvider&& getZ);

private:
    Opt&  _option;
    Step  _stepper;

    std::size_t _numPaths = 0;

    bool _storePaths = false;
    bool _useAV = false;
    bool _useCV = false;

    std::vector<T> _pricePaths_flat;
    std::size_t _stride = 0;
};


// ========================== IMPLEMENTATION ==========================

// -------- Terminal pricing: shared implementation --------
template <typename T, typename Opt, typename Step>
template <typename ZProvider>
MCReturn<T>
MonteCarlo<T, Opt, Step>::priceTerminal_impl_(const std::string& optionType,
                                              std::size_t numPaths,
                                              ZProvider&& getZ)
{
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");
    if (numPaths <= 1)
        throw std::invalid_argument("Number of paths must be >= 2 for priceTerminal().");

    _option.setType(optionType);
    setNumPaths(numPaths);

    const auto g   = _option.payOff();
    const T S0     = _option.getSpotPrice();
    const T r      = _option.getRiskFreeRate();
    const T sigma  = _option.getVolatility();
    const T Tmat   = _option.getMaturityTime();
    const T disc   = std::exp(-r * Tmat);

    if (!_useCV && !_useAV) {
        long N = 0;
        T sumY=0, sumY2=0;

        #if defined(USE_OMP)
        #pragma omp parallel for schedule(static) reduction(+:sumY,sumY2)
        #endif
        for (std::size_t i = 0; i < numPaths; ++i) {
            const T Z  = getZ(i);
            const T ST = _stepper.terminalDrawWithZ(S0, r, sigma, Tmat, Z);
            const T Y  = disc * g(ST);
            sumY += Y; sumY2 += Y*Y;
        }

        N = static_cast<long>(numPaths);
        const T Nf   = static_cast<T>(N);
        const T EY   = sumY / Nf;
        const T varY = (sumY2 - Nf * EY * EY) / (Nf - T(1));
        return { EY, varY, std::sqrt(varY / Nf) };
    }

    if (!_useAV) {
        long N = 0; T sumY=0, sumY2=0, sumX=0, sumX2=0, sumXY=0;

        #if defined(USE_OMP)
        #pragma omp parallel for schedule(static) reduction(+:sumY,sumY2,sumX,sumX2,sumXY)
        #endif
        for (std::size_t i = 0; i < numPaths; ++i) {
            const T Z  = getZ(i);
            const T ST = _stepper.terminalDrawWithZ(S0, r, sigma, Tmat, Z);

            const T Y = disc * g(ST);
            sumY  += Y;     sumY2 += Y*Y;

            if (_useCV) {
                const T X = disc * ST;  // E[X]=S0 exactly
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

    // ---------- AV branch ----------
    if (numPaths % 2 != 0)
        throw std::invalid_argument("priceTerminal with AV=on requires even numPaths.");
    const std::size_t pairs = numPaths / 2;
    long Npairs = static_cast<long>(pairs);

    T sumA=0, sumA2=0, sumXb=0, sumXb2=0, sumAXb=0;

    #if defined(USE_OMP)
    #pragma omp parallel for schedule(static) reduction(+:sumA,sumA2,sumXb,sumXb2,sumAXb)
    #endif
    for (std::size_t i = 0; i < pairs; ++i) {
        const T Z = getZ(i);  // one normal per pair
        const T ST1 = _stepper.terminalDrawWithZ(S0, r, sigma, Tmat, +Z);
        const T ST2 = _stepper.terminalDrawWithZ(S0, r, sigma, Tmat, -Z);

        const T Y1 = disc * g(ST1);
        const T Y2 = disc * g(ST2);
        const T A  = T(0.5) * (Y1 + Y2);

        sumA  += A; sumA2 += A*A;

        if (_useCV) {
            const T X1 = disc * ST1;
            const T X2 = disc * ST2;
            const T Xb = T(0.5) * (X1 + X2);
            sumXb  += Xb; sumXb2 += Xb*Xb; sumAXb += A*Xb;
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

    const T EA_adj  = EA - beta * (EXb - S0);
    const T var_adj = varA - ((varXb > T(0)) ? (covAX * covAX / varXb) : T(0));
    return { EA_adj, var_adj, std::sqrt(var_adj / Nf) };
}

// -------- Terminal pricing: legacy overload (generate inside) --------
template <typename T, typename Opt, typename Step>
MCReturn<T>
MonteCarlo<T, Opt, Step>::priceTerminal(const std::string& optionType, std::size_t numPaths)
{
    // The provider ignores i and just returns a fresh normal each time.
    auto getZ = [&](std::size_t /*i*/) -> T {
        return static_cast<T>(standard_normal_sample());
    };
    return priceTerminal_impl_(optionType, numPaths, getZ);
}

// -------- Terminal pricing: pre-generated normals via span --------
template <typename T, typename Opt, typename Step>
MCReturn<T>
MonteCarlo<T, Opt, Step>::priceTerminal(const std::string& optionType,
                                        std::size_t numPaths,
                                        std::span<const T> preGeneratedZ)
{
    if (!_useAV) {
        if (preGeneratedZ.size() < numPaths)
            throw std::invalid_argument("priceTerminal: preGeneratedZ size must be >= numPaths when AV=off.");
        auto getZ = [&](std::size_t i) -> T { return preGeneratedZ[i]; };
        return priceTerminal_impl_(optionType, numPaths, getZ);
    } else {
        if (numPaths % 2 != 0)
            throw std::invalid_argument("priceTerminal with AV=on requires even numPaths.");
        const std::size_t pairs = numPaths / 2;
        if (preGeneratedZ.size() < pairs)
            throw std::invalid_argument("priceTerminal: preGeneratedZ size must be >= numPaths/2 when AV=on.");
        auto getZ = [&](std::size_t i) -> T { return preGeneratedZ[i]; }; // one per pair
        return priceTerminal_impl_(optionType, numPaths, getZ);
    }
}

// =================== impl ===================

// Compute greeks with Finite Differences (forward, backward and central)
// Uses pre-generated normals (6 independent sets) for strong CRN per estimator.
template <typename T, typename Opt, typename Step>
MCgreeksFD MonteCarlo<T, Opt, Step>::computeGreeksFD(const std::string& optionType,
                                                     std::size_t numPaths,
                                                     std::size_t /*numSteps*/,
                                                     T h)
{
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");
    if (numPaths <= 1) throw std::invalid_argument("numPaths >= 2");
    if (!(h > T(0)))   throw std::invalid_argument("step h must be > 0");

    _option.setType(optionType);

    // Save & restore to avoid side-effects
    Opt opt_backup = _option;

    const T S0  = _option.getSpotPrice();
    const T sig = _option.getVolatility();

    // Number of normals needed for one priceTerminal call
    const std::size_t Nz = _useAV ? (numPaths / 2) : numPaths;

    // ===================== DELTA (Forward) =====================
    const auto Z_deltaFD = generateStandardNormalSamples<T>(Nz);
    _option.setSpotPrice(S0 + h);
    auto Rph = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_deltaFD));

    _option.setSpotPrice(S0);
    auto Rp  = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_deltaFD));

    const T deltaFD = (Rph.price - Rp.price) / h;

    // ===================== DELTA (Backward) =====================
    const auto Z_deltaBD = generateStandardNormalSamples<T>(Nz);
    _option.setSpotPrice(S0);
    Rp  = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_deltaBD));

    _option.setSpotPrice(std::max(S0 - h, T(1e-12)));
    auto Rmh = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_deltaBD));

    const T deltaBD = (Rp.price - Rmh.price) / h;

    // ===================== DELTA (Central) =====================
    const auto Z_deltaCD = generateStandardNormalSamples<T>(Nz);
    _option.setSpotPrice(S0 + h);
    Rph = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_deltaCD));

    _option.setSpotPrice(std::max(S0 - h, T(1e-12)));
    Rmh = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_deltaCD));

    const T deltaCD = (Rph.price - Rmh.price) / (T(2) * h);

    // Reset S0 for vega
    _option.setSpotPrice(S0);

    // ===================== VEGA (Forward) =====================
    const auto Z_vegaFD = generateStandardNormalSamples<T>(Nz);
    _option.setVolatility(sig + h);
    auto Rvh = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_vegaFD));

    _option.setVolatility(sig);
    auto Rv  = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_vegaFD));

    const T vegaFD = (Rvh.price - Rv.price) / h;

    // ===================== VEGA (Backward) =====================
    const auto Z_vegaBD = generateStandardNormalSamples<T>(Nz);
    _option.setVolatility(sig);
    Rv  = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_vegaBD));

    _option.setVolatility(std::max(sig - h, T(1e-8)));
    auto Rvm = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_vegaBD));

    const T vegaBD = (Rv.price - Rvm.price) / h;

    // ===================== VEGA (Central) =====================
    const auto Z_vegaCD = generateStandardNormalSamples<T>(Nz);
    _option.setVolatility(sig + h);
    Rvh = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_vegaCD));

    _option.setVolatility(std::max(sig - h, T(1e-8)));
    Rvm = this->priceTerminal(optionType, numPaths, std::span<const T>(Z_vegaCD));

    const T vegaCD = (Rvh.price - Rvm.price) / (T(2) * h);

    // Restore original option state
    _option = opt_backup;

    return MCgreeksFD{ deltaFD, vegaFD, deltaBD, vegaBD, deltaCD, vegaCD };
}

// Strong error on STATE vs exact GBM (CRN)
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

    T sumAbs = 0, sumAbs2 = 0, sumSq = 0;

    #if defined(USE_OMP)
        #pragma omp parallel for schedule(static) reduction(+:sumAbs,sumSq,sumAbs2)
    #endif

    for (std::size_t i = 0; i < numPaths; ++i) {
        T S  = S0;
        T WT = 0;

        for (std::size_t j = 0; j < numSteps; ++j) {
            const T Z = static_cast<T>(standard_normal_sample());
            S  = _stepper.advanceWithZ(S, r, sig, dt, Z);
            WT += std::sqrt(dt) * Z;
        }
        const T Sexact = S0 * std::exp((r - T(0.5)*sig*sig)*Tmat + sig*WT);

        const T diff  = Sexact - S;
        const T adiff = std::fabs(diff);

        sumAbs  += adiff;
        sumAbs2 += adiff*adiff;
        sumSq   += diff*diff;
    }

    const T Nf     = static_cast<T>(numPaths);
    const T meanAbs= sumAbs / Nf;
    const T varAbs = (sumAbs2 - Nf * meanAbs * meanAbs) / (Nf - T(1));
    const T seAbs  = std::sqrt( varAbs / Nf );
    const T rms    = std::sqrt( sumSq / Nf );

    return { static_cast<double>(meanAbs),
             static_cast<double>(rms),
             static_cast<double>(seAbs) };
}

// Weak error on PAYOFF (bias) vs exact GBM (CRN)
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

    T sumYh = 0, sumYh2 = 0;
    T sumYs = 0, sumYs2 = 0;
    T sumD  = 0, sumD2  = 0;

    #if defined(USE_OMP)
        #pragma omp parallel for schedule(static) reduction(+:sumYh,sumYh2,sumYs,sumYs2,sumD,sumD2)
    #endif

    for (std::size_t i = 0; i < numPaths; ++i) {
        T S  = S0;
        T WT = 0;

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

// Statistical convergence via pricePath (stepper-agnostic)
template <typename T, typename Opt, typename Step>
MCStatConvergence
MonteCarlo<T, Opt, Step>::statisticalConvergenceByPaths(const std::string& optionType,
                                                        const std::vector<std::size_t>& Ns,
                                                        bool useAV, bool useCV)
{
    this->useAntitheticVariates(useAV);
    this->useControlVariates(useCV);

    MCStatConvergence out;
    out.Ns      = Ns;
    out.prices  = std::vector<double>(Ns.size());
    out.stdErrs = std::vector<double>(Ns.size());

    std::vector<double> xN; xN.reserve(Ns.size());
    std::vector<double> se; se.reserve(Ns.size());

    for (std::size_t i = 0; i < Ns.size(); ++i) {
        const std::size_t N = Ns[i];
        auto r = this->pricePath(optionType, N);
        out.prices[i]  = static_cast<double>(r.price);
        out.stdErrs[i] = static_cast<double>(r.stdDev);
        xN.push_back(static_cast<double>(N));
        se.push_back(out.stdErrs[i]);
    }
    out.slope_loglog = (Ns.size() >= 2)
                 ? slope_loglog_xy(xN, se)
                 : std::numeric_limits<double>::quiet_NaN();
    return out;
}

// Full-path (discretized) MC
template <typename T, typename Opt, typename Step>
MCReturn<T> MonteCarlo<T, Opt, Step>::pricePath(const std::string& optionType, 
                                                std::size_t numPaths, 
                                                std::size_t numSteps,
                                                bool storePaths)
{
    if (optionType != "call" && optionType != "put")
        throw std::invalid_argument("Invalid option type");
    if (numPaths <= 1)
        throw std::invalid_argument("Number of paths must be >= 2");

    _option.setType(optionType);
    setNumPaths(numPaths);

    const auto g   = _option.payOff();
    const T S0     = _option.getSpotPrice();
    const T r      = _option.getRiskFreeRate();
    const T sigma  = _option.getVolatility();
    const T Tmat   = _option.getMaturityTime();
    const T disc   = std::exp(-r * Tmat);
    const T dt = Tmat / static_cast<T>(numSteps);
    // const std::size_t numSteps = std::max<std::size_t>(1, static_cast<std::size_t>(TRADING_DAYS_PER_YEAR * Tmat));

    // final decision to store or not (call-site param OR persistent flag)
    const bool store = storePaths || _storePaths;

    if (store) {
        _stride = numSteps + 1;
        _pricePaths_flat.assign(numPaths * _stride, T(0));
    } else {
        _stride = 0;
        _pricePaths_flat.clear();
    }

    if (!_useAV) {
        long N = 0; T sumY=0, sumY2=0, sumX=0, sumX2=0, sumXY=0;

        #if defined(USE_OMP)
            #pragma omp parallel for schedule(static) reduction(+:sumY,sumY2,sumX,sumX2,sumXY)
        #endif

        for (std::size_t i = 0; i < numPaths; ++i) {
            T S = S0, WT = 0;
            // const std::size_t base = i * _stride;
            std::size_t base = 0;
            if (store) {
                base = i * _stride;
                _pricePaths_flat[base] = S;
            }

            for (std::size_t j = 1; j <= numSteps; ++j) {
                const T Z = static_cast<T>(standard_normal_sample());
                S  = _stepper.advanceWithZ(S, r, sigma, dt, Z);
                WT += std::sqrt(dt) * Z;
                if (store) _pricePaths_flat[base + j] = S;
            }
            const T Y = disc * g(S);
            sumY += Y; sumY2 += Y*Y;

            if (_useCV) {
                const T ST_exact = S0 * std::exp((r - T(0.5)*sigma*sigma)*Tmat + sigma*WT);
                const T X = disc * ST_exact; // E[X] = S0
                sumX  += X;  sumX2 += X*X;  sumXY += Y*X;
            }
        }
        N = static_cast<long>(numPaths);
        const T Nf = static_cast<T>(N);
        const T EY = sumY / Nf;
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

    // AV: pairs of antithetic legs
    const std::size_t pairs = (numPaths / 2);
    if (pairs == 0) throw std::invalid_argument("Need at least 2 paths when AV is on.");
    long Npairs = static_cast<long>(pairs);

    T sumA=0, sumA2=0, sumXb=0, sumXb2=0, sumAXb=0;

    #if defined(USE_OMP)
        #pragma omp parallel for schedule(static) reduction(+:sumA,sumA2,sumXb,sumXb2,sumAXb)
    #endif

    for (std::size_t i = 0; i < pairs; ++i) {
        T S1 = S0, WT1 = 0;
        T S2 = S0, WT2 = 0;

        // const std::size_t base1 = (2*i) * _stride;
        // _pricePaths_flat[base1] = S1;
        // const std::size_t base2 = (2*i + 1) * _stride;
        // _pricePaths_flat[base2] = S2;

        std::size_t base1 = 0, base2 = 0;
        if (store) {
            base1 = (2*i) * _stride;
            _pricePaths_flat[base1] = S1;
            base2 = (2*i + 1) * _stride;
            _pricePaths_flat[base2] = S2;
        }

        for (std::size_t j = 1; j <= numSteps; ++j) {
            const T Z = static_cast<T>(standard_normal_sample());
            S1  = _stepper.advanceWithZ(S1, r, sigma, dt, +Z);
            S2  = _stepper.advanceWithZ(S2, r, sigma, dt, -Z);
            WT1 += std::sqrt(dt) * (+Z);
            WT2 += std::sqrt(dt) * (-Z);
            // _pricePaths_flat[base1 + j] = S1;
            // _pricePaths_flat[base2 + j] = S2;
            if (store) {
                _pricePaths_flat[base1 + j] = S1;
                _pricePaths_flat[base2 + j] = S2;
            }
        }

        const T Y1 = disc * g(S1);
        const T Y2 = disc * g(S2);
        const T A  = T(0.5) * (Y1 + Y2);

        sumA  += A; sumA2 += A*A;

        if (_useCV) {
            const T ST1_exact = S0 * std::exp((r - T(0.5)*sigma*sigma)*Tmat + sigma*WT1);
            const T ST2_exact = S0 * std::exp((r - T(0.5)*sigma*sigma)*Tmat + sigma*WT2);
            const T X1 = disc * ST1_exact;
            const T X2 = disc * ST2_exact;
            const T Xb = T(0.5) * (X1 + X2);
            sumXb  += Xb; sumXb2 += Xb*Xb; sumAXb += A*Xb;
        }
    }

    const T Nf = static_cast<T>(Npairs);
    const T EA = sumA / Nf;
    const T varA = (sumA2 - Nf * EA * EA) / (Nf - T(1));

    if (!_useCV) return { EA, varA, std::sqrt(varA / Nf) };

    const T EXb   = sumXb / Nf;
    const T varXb = (sumXb2 - Nf * EXb * EXb) / (Nf - T(1));
    const T covAX = (sumAXb - Nf * EA * EXb) / (Nf - T(1));
    const T beta  = (varXb > T(0)) ? (covAX / varXb) : T(0);

    const T EA_adj  = EA - beta * (EXb - S0);
    const T var_adj = varA - ((varXb > T(0)) ? (covAX * covAX / varXb) : T(0));
    return { EA_adj, var_adj, std::sqrt(var_adj / Nf) };
}


#endif // MC_ENGINE_HPP

