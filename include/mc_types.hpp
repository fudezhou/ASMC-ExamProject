#ifndef MC_TYPES_HPP
#define MC_TYPES_HPP

#include <vector>
#include <cstddef>

template <typename T>
struct MCReturn {
    T price;
    T var;
    T stdDev; // standard error of the estimator (used for CI)
};

struct StrongStats {
    double meanAbs; // E[|ΔS_T|]
    double rms;     // sqrt(E[(ΔS_T)^2])
    double seAbs;   // SE of |ΔS_T|
};

struct WeakStats {
    double bias;       // |E[Y_h] - E[Y]|
    double seBias;     // SE of the difference
    double meanApprox; // E[Y_h]
    double meanExact;  // E[Y]
};

struct MCStatConvergence {
    std::vector<std::size_t> Ns;
    std::vector<double> prices;
    std::vector<double> stdErrs;
    double slope_loglog = 0.0; // ~ -0.5
};

struct StrongConvergenceResult {
    std::vector<std::size_t> M;
    std::vector<double> dt;
    std::vector<double> meanAbs;
    std::vector<double> rms;
    double slope_meanAbs = 0.0;
    double slope_rms     = 0.0;
};

struct WeakConvergenceResult {
    std::vector<std::size_t> M;
    std::vector<double> dt;
    std::vector<double> bias;
    std::vector<double> seBias;
    double slope_bias = 0.0;
};

#endif // MC_TYPES_HPP

