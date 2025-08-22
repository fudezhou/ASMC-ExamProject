#ifndef MATH_HPP
#define MATH_HPP

#include <cmath>
#include <random>
#include <vector>
#include <stdexcept>
#include <cstddef>

// Trading calendar used in your full-path discretization
constexpr std::size_t TRADING_DAYS_PER_YEAR = 252;

// ===== RNG utilities (deterministic, OpenMP-safe) =====
#include <atomic>
#if defined(USE_OMP)
  #include <omp.h>
#endif

// SplitMix64 — cheap mixer to expand seeds
inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Global base seed (default fixed for reproducibility)
inline uint64_t& rng_base_seed_ref() {
    static uint64_t s = 0xD1B54A32D192ED03ULL; // any constant you like
    return s;
}

// Allow user code to set the base seed (call once at program start)
inline void rng_set_seed(uint64_t s) { rng_base_seed_ref() = s; }

// Compute a per-thread, per-stream seed deterministically
inline uint64_t rng_compute_seed(uint64_t stream_id = 0) {
#if defined(USE_OMP)
    unsigned tid = static_cast<unsigned>(omp_get_thread_num());
#else
    unsigned tid = 0u;
#endif
    uint64_t s = rng_base_seed_ref();
    s ^= splitmix64(stream_id + 0x9E3779B97F4A7C15ULL);
    s ^= splitmix64(static_cast<uint64_t>(tid) + 1u);
    return splitmix64(s);
}

// Thread-local engine & normal dist, with lazy deterministic seeding
inline std::mt19937_64& tls_engine(uint64_t stream_id = 0) {
    struct TL { std::mt19937_64 eng; bool seeded = false; };
    thread_local TL tl;
    if (!tl.seeded) { tl.eng.seed(rng_compute_seed(stream_id)); tl.seeded = true; }
    return tl.eng;
}

// Optional: force reseed of the current thread (e.g., before a new “sequence”)
inline void rng_reseed_thread(uint64_t stream_id = 0) {
    tls_engine().seed(rng_compute_seed(stream_id));
}

// -------- Standard Normal draw --------
inline double standard_normal_sample(uint64_t stream_id = 0) {
    thread_local std::normal_distribution<double> dist(0.0, 1.0);
    return dist(tls_engine(stream_id));
}


// -------- Normal CDF --------
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// -------- Normal PDF --------
inline double norm_pdf(double x) {
    return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI);
}

// // -------- Thread-local RNG & N(0,1) --------
// inline double standard_normal_sample() {
//     thread_local std::mt19937_64 rng{
//         0x9E3779B97F4A7C15ULL ^ (std::random_device{}() + (uint64_t)std::hash<void*>{}(&rng))
//     };
//     thread_local std::normal_distribution<double> dist(0.0, 1.0);
//     return dist(rng);
// }

// -------- slope of log(y) vs log(x) --------
inline double slope_loglog_xy(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) throw std::invalid_argument("Need >=2 points.");
    double Sx=0, Sy=0, Sxx=0, Sxy=0; const double n = (double)x.size();
    for (std::size_t i=0;i<x.size();++i) {
        const double lx = std::log(x[i]);
        const double ly = std::log(y[i]);
        Sx  += lx; Sy  += ly;
        Sxx += lx*lx; Sxy += lx*ly;
    }
    const double num = n*Sxy - Sx*Sy;
    const double den = n*Sxx - Sx*Sx;
    return num / den;
}

template <typename T>
inline std::vector<T> generateStandardNormalSamples(std::size_t numSamples) {
    std::vector<T> samples(numSamples);
    for (std::size_t i = 0; i < numSamples; ++i) {
        samples[i] = static_cast<T>(standard_normal_sample());
    }
    return samples;
}

// inline double slope_loglog_xy(const std::vector<double>& x,
//                               const std::vector<double>& y)
// {
//     const int n = static_cast<int>(x.size());
//     double Sx=0, Sy=0, Sxx=0, Sxy=0;
//     for (int i=0;i<n;++i) {
//         const double lx = std::log(x[i]);
//         const double ly = std::log(y[i]);
//         Sx+=lx; Sy+=ly; Sxx+=lx*lx; Sxy+=lx*ly;
//     }
//     return (n*Sxy - Sx*Sy) / (n*Sxx - Sx*Sx);
// }

#endif // MATH_HPP
