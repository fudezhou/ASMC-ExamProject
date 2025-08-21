#ifndef STEPPER_HPP
#define STEPPER_HPP

#include <cmath>
#include "math.hpp"

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

// GBM exact
template <typename T>
class GBM: public Stepper<T> {
public:
    T advance(T S_t, T r, T sigma, T dt) const {
        const T Z = static_cast<T>(standard_normal_sample());
        return S_t * std::exp((r - 0.5*sigma*sigma)*dt + sigma*Z*std::sqrt(dt));
    }
    T advanceWithZ(T S_t, T r, T sigma, T dt, T Z) const {
        return S_t * std::exp((r - 0.5*sigma*sigma)*dt + sigma*Z*std::sqrt(dt));
    }
    T terminalDraw(T S0, T r, T sigma, T Tmat) const {
        const T Z = static_cast<T>(standard_normal_sample());
        return S0 * std::exp((r - 0.5*sigma*sigma)*Tmat + sigma*Z*std::sqrt(Tmat));
    }
    T terminalDrawWithZ(T S0, T r, T sigma, T Tmat, T Z) const {
        return S0 * std::exp((r - 0.5*sigma*sigma)*Tmat + sigma*Z*std::sqrt(Tmat));
    }

    virtual const char *getName() const override { return "GBM"; };
};

// Euler–Maruyama
template <typename T>
class Euler: public Stepper<T> {
public:
    T advance(T S_t, T r, T sigma, T dt) const {
        const T Z = static_cast<T>(standard_normal_sample());
        return S_t * (T(1) + r*dt + sigma*Z*std::sqrt(dt));
    }
    T advanceWithZ(T S_t, T r, T sigma, T dt, T Z) const {
        return S_t * (T(1) + r*dt + sigma*Z*std::sqrt(dt));
    }

    
    virtual const char *getName() const override { return "Euler"; };
};

// Milstein
template <typename T>
class Milstein: public Stepper<T> {
public:
    T advance(T S_t, T r, T sigma, T dt) const {
        const T Z = static_cast<T>(standard_normal_sample());
        return  S_t * (T(1) + r*dt + sigma*Z*std::sqrt(dt)
            + T(0.5)*sigma*sigma*(Z*Z - T(1))*dt);
    }
    T advanceWithZ(T S_t, T r, T sigma, T dt, T Z) const {
        return  S_t * (T(1) + r*dt + sigma*Z*std::sqrt(dt)
            + T(0.5)*sigma*sigma*(Z*Z - T(1))*dt);
    }

    virtual const char *getName() const override { return "Milstein"; };
};

#endif // STEPPER_HPP


