#ifndef CONVERGENCE_HPP
#define CONVERGENCE_HPP

#include <vector>
#include <string>
#include <limits>
#include "mc_types.hpp"
#include "mc_engine.hpp"
#include "math.hpp"

template <typename T, typename Opt, typename Step>
class ConvergenceAnalyzer {
public:
    ConvergenceAnalyzer(Opt& opt, const Step& stepper) : _opt(opt), _stepper(stepper) {}

    MCStatConvergence mcConvergence(const std::string& optionType,
                                    const std::vector<std::size_t>& Ns,
                                    bool useAV, bool useCV) const
    {
        Opt  optCopy = _opt;
        MonteCarlo<T, Opt, Step> mc(optCopy, _stepper);
        return mc.statisticalConvergenceByPaths(optionType, Ns, useAV, useCV);
    }

    StrongConvergenceResult strongConvergence(std::size_t numPaths,
                                              const std::vector<std::size_t>& Ms) const
    {
        Opt  optCopy = _opt;
        MonteCarlo<T, Opt, Step> mc(optCopy, _stepper);

        StrongConvergenceResult R;
        R.M = Ms;
        R.dt.resize(Ms.size());
        R.meanAbs.resize(Ms.size());
        R.rms.resize(Ms.size());

        const T Tmat = _opt.getMaturityTime();

        for (std::size_t i=0;i<Ms.size();++i) {
            const std::size_t M = Ms[i];
            const double dt = static_cast<double>(Tmat) / double(M);
            auto s = mc.strongErrorOnState(numPaths, M);
            R.dt[i] = dt;
            R.meanAbs[i] = s.meanAbs;
            R.rms[i]     = s.rms;
        }
        
        // compute slopes only when we have at least two points
        if (R.M.size() >= 2) {
            R.slope_meanAbs = slope_loglog_xy(R.dt, R.meanAbs);
            R.slope_rms     = slope_loglog_xy(R.dt, R.rms);
        } else {
            R.slope_meanAbs = std::numeric_limits<double>::quiet_NaN();
            R.slope_rms     = std::numeric_limits<double>::quiet_NaN();
        }
        return R;
    }

    WeakConvergenceResult weakConvergence(const std::string& optionType,
                                          std::size_t numPaths,
                                          const std::vector<std::size_t>& Ms) const
    {
        Opt  optCopy = _opt;
        MonteCarlo<T, Opt, Step> mc(optCopy, _stepper);

        WeakConvergenceResult R;
        R.M = Ms;
        R.dt.resize(Ms.size());
        R.bias.resize(Ms.size());
        R.seBias.resize(Ms.size());

        const T Tmat = _opt.getMaturityTime();

        for (std::size_t i=0;i<Ms.size();++i) {
            const std::size_t M = Ms[i];
            const double dt = static_cast<double>(Tmat) / double(M);
            auto w = mc.weakErrorOnPayoff(optionType, numPaths, M);
            R.dt[i]    = dt;
            R.bias[i]  = w.bias;
            R.seBias[i]= w.seBias;
        }
        
        // compute slope only when we have at least two points
        if (R.M.size() >= 2) {
            R.slope_bias = slope_loglog_xy(R.dt, R.bias);
        } else {
            R.slope_bias = std::numeric_limits<double>::quiet_NaN();
        }
        return R;

    }

private:
    Opt  _opt;
    Step _stepper;
};

#endif // CONVERGENCE_HPP