#include "../include/math.hpp"
#include "../include/mc_types.hpp"
#include "../include/option.hpp"
#include "../include/stepper.hpp"
#include "../include/mc_engine.hpp"
#include "../include/convergence.hpp"

// ===== Explicit instantiations for your concrete types =====
template class Option<double>;
template class EuropeanOption<double>;
template class GBM<double>;
template class Euler<double>;
template class Milstein<double>;

template class MonteCarlo<double, EuropeanOption<double>, GBM<double>>;
template class MonteCarlo<double, EuropeanOption<double>, Euler<double>>;
template class MonteCarlo<double, EuropeanOption<double>, Milstein<double>>;

template class ConvergenceAnalyzer<double, EuropeanOption<double>, GBM<double>>;
template class ConvergenceAnalyzer<double, EuropeanOption<double>, Euler<double>>;
template class ConvergenceAnalyzer<double, EuropeanOption<double>, Milstein<double>>;
