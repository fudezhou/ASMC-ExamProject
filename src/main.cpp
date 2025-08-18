#include <iostream>
// #include <iomanip> // for std::setprecision
#include "../include/utilities.hpp"

int main() {

    EuropeanOption<double> europeanOption;
    europeanOption.setSpotPrice(100.0);
    europeanOption.setStrikePrice(100.0);
    europeanOption.setRiskFreeRate(0.05);
    europeanOption.setVolatility(0.2);
    europeanOption.setMaturityTime(1.0);

    europeanOption.info();

    std::cout << "\nEuropean Option Price (Call): " << europeanOption.closedForm("call") << "\n";
    std::cout << "European Option Price (Put): " << europeanOption.closedForm("put") << "\n";

    // npaths default = 100000
    MCReturn nMCGBMcall = europeanOption.naiveMonteCarloGBM("call", 100000);
    std::cout << "European Option Price \n(Call - Naive Monte Carlo GBM): " << nMCGBMcall.price << " ± " << 1.96 * nMCGBMcall.stdDev << " CI at 95%\n";
    MCReturn nMCGBMput = europeanOption.naiveMonteCarloGBM("put", 100000);
    std::cout << "European Option Price \n(Put - Naive Monte Carlo GBM): " << nMCGBMput.price << " ± " << 1.96 * nMCGBMput.stdDev << " CI at 95%\n";

    MCReturn fLVCVMCGBMcall = europeanOption.fastLowVarMCGBM("call", 20000);
    std::cout << "European Option Price \n(Call - Fast Low Variance Monte Carlo GBM): " << fLVCVMCGBMcall.price << " ± " << 1.96 * fLVCVMCGBMcall.stdDev << " CI at 95%\n";
    MCReturn fLVCVMCGBMput = europeanOption.fastLowVarMCGBM("put", 20000);
    std::cout << "European Option Price \n(Put - Fast Low Variance Monte Carlo GBM): " << fLVCVMCGBMput.price << " ± " << 1.96 * fLVCVMCGBMput.stdDev << " CI at 95%\n";

    return 0;
}