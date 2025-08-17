#include <iostream>
// #include <iomanip> // for std::setprecision
#include "../include/utilities.hpp"

int main() {

    // std::cout << "CDF at 0.15: " << norm_cdf(0.15) << std::endl;

    double spotPrice = 100.0;
    double strikePrice = 100.0;
    double riskFreeRate = 0.05;
    double volatility = 0.2;
    double timeToExpiration = 1.0;

    BlackAndScholes<double> bs(spotPrice, strikePrice, riskFreeRate, volatility, timeToExpiration);

    double callPrice = bs.price("call");
    double putPrice = bs.price("put");
    //double error = bs.price("error") - callPrice;

    std::cout << "Call Price: " << callPrice << std::endl;
    std::cout << "Put Price: " << putPrice << std::endl;
    //std::cout << "Error: " << error << std::endl;

    return 0;
}