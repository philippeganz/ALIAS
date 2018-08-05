///
/// \file src/test/WS.cpp
/// \brief Test suite to test the WS method.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.5.0
/// \date 2018-07-21
/// \copyright GPL-3.0
///

#include "test/WS.hpp"

namespace astroqut{
namespace test{
namespace astro{

void Chandra()
{
    std::cout << std::endl << std::endl << "WS test with 512x512 Chandra data." << std::endl;

    Matrix<double> picture(std::string("data/512_chandra/F.data"), 262144, 1, double());

    Matrix<double> sensitivity(std::string("data/512_chandra/E.data"), 262144, 1, double());

    Matrix<double> background(std::string("data/512_chandra/O.data"), 262144, 1, double());

    Matrix<double> expected_result(std::string("data/512_chandra/MATLAB/1000/sol.data"), 263168, 1, double());

    WS::Parameters<double> options;
    options.MC_max = 5000;
    options.MC_quantile_PF = (size_t) (options.MC_max * (1.0 - 1.0/(std::sqrt(PI*std::log(512)))) - 1.0);
    options.MC_quantile_PS = (size_t) (options.MC_max * (1.0 - 1.0/(512*512)) - 1.0);

    Matrix<double> solution = WS::Solve(picture, sensitivity, background, options);

    std::cout << std::endl;
    std::cout << "MATLAB result : " << std::endl;
    std::cout << "FISTA: did not converge in 1000 iterations (nnz(x)= 2329)" << std::endl;
    std::cout << "FISTA: Relative error: -1.93e-06" << std::endl;
    std::cout << "FISTA: Final FLasso value: -1.387380e+06" << std::endl;
    std::cout << "Elapsed time is 9293.529344 seconds." << std::endl << std::endl;

    Compare(expected_result, solution);

    double relative_error = std::abs((solution - expected_result).Norm(two)) / std::abs(expected_result.Norm(two));
    std::cout << relative_error << " relative norm error." << std::endl;
}

} // namespace astro
} // namespace test
} // namespace astroqut
