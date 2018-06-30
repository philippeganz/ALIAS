///
/// \file src/test/WS.cpp
/// \brief Test suite to test the WS method.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.4.1
/// \date 2018-06-16
/// \copyright GPL-3.0
///

#include "test/WS.hpp"

namespace astroqut{
namespace test{
namespace astro{

void Chandra()
{
    std::cout << "WS test with 512x512 Chandra data." << std::endl;

    Matrix<double> picture(std::string("data/512_chandra/F.data"), 262144, 1, double());

    Matrix<double> sensitivity(std::string("data/512_chandra/E.data"), 262144, 1, double());

    Matrix<double> background(std::string("data/512_chandra/O.data"), 262144, 1, double());

    Matrix<double> expected_result(std::string("data/512_chandra/solstatic.data"), 263168, 1, double());

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    Matrix<double> solution = WS::Solve(picture, sensitivity, background, WS::Parameters<double>());

    end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_time = end-start;
    std::cout << "Time for FISTA solver with Chandra data :";
    std::cout << std::defaultfloat << elapsed_time.count() << " seconds" << std::endl << std::endl;

    std::cout << "MATLAB result : " << std::endl;
    std::cout << "FISTA: converged in 1000 iterations (nnz(x)= 2984)" << std::endl;
    std::cout << "FISTA: Relative error: -1.34e-05" << std::endl;
    std::cout << "Elapsed time is 7952.558757 seconds." << std::endl << std::endl;

    Compare(expected_result, solution);

    double relative_error = std::abs((solution - expected_result).Norm(two)) / std::abs(expected_result.Norm(two));
    std::cout << relative_error << " relative norm error." << std::endl;

    "data/512_chandra/computedsol.data" << solution;
}

} // namespace astro
} // namespace test
} // namespace astroqut
