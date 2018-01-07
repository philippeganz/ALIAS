///
/// \file src/test/fista.cpp
/// \brief Implementation of the FISTA class test suite.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.2.0
/// \date 2017-12-28
/// \copyright GPL-3.0
///

#include "test/fista.hpp"

namespace astroqut{
namespace test{
namespace fista{

bool SmallExample1()
{
    std::cout << "FISTA test with small data : " << std::endl << std::endl;

    astroqut::MatMult<double> A(Matrix(new double[12]{1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0}, 3, 4), 3, 4);
    astroqut::Matrix<double> u(new double[3]{3.0,2.0,1.0}, 3, 1);
    astroqut::Matrix<double> b(new double[3]{1.0,1.0,2.0}, 3, 1);
    astroqut::fista::poisson::Parameters options;
    options.tol = 1e-4;
    options.init_value = astroqut::Matrix(new double[4]{}, 4, 1);
    options.log_period = 50;

    astroqut::Matrix<double> expected_result(new double[4]{0.939214191176937, 0.0, 0.0, -0.650697826008454}, 4, 1);
    astroqut::Matrix<double> actual_result = astroqut::fista::poisson::Solve(A, u, b, 1, options);

    expected_result.Print();
    actual_result.Print();

    double relative_error = std::abs((actual_result - expected_result).Norm(two)) / std::abs(expected_result.Norm(two));

    bool fista_test = (relative_error < 100*options.tol);

    std::cout << std::endl << (fista_test ? "Success" : "Failure") << ", achieved ";
    std::cout << relative_error << " relative norm error with a tol of " << options.tol << "." << std::endl;

    return fista_test;
}

bool SmallExample2()
{
    std::cout << "FISTA test with small data : " << std::endl << std::endl;

    astroqut::MatMult<double> A(Matrix(new double[12]{1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0}, 3, 4), 3, 4);
    astroqut::Matrix<double> u(new double[3]{3.0,2.0,1.0}, 3, 1);
    astroqut::Matrix<double> b(new double[3]{1.0,1.0,2.0}, 3, 1);
    astroqut::fista::poisson::Parameters options;
    options.tol = 1e-8;
    options.init_value = astroqut::Matrix(new double[4]{}, 4, 1);
    options.log_period = 50;

    astroqut::Matrix<double> expected_result(new double[4]{0.973320232506235, 0.0, 0.0, -0.674615621567228}, 4, 1);
    astroqut::Matrix<double> actual_result = astroqut::fista::poisson::Solve(A, u, b, 1, options);

    expected_result.Print();
    actual_result.Print();

    double relative_error = std::abs((actual_result - expected_result).Norm(two)) / std::abs(expected_result.Norm(two));

    bool fista_test = (relative_error < 100*options.tol);

    std::cout << std::endl << (fista_test ? "Success" : "Failure") << ", achieved ";
    std::cout << relative_error << " relative norm error with a tol of " << options.tol << "." << std::endl;

    return fista_test;
}

bool SmallExample3()
{
    std::cout << "FISTA test with small data : " << std::endl << std::endl;

    astroqut::MatMult<double> A(Matrix(new double[12]{1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0}, 3, 4), 3, 4);
    astroqut::Matrix<double> u(new double[3]{3.0,2.0,1.0}, 3, 1);
    astroqut::Matrix<double> b(new double[3]{1.0,1.0,2.0}, 3, 1);
    astroqut::fista::poisson::Parameters options;
    options.tol = 1e-12;
    options.init_value = astroqut::Matrix(new double[4]{}, 4, 1);
    options.log_period = 50;

    astroqut::Matrix<double> expected_result(new double[4]{0.973633428618360, 0.0, 0.0, -0.674833246032450}, 4, 1);
    astroqut::Matrix<double> actual_result = astroqut::fista::poisson::Solve(A, u, b, 1, options);

    expected_result.Print();
    actual_result.Print();

    double relative_error = std::abs((actual_result - expected_result).Norm(two)) / std::abs(expected_result.Norm(two));

    bool fista_test = (relative_error < 100*options.tol);

    std::cout << std::endl << (fista_test ? "Success" : "Failure") << ", achieved ";
    std::cout << relative_error << " relative norm error with a tol of " << options.tol << "." << std::endl;

    return fista_test;
}

void Time(size_t length)
{
    std::cout << "FISTA test with big data : " << std::endl << std::endl;

    std::default_random_engine generator;
    generator.seed(123456789);
    std::normal_distribution<double> distribution(100.0,10.0);
    size_t test_height = length*length;
    size_t test_width = length;

    double * A_data = new double[test_height*test_width]; // destroyed when A is destroyed
    #pragma omp parallel for
    for( size_t i = 0; i < test_height*test_width; ++i )
    {
        A_data[i] = distribution(generator);
    }
    astroqut::MatMult<double> A(astroqut::Matrix<double>(A_data, test_height, test_width), test_height, test_width);

    double * u_data = new double[test_height]; // destroyed when u is destroyed
    #pragma omp parallel for
    for( size_t i = 0; i < test_height; ++i )
    {
        u_data[i] = distribution(generator);
    }
    astroqut::Matrix<double> u(u_data, test_height, 1);

    double * b_data = new double[test_height]; // destroyed when b is destroyed
    #pragma omp parallel for
    for( size_t i = 0; i < test_height; ++i )
    {
        b_data[i] = distribution(generator);
    }
    astroqut::Matrix<double> b(b_data, test_height, 1);

    astroqut::fista::poisson::Parameters options;
    options.tol = 1e-8;
    options.init_value = astroqut::Matrix(new double[test_width]{}, test_width, 1);
    options.log_period = 1;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();
    astroqut::fista::poisson::Solve(A, u, b, 1, options);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end-start;
    std::cout << std::endl << std::endl << "Time for FISTA solver with " << test_height;
    std::cout << "x" << test_width << " double matrix : " << std::defaultfloat << elapsed_time.count();
    std::cout << " seconds" << std::endl << std::endl;
}

} // namespace matrix
} // namespace test
} // namespace astroqut
