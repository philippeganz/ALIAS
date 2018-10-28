///
/// \file src/main.cpp
/// \brief Launcher for the ASTROQUT solver.
/// \details Handle the user input, calls the preparation tools and the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.6.0
/// \date 2018-10-27
/// \copyright GPL-3.0
///

#include "const.hpp"
#include "test.hpp"
#include "WS/astroQUT.hpp"

#include <cstdlib>
#include <iostream>
#include <string>


int main( int argc, char **argv )
{
//
//    try
//    {
//        astroqut::test::MatrixTest<double>();
//    }
//    catch (const std::exception& err)
//    {
//        std::cerr << err.what() << std::endl;
//        std::cerr << "Matrix tests failed! Please refer to the individual test results for more details." << std::endl;
//        return EXIT_FAILURE;
//    }
//
//    try
//    {
//        astroqut::test::OperatorTest();
//    }
//    catch (const std::exception& err)
//    {
//        std::cerr << err.what() << std::endl;
//        std::cerr << "Operator tests failed! Please refer to the individual test results for more details." << std::endl;
//        return EXIT_FAILURE;
//    }
//
//    try
//    {
//        astroqut::test::FISTATest();
//    }
//    catch (const std::exception& err)
//    {
//        std::cerr << err.what() << std::endl;
//        std::cerr << "FISTA tests failed! Please refer to the individual test results for more details." << std::endl;
//        return EXIT_FAILURE;
//    }
//
//    astroqut::test::astro::Chandra();
//    astroqut::test::astro::Sym256();

//    if( argc != 6 )
//    {
//        std::cerr << "usage : ASTROQUT <source> <sensitivity> <background> <size> <option file>" << std::endl << std::endl;
//        std::cerr << "  source - Path to the source image;" << std::endl;
//        std::cerr << "  sensitivity - Path to the sensitivity image or an integer >= 1 in which case we consider constant sensitivity;" << std::endl;
//        std::cerr << "  background - Path to the background image or an integer >= 0 in which case we consider constant background;" << std::endl;
//        std::cerr << "  solution - Path to the solution file;" << std::endl;
//        std::cerr << "  size - Width of the picture;" << std::endl;
//        std::cerr << "  option file - Path to the parameters file." << std::endl << std::endl;
//        return EXIT_FAILURE;
//    }

    astroqut::WS::Parameters<double> options;
    options.pic_size = strtoul(argv[5], nullptr, 0);
    options.MC_max = strtoul(argv[6], nullptr, 0);
    options.bootstrap_max = strtoul(argv[7], nullptr, 0);

    astroqut::WS::Solve(std::string(argv[1]),
                        std::string(argv[2]),
                        std::string(argv[3]),
                        std::string(argv[4]),
                        options);

    return EXIT_SUCCESS;
}
