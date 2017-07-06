///
/// \file src/main.cpp
/// \brief Launcher for the ASTROQUT solver.
/// \details Handle the user input, calls the preparation tools and the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-01
/// \copyright GPL-3.0
///

#include <cstdlib>
#include <iostream>
#include <string>

#include "const.hpp"
#include "datacontainer.hpp"
#include "fista.hpp"

int main( int argc, char **argv )
{

//    if( argc != 6 )
//    {
//        std::cerr << "usage : ASTROQUT <source> <sensitivity> <background> <size> <option file>" << std::endl << std::endl;
//        std::cerr << "  source - Path to the source image;" << std::endl;
//        std::cerr << "  sensitivity - Path to the sensitivity image or an integer >= 1 in which case we consider constant sensitivity;" << std::endl;
//        std::cerr << "  background - Path to the background image or an integer >= 0 in which case we consider constant background;" << std::endl;
//        std::cerr << "  option file - Path to the parameters file." << std::endl << std::endl;
//        return EXIT_FAILURE;
//    }

    size_t picture_size = strtol(argv[4], nullptr, 0);
    astroqut::DataContainer<double> source(std::string(argv[1]), picture_size, picture_size);
    astroqut::DataContainer<double> sensitivity(std::string(argv[2]), picture_size, picture_size);
    astroqut::DataContainer<double> background(std::string(argv[3]), picture_size, picture_size);
//    astroqut::Configuration configuration(argv[5]);
//    astroqut::Operators operators(source.width, configuration);

    return EXIT_SUCCESS;
}
