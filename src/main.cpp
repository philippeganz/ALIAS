///
/// \file src/main.cpp
/// \brief Launcher for the ASTROQUT solver.
/// \details Handle the user input, calls the preparation tools and the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.6.0
/// \date 2019-01-11
/// \copyright GPL-3.0
///

#include "const.hpp"
#include "test.hpp"
#include "WS/astroQUT.hpp"

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>


int main( int argc, char **argv )
{
    astroqut::WS::Parameters<double> options;
    std::string source;
    std::string sensitivity;
    std::string background;
    std::string blurring;
    std::string result;

    int c;

    while (1)
    {
        static struct option long_options[] =
        {
            {"source",      required_argument, nullptr, 'f'},
            {"sensitivity", required_argument, nullptr, 'e'},
            {"background",  required_argument, nullptr, 'o'},
            {"blurring",    required_argument, nullptr, 'b'},
            {"result",      required_argument, nullptr, 'r'},
            {"size",        required_argument, nullptr, 's'},
            {"bootstrap",   required_argument, nullptr, 'x'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "f:e:o:b:r:s:x:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'f':
            source = std::string(optarg);
            break;

        case 'e':
            sensitivity = std::string(optarg);
            break;

        case 'o':
            background = std::string(optarg);
            break;

        case 'b':
            blurring = std::string(optarg);
            break;

        case 'r':
            result = std::string(optarg);
            break;

        case 's':
            options.pic_size = strtoul(optarg, nullptr, 0);
            break;

        case 'x':
            options.bootstrap_max = strtoul(optarg, nullptr, 0);
            break;

        default:
            std::cerr << "usage : ASTROQUT -f source -e sensitivity -o background -b blurring -r result -s size -x bootstrap" << std::endl << std::endl;
            std::cerr << "  source - Path to the source image;" << std::endl;
            std::cerr << "  sensitivity - Path to the sensitivity image;" << std::endl;
            std::cerr << "  background - Path to the background image;" << std::endl;
            std::cerr << "  blurring - Path to the blurring filter;" << std::endl;
            std::cerr << "  result - Path to the solution file;" << std::endl;
            std::cerr << "  size - Width of the picture;" << std::endl;
            std::cerr << "  bootstrap - Amount of bootstraps to perform." << std::endl << std::endl;
            return EXIT_FAILURE;
        }
    }

    astroqut::WS::Solve(source,
                        sensitivity,
                        background,
                        blurring,
                        result,
                        options);

    return EXIT_SUCCESS;
}
