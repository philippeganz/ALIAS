///
/// \file src/main.cpp
/// \brief Launcher for the ASTROQUT solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-04-08
/// \copyright GPL-3.0
///

#include <cstdlib>
#include <iostream>

#include "const.hpp"

using namespace std;

int main( int argc, char **argv )
{

    if( argc != 5 )
    {
        cerr << "usage : ASTROQUT <source> <sensitivity> <background> <K matrix> <option file>" << endl << endl;
        cerr << "  source - Image to estimate profile and point sources, must be a square image!" << endl;
        cerr << "  sensitivity - Image of same size as source or a float >= 1 in which case we consider constant sensitivity." << endl;
        cerr << "  background - Image of same size as source or a float >= 0 in which case we consider constant background." << endl;
        cerr << "  K matrix - " << endl;
        cerr << "  option file - Path to the astro parameters option file." << endl << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
