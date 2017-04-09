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
#include <string>

#include "const.hpp"

using namespace std;

int main( int argc, char **argv )
{
    if( argc != 5 )
    {
        cerr << "usage : ASTROQUT <source> <sensitivity> <background> <option file>" << endl << endl;
        cerr << "  source - Path to the source image;" << endl;
        cerr << "  sensitivity - Path to the sensitivity image or an integer >= 1 in which case we consider constant sensitivity;" << endl;
        cerr << "  background - Path to the background image or an integer >= 0 in which case we consider constant background;" << endl;
        cerr << "  option file - Path to the parameters file." << endl << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
