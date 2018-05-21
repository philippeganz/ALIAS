///
/// \file src/test/WS.cpp
/// \brief Test suite to test the WS method.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.3.1
/// \date 2018-05-21
/// \copyright GPL-3.0
///

#include "test/WS.hpp"

namespace astroqut{
namespace test{
namespace astro{

void Chandra()
{
    std::cout << "WS test with 512x512 Chandra data." << std::endl;

    Matrix<double> picture("data/512_chandra/F.data", 512, 512);

    Matrix<double> sensitivity("data/512_chandra/E.data", 512, 512);

    Matrix<double> background("data/512_chandra/O.data", 512, 512);

    Matrix<double> solution = WS::Solve(picture, sensitivity, background);
}

} // namespace astro
} // namespace test
} // namespace astroqut
