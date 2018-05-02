///
/// \file src/test/WS.cpp
/// \brief Test suite to test the WS method.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.3.0
/// \date 2018-05-01
/// \copyright GPL-3.0
///

#include "test/WS.hpp"

namespace astroqut{
namespace test{
namespace astro{

void Chandra()
{
    std::ifstream file("data/512_chandra/F.data", std::ios::binary | std::ios::in | std::ios::ate);
    Matrix<double> picture(512, 512);
    file >> picture;
    file.close();

    std::ifstream file2("data/512_chandra/E.data", std::ios::binary | std::ios::in | std::ios::ate);
    Matrix<double> sensitivity(512, 512);
    file2 >> sensitivity;
    file2.close();

    std::ifstream file3("data/512_chandra/O.data", std::ios::binary | std::ios::in | std::ios::ate);
    Matrix<double> background(512, 512);
    file3 >> background;
    file3.close();

    Matrix<double> solution = WS::Solve(picture, sensitivity, background);
}

} // namespace astro
} // namespace test
} // namespace astroqut
