///
/// \file include/test/container.hpp
/// \brief Test suite to validate the DataContainer class.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-30
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_CONTAINER_HPP
#define ASTROQUT_TEST_CONTAINER_HPP

#include <chrono>
#include <iostream>

#include "datacontainer.hpp"

namespace astroqut{
namespace test{
namespace container{

bool TransposeSquare();
bool TransposeRect();

bool Add();
bool Sub();
bool MultSquare();
bool MultRect();

bool NormOne();
bool NormTwo();
bool NormInf();

void Time();

} // namespace container
} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_CONTAINER_HPP
