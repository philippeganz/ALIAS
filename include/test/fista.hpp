///
/// \file include/test/fista.hpp
/// \brief Test suite to validate the FISTA class.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-30
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_FISTA_HPP
#define ASTROQUT_TEST_FISTA_HPP

#include <chrono>
#include <iostream>
#include <random>

#include "fista/poisson.hpp"

namespace astroqut{
namespace test{
namespace fista{

bool SmallExample1();
bool SmallExample2();
bool SmallExample3();

void Time();

} // namespace datacontainer
} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_FISTA_HPP
