///
/// \file include/test/fista.hpp
/// \brief Test suite to validate the FISTA class.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.2.0
/// \date 2018-01-04
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_FISTA_HPP
#define ASTROQUT_TEST_FISTA_HPP

#include "fista/poisson.hpp"
#include "settings.hpp"

#include <chrono>
#include <iostream>
#include <random>

namespace astroqut{
namespace test{
namespace fista{

bool SmallExample1();
bool SmallExample2();
bool SmallExample3();

void Time(size_t length);

} // namespace fista
} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_FISTA_HPP
