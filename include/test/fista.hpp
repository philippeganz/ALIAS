///
/// \file include/test/fista.hpp
/// \brief Test suite to validate the FISTA class.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.3.0
/// \date 2018-01-12
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

bool SmallExample();

void Time(size_t length);

} // namespace fista
} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_FISTA_HPP
