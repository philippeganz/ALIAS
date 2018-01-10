///
/// \file include/test.hpp
/// \brief Test suites to validate the project code
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.2.0
/// \date 2017-12-28
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_HPP
#define ASTROQUT_TEST_HPP

#include "test/matrix.hpp"
#include "test/fista.hpp"

namespace astroqut{
namespace test{

bool Matrix();
bool FISTA();

} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_HPP
