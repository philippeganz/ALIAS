///
/// \file include/test/matrix.hpp
/// \brief Test suite to validate the Matrix class.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.2.0
/// \date 2017-12-28
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_MATRIX_HPP
#define ASTROQUT_TEST_MATRIX_HPP

#include "utils/functions.hpp"
#include "utils/matrix.hpp"

#include <chrono>
#include <iostream>

namespace astroqut{
namespace test{
namespace matrix{

bool TransposeSquare();
bool TransposeRect();

bool Add();
bool Sub();
bool MultSquare();
bool MultRect();

bool NormOne();
bool NormTwo();
bool NormInf();

bool Shrink();

void Time();

} // namespace matrix
} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_MATRIX_HPP
