///
/// \file include/test/matrix.hpp
/// \brief Test suite to validate the Matrix class.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.2.0
/// \date 2018-01-04
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_MATRIX_HPP
#define ASTROQUT_TEST_MATRIX_HPP

#include "utils/matrix.hpp"
#include "settings.hpp"

#include <chrono>
#include <iostream>
#include <omp.h>

namespace astroqut{
namespace test{
namespace matrix{

bool TransposeSquare();
bool TransposeRect();

bool Add();
bool Sub();
bool MultSquare();
bool MultRect();
bool MultVectMat();
bool MultMatVect();

bool NormOne();
bool NormTwo();
bool NormInf();

bool Sum();

bool Shrink();

enum TimeTestType{integer, long_integer, floating, double_floating};
void Time(size_t length, TimeTestType type);

void Optimizations(size_t length);

} // namespace matrix
} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_MATRIX_HPP
