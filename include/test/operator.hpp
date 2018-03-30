///
/// \file include/test/operator.hpp
/// \brief Test suite to validate the Operator class.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.3.0
/// \date 2018-02-25
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_OPERATOR_HPP
#define ASTROQUT_TEST_OPERATOR_HPP

#include "utils/linearop/operator/convolution.hpp"
#include "utils/linearop/operator/abeltransform.hpp"

#include <chrono>
#include <iostream>
#include <random>

namespace astroqut{
namespace test{
namespace oper{

bool ConvolutionTest();
void ConvolutionTime(size_t data_length, size_t filter_length);

bool AbelTestBuild();
bool AbelTestApply();
bool AbelTime();

} // namespace oper
} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_OPERATOR_HPP
