///
/// \file include/test.hpp
/// \brief Test suite to validate the DataContainer class.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-29
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_DATACONTAINER_HPP
#define ASTROQUT_TEST_DATACONTAINER_HPP

#include <chrono>
#include <iostream>

#include <datacontainer.hpp>

namespace astroqut{
namespace test{
namespace datacontainer{

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

} // namespace datacontainer
} // namespace test
} // namespace astroqut

#endif // ASTROQUT_TEST_DATACONTAINER_HPP
