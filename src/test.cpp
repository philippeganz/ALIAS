///
/// \file src/test.cpp
/// \brief Implementation of the test suites.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-01-12
/// \copyright GPL-3.0
///

#include "test.hpp"

namespace astroqut{
namespace test{

bool FISTA()
{
    bool fista_small = fista::SmallExample();

    fista::Time(1024);

    return fista_small;
}

} // namespace test
} // namespace astroqut
