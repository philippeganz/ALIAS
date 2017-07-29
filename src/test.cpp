///
/// \file src/test.cpp
/// \brief Implementation of the test suites.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-20
/// \copyright GPL-3.0
///

#include "test.hpp"

namespace astroqut{
namespace test{

bool dataContainer()
{
    bool transpose_square = datacontainer::TransposeSquare();
    bool transpose_rect = datacontainer::TransposeRect();

    bool add = datacontainer::Add();
    bool sub = datacontainer::Sub();
    bool mult_square = datacontainer::MultSquare();
    bool mult_rect = datacontainer::MultRect();

    bool norm_one = datacontainer::NormOne();
    bool norm_two = datacontainer::NormTwo();
    bool norm_inf = datacontainer::NormInf();

    datacontainer::Time();

    return transpose_square && transpose_rect && add && sub && mult_square && mult_rect && norm_one && norm_two && norm_inf;
}

} // namespace test
} // namespace astroqut
