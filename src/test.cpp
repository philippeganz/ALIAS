///
/// \file src/test.cpp
/// \brief Implementation of the test suites.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.2.0
/// \date 2017-12-28
/// \copyright GPL-3.0
///

#include "test.hpp"

namespace astroqut{
namespace test{

bool Matrix()
{
    bool transpose_square = matrix::TransposeSquare();
    bool transpose_rect = matrix::TransposeRect();

    bool add = matrix::Add();
    bool sub = matrix::Sub();
    bool mult_square = matrix::MultSquare();
    bool mult_rect = matrix::MultRect();

    bool norm_one = matrix::NormOne();
    bool norm_two = matrix::NormTwo();
    bool norm_inf = matrix::NormInf();

    bool shrink = matrix::Shrink();

//    matrix::Time();

    return transpose_square && transpose_rect && add && sub && mult_square && mult_rect && norm_one && norm_two && norm_inf && shrink;
}

bool FISTA()
{
    bool fista_small1 = fista::SmallExample1();
    bool fista_small2 = fista::SmallExample2();
    bool fista_small3 = fista::SmallExample3();

//    fista::Time();

    return fista_small1 && fista_small2 && fista_small3;
}

} // namespace test
} // namespace astroqut
