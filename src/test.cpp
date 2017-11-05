///
/// \file src/test.cpp
/// \brief Implementation of the test suites.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-30
/// \copyright GPL-3.0
///

#include "test.hpp"

namespace astroqut{
namespace test{

bool DataContainer()
{
    bool transpose_square = container::TransposeSquare();
    bool transpose_rect = container::TransposeRect();

    bool add = container::Add();
    bool sub = container::Sub();
    bool mult_square = container::MultSquare();
    bool mult_rect = container::MultRect();

    bool norm_one = container::NormOne();
    bool norm_two = container::NormTwo();
    bool norm_inf = container::NormInf();

//    container::Time();

    return transpose_square && transpose_rect && add && sub && mult_square && mult_rect && norm_one && norm_two && norm_inf;
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
