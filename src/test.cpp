///
/// \file src/test.cpp
/// \brief Implementation of the test suites.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.2.0
/// \date 2018-01-03
/// \copyright GPL-3.0
///

#include "test.hpp"

namespace astroqut{
namespace test{

bool Matrix()
{
    matrix::Time(2048, matrix::double_floating);

    matrix::Optimizations(1000000000);

    bool transpose_square = matrix::TransposeSquare();
    bool transpose_rect = matrix::TransposeRect();

    bool add = matrix::Add();
    bool sub = matrix::Sub();
    bool mult_square = matrix::MultSquare();
    bool mult_rect = matrix::MultRect();
    bool vect_mat = matrix::MultVectMat();
    bool mat_vect = matrix::MultMatVect();

    bool norm_one = matrix::NormOne();
    bool norm_two = matrix::NormTwo();
    bool norm_inf = matrix::NormInf();

    bool sum = matrix::Sum();

    bool shrink = matrix::Shrink();

    return transpose_square && transpose_rect && add && sub && mult_square && mult_rect && vect_mat && mat_vect && norm_one && norm_two && norm_inf && sum && shrink;
}

bool FISTA()
{
    bool fista_small1 = fista::SmallExample1();
    bool fista_small2 = fista::SmallExample2();
    bool fista_small3 = fista::SmallExample3();

    fista::Time(1024);

    return fista_small1 && fista_small2 && fista_small3;
}

} // namespace test
} // namespace astroqut
