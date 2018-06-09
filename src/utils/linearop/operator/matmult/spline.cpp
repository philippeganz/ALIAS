///
/// \file include/utils/linearop/operator/spline/generatebasis.cpp
/// \brief Create a basis for the King's spline function
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-06-02
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/matmult/spline.hpp"

namespace astroqut
{

/** Generate a full King spline matrix
 *  \brief Builds an Abel transform matrix with diagonal radius, without duplicating data or inserting zeros
 *  \param pic_side Width of the square picture
 */
Matrix<long double> Spline::Generate( size_t pic_side )
{
    size_t height = 2*pic_side;
    size_t height_sqrt = (size_t) std::sqrt(height);

    Matrix<long double> temp(0.0, height, pic_side);
    Matrix<long double> result(pic_side, pic_side);
    long double* bi = new long double[height_sqrt];
    long double* ri = new long double[height_sqrt];

    for( size_t i = 0; i < height_sqrt; ++i )
    {
        bi[i] = 1.0L/6.0L + i*(10.0L - 1.0L / 6.0L) / ((long double)height_sqrt - 1.0L);
    }

    for( size_t i = 0; i < height_sqrt; ++i )
    {
        ri[i] = 1.0L + i*((long double)pic_side / 2.0L - 1.0L) / ((long double)height_sqrt - 1.0L);
    }

    size_t row = 0;
    for( size_t subrow1 = 0; subrow1 < height_sqrt; ++subrow1 )
    {
        for( size_t subrow2 = 0; subrow2 < height_sqrt; ++subrow2 )
        {
            for( size_t col = 0; col < pic_side/2; ++col )
            {
                double local_result = std::pow(1.0L + std::pow( ((long double)col + 1.0L) / ri[subrow2], 2.0L), -3.0L * bi[subrow1]);
                temp[row*pic_side + pic_side/2 + col] = local_result;
                temp[row*pic_side + pic_side/2 - col - 1] = local_result;
            }
            long double min_value = *std::min_element(&temp[row*pic_side], &temp[(row+1)*pic_side]);
            long double accumulator = 0.0L;
            for( size_t i = 0; i < pic_side; ++i )
            {
                temp[row*pic_side + i] -= min_value;
                accumulator += temp[row*pic_side + i];
            }
            for( size_t i = 0; i < pic_side; ++i )
            {
                temp[row*pic_side + i] /= accumulator;
            }
            ++row;
        }
    }

    for( size_t i = 0; i < pic_side/2; ++i )
    {
        for( size_t j = pic_side/2; j < pic_side; ++j )
        {
            temp[i*pic_side + j] = 0.0L;
        }
    }

    for( size_t i = pic_side/2; i < height; ++i )
    {
        for( size_t j = 0; j < pic_side/2; ++j )
        {
            temp[i*pic_side + j] = 0.0L;
        }
    }

    for( size_t i = 0; i < pic_side; ++i )
    {
        size_t index = std::floor(i*(height_sqrt*height_sqrt-1.0)/(pic_side-1.0));
        std::copy( &temp[index*pic_side], &temp[(index+1)*pic_side], &result[i*pic_side]);
    }

    delete[] bi;
    delete[] ri;
    return result.Transpose();
}

} // namespace astroqut
