///
/// \file include/utils/linearop/operator/matmult/spline.hpp
/// \brief Spline operator class header
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-04-22
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_SPLINE_HPP
#define ASTROQUT_UTILS_OPERATOR_SPLINE_HPP

#include "utils/linearop/operator/matmult.hpp"

namespace astroqut
{

namespace spline
{

Matrix<double> Generate(size_t pic_side);

} // namespace spline

class Spline : public MatMult<double>
{
private:

public:

    /** Default constructor
     */
    Spline()
    {}

    /** Full member constructor
     *  \param low_pass_filter QMF matrix for low pass filtering
     *  \param low_pass_filter Mirrored QMF matrix for high pass filtering
     *  \param height Height of the full Abel matrix
     *  \param width Width of the full Abel matrix
     */
    Spline(Matrix<double>&& data, size_t height, size_t width)
        : MatMult<double>(std::forward<Matrix<double>>(data), height, width)
    {}

    /** Build constructor
     *  \brief Builds the Spline operator with the qmf matrix corresponding to type and parameter
     *  \param type Spline type, can be one of haar, beylkin, coiflet, daubechies, symmlet, vaidyanathan, battle
     *  \param parameter Integer parameter specific to each wavelet type
     */
    Spline(size_t pic_size)
        : MatMult<double>(spline::Generate(pic_size), pic_size, pic_size)
    {}

    /** Default destructor
     */
    virtual ~Spline()
    {}
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_SPLINE_HPP
