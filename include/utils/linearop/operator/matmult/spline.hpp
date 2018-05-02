///
/// \file include/utils/linearop/operator/matmult/spline.hpp
/// \brief Spline operator class header
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-05-01
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
     *  \param pic_size Side size of the picture in pixel
     */
    Spline(size_t pic_size)
        : MatMult<double>(spline::Generate(pic_size), pic_size, pic_size)
    {}

    /** Default destructor
     */
    virtual ~Spline()
    {}

    /** Transpose in-place
     *   \return A reference to this
     */
    Spline& Transpose() override final
    {
        std::swap(this->height_, this->width_);
        this->transposed_ = !this->transposed_;
        this->data_ = std::move(this->data_.Transpose());
        return *this;
    }
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_SPLINE_HPP
