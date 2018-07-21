///
/// \file include/utils/linearop/operator/convolution.hpp
/// \brief Convolution class header
/// \details Provide a convolution operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.5.0
/// \date 2018-04-22
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_CONVOLUTION_HPP
#define ASTROQUT_UTILS_OPERATOR_CONVOLUTION_HPP

#include "utils/linearop/operator.hpp"

namespace astroqut
{

template <class T>
class Convolution : public Operator<T>
{
public:

    /** Default constructor
     */
    Convolution()
        : Operator<T>(0, 0)
    {}

    /** Full member constructor
     *  \param data Matrix containing the filter's data. Needs to be already inverted if not symmetrical.
     */
    Convolution(Matrix<T>&& data)
        : Operator<T>(std::forward<Matrix<T>>(data), data.Height(), data.Width(), false)
    {
        if( this->height_ % 2 != 1 || this->width_ % 2 != 1 )
        {
            throw std::invalid_argument("The filter size for the convolution needs to be odd in both dimensions.");
        }
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    Convolution* Clone() const override final
    {
        return new Convolution(*this);
    }

    /** Default destructor
     */
    virtual ~Convolution()
    {}

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    bool IsValid() const override final
    {
        if( this->height_ != 0 && this->width_ != 0 &&
            !this->data_.IsEmpty() &&
            this->height_ % 2 == 1 && this->width_ % 2 == 1 )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument("Convolution dimensions must be non-zero and odd, and data shall not be empty!");
        }
    }

    Matrix<T> operator*(const Matrix<T>& other) const override final
    {
#ifdef DO_ARGCHECKS
        if( !IsValid() || !other.IsValid() )
        {
            throw std::invalid_argument("Can not perform a convolution with these Matrices.");
        }
#endif // DO_ARGCHECKS

        Matrix<T> result((T) 0, other.Height(), other.Width());

        size_t height_dist_from_center = (this->height_ - 1) / 2;
        size_t width_dist_from_center = (this->width_ - 1) / 2;

#ifdef DEBUG
        int progress_step = std::max(1, (int)other.Height()/100);
#endif // DEBUG

        #pragma omp parallel for
        for( size_t row = 0; row < other.Height(); ++row )
        {
#ifdef DEBUG
            if( row % progress_step == 0 )
                std::cout << "*";
#endif // DEBUG
            int relative_dist_row = row - height_dist_from_center;
            int filter_start_row = relative_dist_row < 0 ? -relative_dist_row : 0;
            int matrix_start_row = relative_dist_row < 0 ? 0 : relative_dist_row;
            for( size_t col = 0; col < other.Width(); ++col )
            {
                int relative_dist_col = col - width_dist_from_center;
                int filter_start_col = relative_dist_col < 0 ? - relative_dist_col : 0;
                int matrix_start_col = relative_dist_col < 0 ? 0 : relative_dist_col;
                for(size_t filter_row = filter_start_row, matrix_row = matrix_start_row;
                    filter_row < this->Height() && matrix_row < other.Height();
                    ++filter_row, ++matrix_row)
                {
                    for(size_t filter_col = filter_start_col, matrix_col = matrix_start_col;
                        filter_col < this->Width() && matrix_col < other.Width();
                        ++filter_col, ++matrix_col)
                    {
                        result[row * other.Width() + col] += other[matrix_row * other.Width() + matrix_col] * this->data_[filter_row * this->width_ + filter_col];
                    }
                }
            }
        }

        return result;
    }
};



} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_CONVOLUTION_HPP
