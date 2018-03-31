///
/// \file include/utils/linearop/operator/abeltransform.hpp
/// \brief Abel transform class header
/// \details Provide the Abel transform operator
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-03-30
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_ABELTRANSFORM_HPP
#define ASTROQUT_UTILS_OPERATOR_ABELTRANSFORM_HPP

#include "utils/linearop/operator.hpp"

#ifdef DEBUG
    #include <sstream>
#endif // DEBUG

namespace astroqut
{

template <class T>
class AbelTransform : public Operator<T>
{
private:
    size_t pic_side_;
    size_t wavelet_amount_;

public:

    /** Default constructor
     */
    AbelTransform()
        : Operator<T>(0, 0)
    {}

    /** Full member constructor
     *  \param data Matrix containing the upper left component of the Abel transform.
     *  \param height Height of the full Abel matrix
     *  \param width Width of the full Abel matrix
     */
    AbelTransform(Matrix<T>&& data, size_t height, size_t width)
        : Operator<T>(std::forward<Matrix<T>>(data), height, width, false)
    {}

    /** Build constructor
     *  \brief Builds an Abel transform matrix with diagonal radius
     *  \param wavelet_amount Power of 2 amount of wavelets, typically (pixel_amount/2)^2
     *  \param pixel_amount Total amount of pixels of the target picture
     *  \param radius Amount of pixels from centre to border of galaxy, typically pixel_amount/2
     */
    AbelTransform(unsigned int wavelets_amount, unsigned int pixel_amount, unsigned int radius)
        : Operator<T>(Matrix<T>((T) 0, pixel_amount/4, wavelets_amount/2), pixel_amount, wavelets_amount, false)
        , pic_side_(std::sqrt(pixel_amount))
        , wavelet_amount_(wavelets_amount)
    {
        size_t pic_side_half = pic_side_/2;
        size_t pic_side_extended = std::floor(pic_side_half*std::sqrt(2.0));
        size_t wavelet_amount_half = wavelet_amount_/2;
        double radius_extended = radius * std::sqrt(2.0);
        double radius_extended_to_pic_side_extended_ratio = radius_extended/pic_side_extended;
        double radius_to_pic_side_ratio = radius/pic_side_half;
        double radius_extended_to_wavelet_amount_half_ratio = radius_extended/wavelet_amount_half;
        double* x_axis = new double[wavelet_amount_half];
        for( size_t i = 0; i < wavelet_amount_half; ++i )
        {
            x_axis[i] = (i+1) * radius_extended_to_wavelet_amount_half_ratio;
        }

        #pragma omp parallel for
        for( size_t i = 0; i < pic_side_half; ++i )
        {
            double z = i * radius_to_pic_side_ratio;
            for( size_t j = 0; j < pic_side_half; ++j )
            {
                double y = j * radius_extended_to_pic_side_extended_ratio;
                double s = std::sqrt(y*y + z*z);
                for(unsigned int k = 0; k < wavelet_amount_half; ++k)
                {
                    if( x_axis[k] <= s )
                    {
                        continue;
                    }

                    double ri0 = s;
                    if( k != 0 && x_axis[k-1] >= s )
                    {
                        ri0 = x_axis[k-1];
                    }

                    double ri1 = x_axis[k];
                    if( x_axis[k] < s )
                    {
                        ri1 = s;
                    }

                    size_t index = (wavelet_amount_half-i-1)*pic_side_half*wavelet_amount_half + (pic_side_half-j-1)*wavelet_amount_half + wavelet_amount_half-k-1;
                    this->data_[index] = 2*(std::sqrt(ri1*ri1 - s*s) - std::sqrt(ri0*ri0 - s*s));
                }
            }
        }
        delete[] x_axis;
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    AbelTransform* Clone() const override final
    {
        return new AbelTransform(*this);
    }

    /** Default destructor
     */
    virtual ~AbelTransform()
    {}

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    bool IsValid() const override final
    {
        if( this->height_ != 0 && this->width_ != 0 &&
            !this->data_.IsEmpty() )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument("Operator dimensions must be non-zero and function shall not be nullptr!");
        }
    }

    Matrix<T> operator*(const Matrix<T>& other) const override final
    {
#ifdef DO_ARGCHECKS
        try
        {
            this->ArgTest(other, mult);
        }
        catch (const std::exception&)
        {
            throw;
        }
#endif // DO_ARGCHECKS

        size_t pic_side_half = pic_side_/2;
        size_t wavelet_amount_half = wavelet_amount_/2;
        Matrix<T> result( (T) 0, this->Height(), other.Width() );

#ifdef DEBUG
        int progress_step = std::max(1, (int)(pic_side_half*pic_side_half)/100);
        int step = 0;
        std::cout << std::endl;
#endif // DEBUG

        // iterating over blocks
        #pragma omp parallel for
        for( size_t block = 0; block < pic_side_half; ++block )
        {
            // iterating over rows
            for( size_t i = 0; i < pic_side_half; ++i )
            {
#ifdef DEBUG
                if( (block*pic_side_half+i) % progress_step == 0 )
                {
                    std::stringstream output;
                    output << "\r" << step++;
                    std::cout << output.str();
                }
#endif // DEBUG

                // iterating over matrix multiplication vectors
                for( size_t k = 0; k < wavelet_amount_half; ++k )
                {
                    size_t abel_left_index = (block*pic_side_half + i)*wavelet_amount_half + k;
                    size_t abel_right_index = (pic_side_half*(pic_side_half - block - 1) + i)*wavelet_amount_half + wavelet_amount_half - k - 1;

                    // iterating over columns
                    for( size_t j = 0; j < other.Width(); ++j )
                    {

                        // left part
                        size_t target_left_index = k*other.Width() + j;

                        T left_temp = this->data_[abel_left_index] * other[target_left_index];

                        size_t result_left_upper_index = (block*pic_side_ + i)*other.Width() + j;
                        result[result_left_upper_index] += left_temp;

                        size_t result_left_lower_index = (block*pic_side_ + pic_side_ - i - 1)*other.Width() + j;
                        result[result_left_lower_index] += left_temp;


                        // right part
                        size_t target_right_index = (wavelet_amount_half + k)*other.Width() + j;

                        T right_temp = this->data_[abel_right_index] * other[target_right_index];

                        size_t result_right_upper_index = ((block + pic_side_half)*pic_side_ + i)*other.Width() + j;
                        result[result_right_upper_index] += right_temp;

                        size_t result_right_lower_index = ((block + pic_side_half)*pic_side_ + pic_side_ - i - 1)*other.Width() + j;
                        result[result_right_lower_index] += right_temp;
                    }
                }
            }
        }

        return result;
    }
};



} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_ABELTRANSFORM_HPP
