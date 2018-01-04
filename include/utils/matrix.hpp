///
/// \file include/utils/matrix.hpp
/// \brief Matrix class header
/// \details Provide matrix container with multiple matrix operations used in the whole project.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.2.0
/// \date 2018-01-04
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_MATRIX_HPP
#define ASTROQUT_UTILS_MATRIX_HPP

#include "utils/linearop.hpp"
#include "utils/matrix/functions.hpp"

//#include <fstream>
#include <iomanip>
#include <iostream>
//#include <sstream>
#include <string>

namespace astroqut
{

/** Types of norm currently implemented
 */
enum NormType {one, two, two_squared, inf};

template <class T>
class Matrix : public LinearOp<T>
{
private:
    T* data_; //!< Member variable "data_"

public:

    /** Default constructor
     *  Create an empty container
     */
    Matrix() noexcept
        : data_(nullptr)
    {}

    /** Empty constructor
     *  Create a container with default-initialized data
     *  \param height Height of the data
     *  \param width Width of the data
     */
    Matrix( const size_t height, const size_t width)
        : LinearOp<T>(height, width)
        , data_(nullptr)
    {
        try
        {
            data_ = new T[this->length_];
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for new array!" << std::endl;
            throw;
        }
    }

    /** Full member constructor
     *  \param data 2D array containing the pixels
     *  \param height Height of the data
     *  \param width Width of the data
     */
    Matrix( T* const data, const size_t height, const size_t width) noexcept
        : LinearOp<T>(height, width)
        , data_(data)
    {}

    /** Constant number constructor
     *  \param number Number to fill data with
     *  \param height Height of the data
     *  \param width Width of the data
     */
    Matrix( const T number, const size_t height, const size_t width)
        : Matrix(height, width)
    {
        // TODO std::fill not yet parallelized : change as soon as available
        // std::fill( data_, data_ + (this->length_), number );
        #pragma omp parallel for
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] = number;
        }
    }

// TODO
//    /** File constructor
//     *  \param file_path Path to the data file
//     *  \param height Height of the data
//     *  \param width Width of the data
//     */
//
//    Matrix(std::string& file_path, const size_t height, const size_t width)
//        : Matrix(height, width)
//    {
//        std::ifstream file;
//        file.exceptions( std::ifstream::failbit | std::ifstream::badbit );
//        try
//        {
//            file.open(file_path, std::ifstream::in);
//        }
//        catch (const std::ifstream::failure& e)
//        {
//            std::cerr << "Could not open " << file_path << "! Please verify that the file exists." << std::endl;
//            throw;
//        }
//        std::string line;
//        T number;
//        for( size_t i = 0; getline(file, line); ++i )
//        {
//            std::stringstream line_stream(line);
//            for( size_t j = 0; line_stream >> number; ++j )
//            {
//                data_[i*width_ + j] = number;
//            }
//        }
//        file.close();
//    }


    /** Copy constructor
     *  \param other Object to copy from
     */
    Matrix(const Matrix<T>& other)
        : Matrix(other.height_, other.width_)
    {
        *this = other;
    }

    /** Move constructor
     *  \param other Object to move from
     */
    Matrix(Matrix<T>&& other) noexcept
        : Matrix(other.data_, other.height_, other.width_)
    {
        other.Data(nullptr);
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    Matrix* Clone() const
    {
        return new Matrix(*this);
    }

    /** Default destructor */
    virtual ~Matrix()
    {
        if( data_ != nullptr )
        {
            delete[] data_;
        }
        data_ = nullptr;
    }

    /** Access data_
     * \return The current value of data_
     */
    T* Data() const noexcept
    {
        return data_;
    }
    /** Set data_
     * \param data New value to set
     */
    void Data(T* const data) noexcept
    {
        data_ = data;
    }

    /** Empty test operator
     *  \return True if empty.
     */
    bool IsEmpty() const noexcept
    {
        if( this->height_ == 0 &&
            this->width_ == 0 &&
            data_ == nullptr)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    bool IsValid() const override final
    {
        if( this->height_ != 0 &&
            this->width_ != 0 &&
            data_ != nullptr )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must be non-zero and data shall not be empty!");
        }
    }

    /** Negativity test
     *  \return True is a negative number is found in the array, False otherwise.
     */
    bool ContainsNeg() const noexcept
    {
        if( this->height_ == 0 ||
            this->width_ == 0 ||
            data_ == nullptr)
        {
            return false;
        }

        for( size_t i = 0; i < this->length_; ++i )
        {
            if( data_[i] < 0 )
            {
                return true;
            }
        }

        return false;
    }

    /** Comparison operator equal
     *  \param other Object to compare to
     *  \return True if both object are the same element-wise, False else
     */
    bool operator==(const Matrix& other) const noexcept
    {
        if( this->height_ != other.height_ ||
            this->width_ != other.width_ )
        {
            return false;
        }

        if( data_ == other.data_ )
        {
            return true;
        }
        else
        {
            // quick check for first value
            if( !IsEqual(data_[0], other.data_[0]) )
            {
                return false;
            }

            // if first values are the same, check for the rest in parallel
            bool are_they_equal = true;
            #pragma omp parallel
            {
                #pragma omp for
                for(size_t i = 0; i < this->length_; ++i)
                {
                    if( !IsEqual(data_[i], other.data_[i]) )
                    {
                        #pragma omp critical
                        {
                            are_they_equal = false;
                        }
                        #pragma omp cancel for
                    }
                }
            }
            return are_they_equal;

            // TODO std::equal not yet parallelized : change as soon as available
            // std::equal(data_, data_+(this->length_), other.data_);
        }
    }

    /** Comparison operator not-equal
     *  \param other Object to compare to
     *  \return False if both object are the same element-wise, True else
     */
    bool operator!=(const Matrix& other) const noexcept
    {
        return !(*this == other);
    }

    /** Assignment operator
     *  \param other Object to assign to current object
     *  \return A reference to this
     */
    Matrix& operator=(const Matrix& other)
    {
        if( this != &other )
        {
            // we need to allocate memory if data_ has not yet been allocated
            // or reallocate if the array size is not the same
            if( data_ == nullptr || this->length_ != other.length_ )
            {
                if( data_ != nullptr )
                {
                    delete[] data_;
                }

                try
                {
                    data_ = new T[other.length_];
                }
                catch (const std::bad_alloc& ba)
                {
                    std::cerr << "Could not allocate memory for new array!" << std::endl;
                    throw;
                }
            }

            this->height_ = other.height_;
            this->width_ = other.width_;
            this->length_ = other.length_;

            #pragma omp parallel for simd
            for(size_t i = 0; i < other.length_; ++i)
            {
                data_[i] = other.data_[i];
            }

            // TODO std::copy not yet parallelized : change as soon as available
            // std::copy(other.data_, other.data_ + other.length_, data_);
        }
        return *this;
    }

    /** Move operator
     *  \param other Object to move to current object
     *  \return A reference to this
     */
    Matrix& operator=(Matrix&& other) noexcept
    {
        if( this != &other )
        {
            this->height_ = other.height_;
            this->width_ = other.width_;
            this->length_ = other.length_;
            std::swap(data_, other.data_);
        }
        return *this;
    }

    /** Additive operator
     *  \param other Object to add to current object
     *  \return A new instance containing the result
     */
    Matrix operator+(const Matrix& other) const &
    {
        try
        {
            this->ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(other.height_, other.width_);

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            result.data_[i] = data_[i] + other.data_[i];
        }

        return result;
    }

    /** Additive operator of temporary instance
     *  \param other Object to add to current object
     *  \return Current object overwritten with the result
     */
    Matrix operator+(const Matrix& other) &&
    {
        *this += other; // if object is temporary, no need to allocate new memory for the result
        return *this;
    }

    /** Additive operator in-place
     *  \param other Object to add from current object
     *  \return A reference to this
     */
    Matrix& operator+=(const Matrix& other)
    {
        try
        {
            this->ArgTest(other, add);
        }
        catch (const std::exception&)
        {
            throw;
        }

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] += other.data_[i];
        }

        return *this;
    }

    /** Subtractive operator
     *  \param other Object to subtract from current object
     *  \return A new instance containing the result
     */
    Matrix operator-(const Matrix& other) const &
    {
        try
        {
            this->ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(other.height_, other.width_);

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            result.data_[i] = data_[i] - other.data_[i];
        }

        return result;
    }

    /** Subtractive operator of temporary instance
     *  \param other Object to subtract from current object
     *  \return Current object overwritten with the result
     */
    Matrix operator-(const Matrix& other) &&
    {
        *this -= other; // if object is temporary, no need to allocate new memory for the result
        return *this;
    }

    /** Subtractive operator in-place
     *  \param other Object to subtract current object to
     *  \return A reference to this
     */
    Matrix& operator-=(const Matrix& other)
    {
        try
        {
            this->ArgTest(other, add);
        }
        catch (const std::exception&)
        {
            throw;
        }

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] -= other.data_[i];
        }

        return *this;
    }

    /** Multiplicative operator
     *  \param other Object to multiply current object to
     *  \return A new instance containing the result
     */
    Matrix operator*(const Matrix& other) const
    {
        try
        {
            this->ArgTest(other, mult);
        }
        catch (const std::exception&)
        {
            throw;
        }

        // other is a vector
        if( other.height_ == 1 || other.width_ == 1 )
        {
            // this is a vector
            if( this->height_ == 1 || this->width_ == 1 )
            {
                return Matrix(Inner(other), 1, 1);
            }
            // this is a matrix
            else
            {
                Matrix result((T) 0, this->height_, 1); // init to zero
                MatrixVectorMult(*this, other, result);
                return result;
            }
        }
        // other is a matrix
        else
        {
            // this is a vector
            if( this->height_ == 1 || this->width_ == 1 )
            {
                Matrix result((T) 0, 1, other.width_); // init to zero
                VectorMatrixMult(*this, other, result);
                return result;
            }
            // this is a matrix
            else
            {
                Matrix result((T) 0, this->height_, other.width_); // init to zero
                MatrixMatrixMult(*this, other, result);
                return result;
            }
        }
    }

    /** Multiplicative operator with single number
     *  \param number Number to multiply the current object with
     *  \return A new instance containing the result
     */
    Matrix operator*(const T number) const &
    {
        try
        {
            IsValid();
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(this->height_, this->width_);

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            result.data_[i] = data_[i] * number;
        }

        return result;
    }

    /** Multiplicative operator of temporary instance with single number
     *  \param number Number to multiply the current object with
     *  \return A reference to this
     */
    Matrix operator*(const T number) &&
    {
        try
        {
            IsValid();
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] *= number;
        }

        return *this;
    }

    /** Multiplicative operator in-place
     *  \param other Object to multiply to current object
     *  \return A reference to this
     */
    Matrix& operator*=(const Matrix& other)
    {
        *this = *this * other;
        return *this;
    }

    /** Element-wise multiply
     *  \param other Object to multiply to current object, element-wise
     *  \return A new instance containing the result
     */
    Matrix operator->*(const Matrix& other) const &
    {
        try
        {
            this->ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(this->height_, this->width_);

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            result.data_[i] = data_[i] * other.data_[i];
        }

        return result;
    }

    /** Element-wise multiply of temporary instance
     *  \param other Object to multiply to current object, element-wise
     *  \return A reference to this
     */
    Matrix operator->*(const Matrix& other) &&
    {
        try
        {
            this->ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] *= other.data_[i];
        }

        return *this;
    }

    /** Element-wise divide operator
     *  \param other Object to divide from current object, element-wise
     *  \return A new instance containing the result
     */
    Matrix operator/(const Matrix& other) const &
    {
        try
        {
            this->ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(this->height_, this->width_);

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            result.data_[i] = (double) data_[i] / (double) other.data_[i];
        }

        return result;
    }

    /** Element-wise divide operator of temporary instance
     *  \param other Object to divide from current object, element-wise
     *  \return A reference to this
     */
    Matrix operator/(const Matrix& other) &&
    {
        try
        {
            this->ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] = (double) data_[i] / (double) other.data_[i];
        }

        return *this;
    }

    /** Divide operator with single number
     *  \param number Number to divide the current object by
     *  \return A new instance containing the result
     */
    Matrix operator/(const T number) const &
    {
        try
        {
            IsValid();
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(this->height_, this->width_);

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            result.data_[i] = (double)data_[i] / (double)number;
        }

        return result;
    }

    /** Divide operator with single number of temporary instance
     *  \param number Number to divide the current object by
     *  \return A reference to this
     */
    Matrix operator/(const T number) &&
    {
        try
        {
            IsValid();
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] = (double)data_[i] / (double)number;
        }

        return *this;
    }

    /** Transpose
    *   \return A new instance with the transposed matrix
    */
    Matrix Transpose() const &
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        Matrix result(this->width_, this->height_); // width <--> height

        // vector
        if( this->height_ == 1 || this->width_ == 1 )
        {
            std::copy(data_, data_ + (this->length_), result.data_);
        }
        // matrix
        else
        {
            #pragma omp parallel for
            for ( size_t i = 0; i < this->height_; ++i )
			{
			    #pragma omp simd
				for ( size_t j = 0 ; j < this->width_; ++j )
				{
					result.data_[ j * this->height_ + i ] = data_[ i * this->width_ + j ];
				}
			}
        }
        return result;
    }

    /** Transpose in-place
    *   \return A reference to this
    */
    Matrix& Transpose() &&
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        // vector
        if( this->height_ == 1 || this->width_ == 1 )
        {
            std::swap(this->height_, this->width_);
        }
        // matrix
        else
        {
            // square matrix, simple swap
            if( this->height_ == this->width_ )
            {

                #pragma omp parallel for
                for ( size_t i = 0; i < this->height_; ++i )
                {
                    #pragma omp simd
                    for ( size_t j = i+1 ; j < this->width_; ++j )
                    {
                        std::swap(data_[i * this->width_ + j], data_[j * this->width_ + i]);
                    }
                }
                std::swap(this->height_, this->width_); // width <--> height
            }

            // rect matrix, cyclic permutations
            // TODO in-place algorithm
            // for now, simple copy to new array
            else
            {
                Matrix result(this->width_, this->height_); // width <--> height

                #pragma omp parallel for
                for ( size_t i = 0; i < this->height_; ++i )
                {
                    #pragma omp simd
                    for ( size_t j = 0 ; j < this->width_; ++j )
                    {
                        result.data_[ j * this->height_ + i ] = data_[ i * this->width_ + j ];
                    }
                }
                *this = std::move(result);
            }
        }
        return *this;
    }

    /** Log operator
     *  Applies the log function to all positive elements, replace by zero otherwise
     *  \return A new instance containing the result
     */
    Matrix Log() const &
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        Matrix result(this->height_, this->width_);

        #pragma omp parallel for
        for(size_t i = 0; i < this->length_; ++i)
        {
            // real part of log of negative numbers is 0
            result.data_[i] = ((data_[i] >= 0) ? std::log(data_[i]) : 0);
        }

        return result;
    }

    /** Log operator in-place
     *  Applies the log function to all positive elements, replace by zero otherwise
     *  \return Current object overwritten with the result
     */
    Matrix Log() &&
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        #pragma omp parallel for
        for(size_t i = 0; i < this->length_; ++i)
        {
            // real part of log of negative numbers is 0
            data_[i] = ((data_[i] >= 0) ? std::log(data_[i]) : 0);
        }

        return *this;
    }

    /** Shrinkage
    *   Apply the shrinkage algorithm and returns a new instance if constant
    *   \param thresh_factor The threshold factor to be used on the data
    *   \return A new instance containing the result
    */
    Matrix Shrink(const T thresh_factor) const &
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        Matrix result(this->height_, this->width_);

        #pragma omp parallel for
        for( size_t i = 0; i < this->length_; ++i )
        {
            result.data_[i] = data_[i] * std::max(1 - thresh_factor / std::abs((double)data_[i]), 0.0);
        }

        return result;
    }

    /** Shrinkage in-place
    *   Apply the shrinkage algorithm to the current instance
    *   \param thresh_factor The threshold factor to be used on the data
    *   \return A reference to this
    */
    Matrix Shrink(const T thresh_factor) &&
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        #pragma omp parallel for
        for( size_t i = 0; i < this->length_; ++i )
        {
            data_[i] *= (T) std::max(1 - thresh_factor / std::abs((double)data_[i]), 0.0);
        }
        return *this;
    }

    /** Inner product
     *  \param other Object to divide from current object, element-wise
     *  \return The result of type T
     */
    T Inner(const Matrix& other) const
    {
        try
        {
            this->ArgTest(other, add);
        }
        catch (const std::exception&)
        {
            throw;
        }

        T result = 0;
        #pragma omp parallel for reduction(+:result)
        for(size_t i = 0; i < this->length_; ++i)
        {
            result += this->data_[i] * other.data_[i];
        }

        return result;

        // TODO std::inner_product not yet parallelized : change as soon as available
        // std::inner_product( first.Data(), first.Data() + first.Length(), other.Data(), (T) 0 );
    }

    /** Norm
     *  Norm of all elements considered as a one dimensional vector
     * \return The result of type T
     */
    T Norm(const NormType l_norm) const
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        switch(l_norm)
        {
        case one:
        {
            T result = 0;
            #pragma omp parallel for reduction(+:result)
            for(size_t i = 0; i < this->length_; ++i)
            {
                result += std::abs(data_[i]);
            }
            return result;

            // TODO std::accumulate not yet parallelized : change as soon as available
            // std::accumulate(data_, data_ + (this->length_), (T) 0, [](const T acc, const T next){return acc + std::abs(next);});
        }
        case two:
        {
            return std::sqrt(Norm(two_squared));
        }
        case two_squared:
        {
            return Inner(*this);
        }
        case inf:
        {
            T result = 0;
            #pragma omp parallel for reduction(max:result)
            for(size_t i = 0; i < this->length_; ++i)
            {
                if( data_[i] > result )
                {
                    result = data_[i];
                }
            }
            return result;

            // TODO std::max_element not yet parallelized : change as soon as available
            // *std::max_element(data_, data_ + (this->length_));
        }
        default:
        {
            return 0;
        }
        }
    }

    /** Sum
     *  Sum of all elements, considered as a one dimensional vector
     *  \return The result of type T
     */
    T Sum() const
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        T result = 0;
        #pragma omp parallel for reduction(+:result)
        for(size_t i = 0; i < this->length_; ++i)
        {
            result += data_[i];
        }
        return result;

        // TODO std::accumulate not yet parallelized : change as soon as available
        // std::accumulate(data_, data_ + this->length_, (T) 0);
    }

    /** Print
    *   Prints the data to the console
    */
    void Print() const
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        std::cout << std::endl;
        for( size_t i = 0; i < this->height_; ++i )
        {
            std::cout << std::endl;
            for( size_t j = 0; j < this->width_; ++j )
            {
                std::cout << std::setw(10) << data_[i*this->width_ + j] << " ";
            }
        }
        std::cout << std::endl;
    }

    void PrintRefQual() &
    {
        std::cout << "I'm an lvalue !" << std::endl;
    }
    void PrintRefQual() &&
    {
        std::cout << "I'm an rvalue !" << std::endl;
    }
};
} // namespace astroqut

#endif // ASTROQUT_UTILS_MATRIX_HPP
