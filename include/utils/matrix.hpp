///
/// \file include/utils/matrix.hpp
/// \brief Matrix class header
/// \details Provide matrix container with multiple matrix operations used in the whole project.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.2.0
/// \date 2018-01-02
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_MATRIX_HPP
#define ASTROQUT_UTILS_MATRIX_HPP

#include "utils/linearop.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
//#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
//#include <sstream>
#include <string>

namespace astroqut
{

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
        std::fill( data_, data_ + (this->length_), number );
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
        std::copy(other.data_, other.data_ + (this->length_), data_);
    }

    /** Move constructor
     *  \param other Object to move from
     */
    Matrix(Matrix<T>&& other) noexcept
        : Matrix(other.data_, other.height_, other.width_)
    {
        other.Data(nullptr);
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
    bool IsValid() const
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

    /** Argument test
     *  \param other Other object to test
     *  \return Throws an error message if instance is not valid.
     */
    bool ArgTest(const Matrix& other, ArgTestType type) const
    {
        bool test_result = false;
        switch(type)
        {
        case mult:
            {
                if( IsValid() && other.IsValid() &&
                    this->width_ == other.height_ )
                {
                    test_result = true;
                }
                break;
            }
        case add:
            {
                if( IsValid() && other.IsValid() &&
                    this->height_ == other.height_ &&
                    this->width_ == other.width_ )
                {
                    test_result = true;
                }
                break;
            }
        default:
            {
                break;
            }
        }
        if( test_result )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
        }
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
            return std::equal(data_, data_+(this->length_), other.data_);
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
            // if the array size is not the same, we need to reallocate memory
            if( this->length_ != other.length_ )
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

            std::copy(other.data_, other.data_ + other.length_, data_);
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
            ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(other.height_, other.width_);

        #pragma GCC ivdep
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
            ArgTest(other, add);
        }
        catch (const std::exception&)
        {
            throw;
        }

        #pragma GCC ivdep
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
            ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(other.height_, other.width_);

        #pragma GCC ivdep
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
            ArgTest(other, add);
        }
        catch (const std::exception&)
        {
            throw;
        }

        #pragma GCC ivdep
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
    Matrix operator*(const Matrix& other) const override final
    {
        try
        {
            ArgTest(other, mult);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result((T) 0, this->height_, other.width_); // initialize to zero

        #pragma omp parallel for
        for(size_t i = 0; i < this->height_; ++i)
        {
            for(size_t k = 0; k < this->width_; ++k)
            {
                for(size_t j = 0; j < other.width_; ++j)
                {
                    result.data_[i*other.width_ + j] += data_[i*this->width_ + k] * other.data_[k*other.width_ + j];
                }
            }
        }

        return result;
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

        #pragma GCC ivdep
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

        #pragma GCC ivdep
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
    Matrix<T>& operator*=(const Matrix& other)
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
            ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(this->height_, this->width_);

        #pragma GCC ivdep
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
            ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        #pragma GCC ivdep
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
            ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        Matrix result(this->height_, this->width_);

        #pragma GCC ivdep
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
            ArgTest(other, add);
        }
        catch (const std::exception& e)
        {
            throw e;
        }

        #pragma GCC ivdep
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

        #pragma GCC ivdep
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

        #pragma GCC ivdep
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
			    #pragma GCC ivdep
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
    Matrix Transpose() &&
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
            #pragma omp parallel for
            for ( size_t i = 0; i < this->height_; ++i )
			{
			    #pragma GCC ivdep
				for ( size_t j = i+1 ; j < this->width_; ++j )
				{
					std::swap(data_[i * this->width_ + j], data_[j * this->width_ + i]);
				}
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

        #pragma GCC ivdep
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

        #pragma GCC ivdep
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

        #pragma GCC ivdep
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

        #pragma GCC ivdep
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
            ArgTest(other, add);
        }
        catch (const std::exception&)
        {
            throw;
        }

        return  std::inner_product( data_,
                                    data_ + (this->length_),
                                    other.data_,
                                    (T) 0 );
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
            return std::accumulate(data_,
                                   data_ + (this->length_),
                                   (T) 0,
                                   [](const T acc, const T next){return acc + std::abs(next);});
        }
        case two:
        {
            return std::sqrt(Norm(two_squared));
        }
        case two_squared:
        {
            return std::inner_product(  data_,
                                        data_ + (this->length_),
                                        data_,
                                        (T) 0 );
        }
        case inf:
        {
            return *std::max_element(data_, data_ + (this->length_));
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

        return std::accumulate(data_, data_ + this->length_, (T) 0);
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
