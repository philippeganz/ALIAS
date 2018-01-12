///
/// \file include/utils/matrix.hpp
/// \brief Matrix class header
/// \details Provide matrix container with multiple matrix operations used in the whole project.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-01-12
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_MATRIX_HPP
#define ASTROQUT_UTILS_MATRIX_HPP

#include "utils/linearop.hpp"
#include "settings.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <string>
#include <type_traits>

#pragma GCC system_header
#include <Eigen/Dense>


namespace astroqut
{

/** Types of norm currently implemented
 */
enum NormType {one, two, two_squared, inf};

/** Function to verify if data is aligned
 *  \param ptr Pointer to evaluate
 *  \param align_byte_size The byte boundary size
 */
static inline bool IsAligned(const void* ptr, size_t align_byte_size)
{
    return (uintptr_t)ptr % align_byte_size == 0;
}

/** Floating point type comparison function
 *  \param first First number to compare
 *  \param second Second number to compare
 *  \param error Amount of ULPs to use as threshold for equality
 */
template <class T>
inline typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type IsEqual(T first, T second, T error=1.0)
{
    return std::abs(first-second) < std::abs(first+second)*std::numeric_limits<T>::epsilon()*error ||
        std::abs(first-second) < std::numeric_limits<T>::min();
}

/** Integral type comparison function
 *  \param first First number to compare
 *  \param second Second number to compare
 */
template <class T>
inline typename std::enable_if<std::numeric_limits<T>::is_integer, bool>::type IsEqual(T first, T second)
{
    return first == second;
}



template <class T>
class Matrix : public LinearOp
{
public:
    typedef T matrix_t __attribute__((aligned (sizeof(T))));

private:
    matrix_t* data_; //!< Member variable "data_"

public:
    /** Default constructor
     *  Create an empty container
     */
    Matrix() noexcept
        : data_(nullptr)
    {
#ifdef DEBUG
        std::cout << "Matrix : Default constructor called" << std::endl;
#endif // DEBUG
    }

    /** Empty constructor
     *  Create a container with default-initialized data
     *  \param height Height of the data
     *  \param width Width of the data
     */
    Matrix( size_t height, size_t width)
        : LinearOp(height, width)
        , data_(nullptr)
    {
#ifdef DEBUG
        std::cout << "Matrix : Empty constructor called" << std::endl;
#endif // DEBUG
        if(this->length_ != 0)
        {
            try
            {
                // allocate aligned memory
                data_ = (T*) _mm_malloc (sizeof(T)*this->length_, sizeof(T));
            }
            catch (const std::bad_alloc&)
            {
                std::cerr << "Could not allocate memory for new array!" << std::endl;
                throw;
            }
        }
    }

    /** Full member constructor
     *  \param data Dynamic 2D array containing the pixels
     *  \param height Height of the data
     *  \param width Width of the data
     */
    Matrix( matrix_t* data, size_t height, size_t width)
        : LinearOp(height, width)
        , data_(data)
    {
        if( !IsAligned(data, sizeof(T)) )
        {
            std::cerr << "Please use only " << sizeof(T) << " bytes aligned data.";
            throw;
        }
#ifdef DEBUG
        std::cout << "Matrix : Full member constructor called" << std::endl;
#endif // DEBUG
    }

    /** Full member constructor
     *  \param data Static 2D array containing the pixels
     *  \param length Length of the static array
     *  \param height Height of the data
     *  \param width Width of the data
     */
    Matrix( const T data[], size_t length, size_t height, size_t width)
        : Matrix(height, width)
    {
        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] = data[i];
        }
#ifdef DEBUG
        std::cout << "Matrix : Full member constructor called" << std::endl;
#endif // DEBUG
    }

    /** Constant number constructor
     *  \param number Number to fill data with
     *  \param height Height of the data
     *  \param width Width of the data
     */
    Matrix( T number, size_t height, size_t width)
        : Matrix(height, width)
    {
#ifdef DEBUG
        std::cout << "Matrix : Constant number constructor called" << std::endl;
#endif // DEBUG
        // TODO std::fill not yet parallelized : change as soon as available
        // std::fill( data_, data_ + (this->length_), number );
        #pragma omp parallel for simd
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
//    Matrix(std::string& file_path, size_t height, size_t width)
//        : Matrix(height, width)
//    {
//        std::ifstream file;
//        file.exceptions( std::ifstream::failbit | std::ifstream::badbit );
//        try
//        {
//            file.open(file_path, std::ifstream::in);
//        }
//        catch (const std::ifstream::failure&)
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
    Matrix(const Matrix& other)
        : Matrix(other.height_, other.width_)
    {
#ifdef DEBUG
        std::cout << "Matrix : Copy constructor called" << std::endl;
#endif // DEBUG
        *this = other;
    }

    /** Move constructor
     *  \param other Object to move from
     */
    Matrix(Matrix&& other) noexcept
        : Matrix(other.data_, other.height_, other.width_)
    {
#ifdef DEBUG
        std::cout << "Matrix : Move constructor called" << std::endl;
#endif // DEBUG
        other.Data(nullptr);
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    Matrix Clone() const
    {
        return Matrix(*this);
    }

    /** Default destructor */
    virtual ~Matrix()
    {
#ifdef DEBUG
        std::cout << "Matrix : Destructor called" << std::endl;
#endif // DEBUG
        if( data_ != nullptr )
        {
            // deallocate aligned memory
            _mm_free(data_);
        }
        data_ = nullptr;
    }

    /** Access data_
     * \return The current value of data_
     */
    matrix_t* Data() const noexcept
    {
        return data_;
    }
    /** Set data_
     * \param data New value to set
     */
    void Data(matrix_t* const data)
    {
        if( !IsAligned(data, sizeof(T)) )
        {
            std::cerr << "Please use only " << sizeof(T) << " bytes aligned data.";
            throw;
        }
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

//            TODO std::equal not yet parallelized : change as soon as available
//            std::equal(data_, data_+(this->length_), other.data_);
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
#ifdef DEBUG
        std::cout << "Matrix : Copy assignment operator called" << std::endl;
#endif // DEBUG

        // we need to deallocate data_ if the array size is not the same
        if( this->length_ != other.length_ )
        {
            if( data_ != nullptr )
            {
                // deallocate aligned memory
                _mm_free(data_);
                data_ = nullptr;
            }
            // and we need to reallocate if there is something to store
            if( other.data_ != nullptr )
            {
                try
                {
                    // allocate aligned data
                    data_ = (T*) _mm_malloc (sizeof(T)*this->length_, sizeof(T));
                }
                catch (const std::bad_alloc&)
                {
                    std::cerr << "Could not allocate memory for new array!" << std::endl;
                    throw;
                }
            }
        }

        // copy data if need be
        if( data_ != nullptr && data_ != other.data_ )
        {
            #pragma omp parallel for simd
            for(size_t i = 0; i < other.length_; ++i)
            {
                data_[i] = other.data_[i];
            }

            // TODO std::copy not yet parallelized : change as soon as available
            // std::copy(other.data_, other.data_ + other.length_, data_);
        }

        // finally, update the size values
        this->height_ = other.height_;
        this->width_ = other.width_;
        this->length_ = other.length_;

        return *this;
    }

    /** Move operator
     *  \param other Object to move to current object
     *  \return A reference to this
     */
    Matrix& operator=(Matrix&& other) noexcept
    {
#ifdef DEBUG
        std::cout << "Matrix : Move assignment operator called" << std::endl;
#endif // DEBUG

        this->height_ = other.height_;
        this->width_ = other.width_;
        this->length_ = other.length_;
        std::swap(data_, other.data_);

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
            this->data_[i] += other.data_[i];
        }

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
            this->data_[i] -= other.data_[i];
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

    /** Multiplicative operator
     *  \param other Object to multiply current object with
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

        Matrix result(this->height_, other.width_);

        // other is a vector
        if( other.height_ == 1 || other.width_ == 1 )
        {
            // this is a vector
            if( this->height_ == 1 || this->width_ == 1 )
            {
                std::cerr << "Vector times a vector is a single number, consider using the Inner function instead" << std::endl;
                result.data_[0] = Inner(other);
            }
            // this is a matrix
            else
            {
                MatrixVectorMult(*this, other, result);
            }
        }
        // other is a matrix
        else
        {
            // this is a vector
            if( this->height_ == 1 || this->width_ == 1 )
            {
                VectorMatrixMult(*this, other, result);
            }
            // this is a matrix
            else
            {
                MatrixMatrixMult(*this, other, result);
            }
        }

        return result;
    }

    /** Multiplicative operator with single number in-place
     *  \param number Number to multiply current object with
     *  \return A reference to this
     */
    Matrix& operator*=(T number)
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] *= number;
        }

        return *this;
    }

    /** Element-wise multiply in-place
     *  \param other Object to multiply to current object, element-wise
     *  \return A reference to this
     */
    Matrix& operator&=(const Matrix& other)
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
            data_[i] *= other.data_[i];
        }

        return *this;
    }

    /** Divide operator with single number in-place
     *  \param number Number to divide the current object with
     *  \return A reference to this
     */
    Matrix& operator/=(T number)
    {
        try
        {
            IsValid();
        }
        catch (const std::exception&)
        {
            throw;
        }

        #pragma omp parallel for simd
        for(size_t i = 0; i < this->length_; ++i)
        {
            data_[i] /= number;
        }

        return *this;
    }

    /** Element-wise divide operator in-place
     *  \param other Object to divide from current object, element-wise
     *  \return A reference to this
     */
    Matrix& operator/=(const Matrix& other)
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
            data_[i] /= other.data_[i];
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
    Matrix&& Transpose() &&
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
        return std::move(*this);
    }

    /** Log operator in-place
     *  Applies the log function to all positive elements, replace by zero otherwise
     *  \return Current object overwritten with the result
     */
    Matrix&& Log() &&
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

        return std::move(*this);
    }

    /** Log operator
     *  Applies the log function to all positive elements, replace by zero otherwise
     *  \return A new instance containing the result
     */
    Matrix Log() const &
    {
        return Matrix(*this).Log();
    }

    /** Shrinkage in-place
     *   Apply the shrinkage algorithm
     *   \param thresh_factor The threshold factor to be used on the data
     *   \return A reference to this
     */
    Matrix&& Shrink(double thresh_factor) &&
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
            data_[i] *= std::max(1 - thresh_factor / std::abs((double)data_[i]), 0.0);
        }

        return std::move(*this);
    }

    /** Shrinkage
     *   Apply the shrinkage algorithm
     *   \param thresh_factor The threshold factor to be used on the data
     *   \return A new instance containing the result
     */
    Matrix Shrink(double thresh_factor) const &
    {
        return Matrix(*this).Shrink(thresh_factor);
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

    void PrintRefQual() const &
    {
        std::cout << "I'm an lvalue !" << std::endl;
    }
    void PrintRefQual() &&
    {
        std::cout << "I'm an rvalue !" << std::endl;
    }
};

/** Additive operator, both Matrix are lvalues
 *  \param first Matrix, lvalue ref
 *  \param second Matrix to add to current object, lvalue ref
 *  \return A new instance containing the result
 */
template<class T>
Matrix<T> operator+(const Matrix<T>& first, const Matrix<T>& second)
{
    Matrix<T> result(second);
    return result += first;
}

/** Additive operator, first Matrix is an rvalue
 *  \param first Matrix, rvalue ref
 *  \param second Matrix to add to current object, lvalue ref
 *  \return A reference to second
 */
template<class T>
Matrix<T>&& operator+(Matrix<T>&& first, const Matrix<T>& second)
{
    return std::move(first += second); // if first is temporary, no need to allocate new memory for the result
}

/** Additive operator, second Matrix is an rvalue
 *  \param first Matrix, lvalue ref
 *  \param second Matrix to add to current object, rvalue ref
 *  \return A reference to second
 */
template<class T>
Matrix<T>&& operator+(const Matrix<T>& first, Matrix<T>&& second)
{
    return std::move(second += first); // if second is temporary, no need to allocate new memory for the result
}

/** Subtractive operator, both Matrix are lvalues
 *  \param first Matrix, lvalue ref
 *  \param second Matrix to subtract from current object, lvalue ref
 *  \return A new instance containing the result
 */
template<class T>
Matrix<T> operator-(const Matrix<T>& first, const Matrix<T>& second)
{
    Matrix<T> result(first);
    return result -= second;
}

/** Subtractive operator, first Matrix is an rvalue, second an lvalue
 *  \param first Matrix, rvalue ref
 *  \param second Matrix to subtract from current object, lvalue ref
 *  \return A reference to first
 */
template<class T>
Matrix<T>&& operator-(Matrix<T>&& first, const Matrix<T>& second)
{
    return std::move(first -= second); // if first is temporary, no need to allocate new memory for the result
}

/** Subtractive operator, first Matrix is an lvalue, second an rvalue
 *  \param first Matrix, lvalue ref
 *  \param second Matrix to subtract from current object, rvalue ref
 *  \return A reference to second
 */
template<class T>
Matrix<T>&& operator-(const Matrix<T>& first, Matrix<T>&& second)
{
    // Need to implement it fully again since `-` is not commutative
    try
    {
        first.ArgTest(second, add);
    }
    catch (const std::exception&)
    {
        throw;
    }

    #pragma omp parallel for simd
    for(size_t i = 0; i < first.Length(); ++i)
    {
        second.Data()[i] = first.Data()[i] - second.Data()[i];
    }

    return std::move(second); // if second is temporary, no need to allocate new memory for the result
}

/** Multiplicative operator element-wise
 *  \param first Matrix, lvalue ref
 *  \param second Matrix to multiply the current object with
 *  \return A new instance containing the result
 */
template <class T>
Matrix<T> operator&(const Matrix<T>& first, const Matrix<T>& second)
{
    Matrix<T> result(first);
    return result &= second;
}

/** Multiplicative operator element-wise of temporary instance
 *  \param first Matrix, rvalue ref
 *  \param second Matrix to multiply the current object with, lvalue ref
 *  \return A reference to first
 */
template <class T>
Matrix<T>&& operator&(Matrix<T>&& first, const Matrix<T>& second)
{
    return std::move(first &= second);
}

/** Multiplicative operator element-wise of temporary instance
 *  \param first Matrix, lvalue ref
 *  \param second Matrix to multiply the current object with, rvalue ref
 *  \return A reference to second
 */
template <class T>
Matrix<T>&& operator&(const Matrix<T>& first, Matrix<T>&& second)
{
    return std::move(second &= first);
}

/** Multiplicative operator with single number
 *  \param mat Matrix, lvalue ref
 *  \param number Number to multiply the current object with
 *  \return A new instance containing the result
 */
template<class T>
Matrix<T> operator*(const Matrix<T>& mat, T number)
{
    Matrix<T> result(mat);
    return result *= number;
}

/** Multiplicative operator with single number of temporary instance
 *  \param mat Matrix, rvalue ref
 *  \param number Number to multiply the current object with
 *  \return A reference to this
 */
template<class T>
Matrix<T>&& operator*(Matrix<T>&& mat, T number)
{
    return std::move(mat *= number);
}

/** Divide operator with single number
 *  \param mat Matrix, lvalue ref
 *  \param number Number to divide the current object with
 *  \return A new instance containing the result
 */
template <class T>
Matrix<T> operator/(const Matrix<T>& mat, T number)
{
    Matrix<T> result(mat);
    return result /= number;
}

/** Divide operator with single number of temporary instance
 *  \param mat Matrix, rvalue ref
 *  \param number Number to divide the current object with
 *  \return A reference to mat
 */
template <class T>
Matrix<T>&& operator/(const Matrix<T>&& first, T number)
{
    return std::move(first /= number);
}

/** Divide operator element-wise
 *  \param first Matrix, lvalue ref
 *  \param second Matrix to divide the current object with
 *  \return A new instance containing the result
 */
template <class T>
Matrix<T> operator/(const Matrix<T>& first, const Matrix<T>& second)
{
    Matrix<T> result(first);
    return result /= second;
}

/** Divide operator element-wise of temporary instance
 *  \param first Matrix, rvalue ref
 *  \param second Matrix to divide the current object with
 *  \return A reference to first
 */
template <class T>
Matrix<T>&& operator/(Matrix<T>&& first, const Matrix<T>& second)
{
    return std::move(first /= second);
}


using namespace settings;

/** Matrix Matrix multiplication
 *  Performs a matrix-matrix multiplication : result = first * second
 *  \param first First matrix of size l by m
 *  \param second Second matrix of size m by n
 *  \param result Resulting matrix of size l by n
 */
template <class T>
inline void MatrixMatrixMult(const Matrix<T>& first, const Matrix<T>& second, const Matrix<T>& result, MMMultType type = default_MMType)
{
    switch(type)
    {
    // Naive implementation, making use of multi-threading when available
    case MM_naive:
        {
            // Init result to zero
            for(size_t i = 0; i < result.Length(); ++i)
            {
                result.Data()[i] = 0;
            }
            #pragma omp parallel for
            for(size_t i = 0; i < first.Height(); ++i)
            {
                for(size_t k = 0; k < first.Width(); ++k)
                {
                    for(size_t j = 0; j < second.Width(); ++j)
                    {
                        result.Data()[i*second.Width() + j] += first.Data()[i*first.Width() + k] * second.Data()[k*second.Width() + j];
                    }
                }
            }

            break;
        }
    case MM_pure_eigen:
        {
            MatrixMatrixMultEigen(first, second, result);

            break;
        }
    default:
        {}
    }
}

/** Matrix Matrix multiplication using the Eigen library
 *  \param first First matrix of size l by m
 *  \param second Second matrix of size m by n
 *  \param result Resulting matrix of size l by n
 */
template <class T>
inline void MatrixMatrixMultEigen(const Matrix<T>& first, const Matrix<T>& second, const Matrix<T>& result)
{
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> eigenMat(first.Data(), first.Height(), first.Width());
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> eigenVect(second.Data(), second.Height(), second.Width());
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> eigenResult(result.Data(), result.Height(), result.Width());
    eigenResult = eigenMat * eigenVect;
}
template <>
inline void MatrixMatrixMultEigen(const Matrix<double>& first, const Matrix<double>& second, const Matrix<double>& result)
{
    // loading the data with 8 bytes alignment
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenMat(first.Data(), first.Height(), first.Width());
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenVect(second.Data(), second.Height(), second.Width());
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenResult(result.Data(), result.Height(), result.Width());
    eigenResult = eigenMat * eigenVect;
}
template <>
inline void MatrixMatrixMultEigen(const Matrix<long long>& first, const Matrix<long long>& second, const Matrix<long long>& result)
{
    // loading the data with 8 bytes alignment
    Eigen::Map<Eigen::Matrix<long long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenMat(first.Data(), first.Height(), first.Width());
    Eigen::Map<Eigen::Matrix<long long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenVect(second.Data(), second.Height(), second.Width());
    Eigen::Map<Eigen::Matrix<long long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenResult(result.Data(), result.Height(), result.Width());
    eigenResult = eigenMat * eigenVect;
}


/** Matrix Vector multiplication
 *  Performs a matrix-vector multiplication : result = mat * vect
 *  \param mat Matrix of size m by n
 *  \param vect Vector of size n by 1
 *  \param result Resulting vector of size m by 1
 *  \param type Computation strategy to use
 */
template <class T>
inline void MatrixVectorMult(const Matrix<T>& mat, const Matrix<T>& vect, const Matrix<T>& result, MVMultType type = default_MVType)
{
    switch(type)
    {
    // Naive implementation, making use of multi-threading when available
    case MV_naive:
        {
            #pragma omp parallel for
            for(size_t i = 0; i < mat.Height(); ++i)
            {
                result.Data()[i] = 0;
                for(size_t j = 0; j < mat.Width(); ++j)
                {
                    result.Data()[i] += mat.Data()[i*mat.Width() + j] * vect.Data()[j];
                }
            }

            break;
        }
    // Eigen implementation, that does currently not make use of multi-threading
    case MV_pure_eigen:
        {
            MatrixMatrixMultEigen(mat, vect, result);

            break;
        }
    // Mixed Eigen-OpenMP implementation, that does make use of multi-threading
    case MV_par_eigen:
        {
            MatrixVectMultEigen(mat, vect, result);

            break;
        }
    default:
        {}
    }
}

/** Matrix Vector multiplication using the Eigen library
 *  \param mat Matrix of size m by n
 *  \param vect Vector of size n by 1
 *  \param result Resulting vector of size m by 1
 */
template <class T>
inline void MatrixVectMultEigen(const Matrix<T>& mat, const Matrix<T>& vect, const Matrix<T>& result)
{
    Eigen::Map<Eigen::Matrix<T, 1, Eigen::Dynamic, Eigen::RowMajor>> eigenVect(vect.Data(), vect.Height());
    #pragma omp parallel for
    for(size_t i = 0; i < mat.Height(); ++i)
    {
        Eigen::Map<Eigen::Matrix<T, 1, Eigen::Dynamic, Eigen::RowMajor>> eigenMat(mat.Data() + i*mat.Width(), mat.Width());
        result.Data()[i] = eigenMat.dot(eigenVect);
    }
}
template <>
inline void MatrixVectMultEigen(const Matrix<double>& mat, const Matrix<double>& vect, const Matrix<double>& result)
{
    Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenVect(vect.Data(), vect.Height());
    #pragma omp parallel for
    for(size_t i = 0; i < mat.Height(); ++i)
    {
        Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenMat(mat.Data() + i*mat.Width(), mat.Width());
        result.Data()[i] = eigenMat.dot(eigenVect);
    }
}
template <>
inline void MatrixVectMultEigen(const Matrix<long long>& mat, const Matrix<long long>& vect, const Matrix<long long>& result)
{
    Eigen::Map<Eigen::Matrix<long long, 1, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenVect(vect.Data(), vect.Height());
    #pragma omp parallel for
    for(size_t i = 0; i < mat.Height(); ++i)
    {
        Eigen::Map<Eigen::Matrix<long long, 1, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Aligned8> eigenMat(mat.Data() + i*mat.Width(), mat.Width());
        result.Data()[i] = eigenMat.dot(eigenVect);
    }
}


/** Vector Matrix multiplication
 *  Performs a matrix-vector multiplication : result = vect * mat
 *  \param vect Vector of size 1 by m
 *  \param mat Matrix of size m by n
 *  \param result Resulting vector of size 1 by n
 */
template <class T>
inline void VectorMatrixMult(const Matrix<T>& vect, const Matrix<T>& mat, const Matrix<T>& result, VMMultType type = default_VMType)
{
    switch(type)
    {
    // Naive implementation, making use of multi-threading when available
    case VM_naive:
        {
            // Init result to zero
            for(size_t i = 0; i < result.Length(); ++i)
            {
                result.Data()[i] = 0;
            }
            #pragma omp parallel for
            for(size_t i = 0; i < mat.Height(); ++i)
            {
                for(size_t j = 0; j < mat.Width(); ++j)
                {
                    result.Data()[j] += vect.Data()[i] * mat.Data()[i*mat.Width() + j];
                }
            }

            break;
        }
    case VM_pure_eigen:
        {
            MatrixMatrixMultEigen(vect, mat, result);

            break;
        }
    default:
        {}
    }
}

} // namespace astroqut

#endif // ASTROQUT_UTILS_MATRIX_HPP
