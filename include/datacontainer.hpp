///
/// \file include/datacontainer.hpp
/// \brief DataContainer class header
/// \details Provide matrix container with multiple matrix operations used throughout the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-29
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_DATACONTAINER_HPP
#define ASTROQUT_DATACONTAINER_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>

namespace astroqut
{

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

/** Types of norm currently implemented
 */
enum NormType {one, two, two_squared, inf};

template <class T> class DataContainer
{
public:
    /** Default constructor
     *  Create an empty container
     */
    DataContainer() noexcept
        : height_(0)
        , width_(0)
        , raw_data_(nullptr)
    {}

    /** Full member constructor
     *  \param data 2D array containing the pixels
     *  \param width Width of the data
     *  \param height Height of the data
     */
    DataContainer( T* const data, const size_t height, const size_t width) noexcept
        : height_(height)
        , width_(width)
        , raw_data_(data)
    {}

    /** Constant number constructor
     *  \param number Number to fill data with
     *  \param width Width of the data
     *  \param height Height of the data
     */
    DataContainer( const T number, const size_t height, const size_t width)
        : height_(height)
        , width_(width)
        , raw_data_(nullptr)
    {
        try
        {
            raw_data_ = new T[height_*width_];
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for new array!" << std::endl;
            throw;
        }
        std::fill( raw_data_, raw_data_ + (height_*width_), number );
    }

    /** File constructor
     *  \param file_path Path to the data file
     *  \param width Width of the data
     *  \param height Height of the data
     */
    DataContainer(std::string& file_path, const size_t height, const size_t width)
        : height_(height)
        , width_(width)
        , raw_data_(nullptr)
    {
        try
        {
            raw_data_ = new T[height_*width_];
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for new array!" << std::endl;
            throw;
        }
        std::ifstream file;
        file.exceptions( std::ifstream::failbit | std::ifstream::badbit );
        try
        {
            file.open(file_path, std::ifstream::in);
        }
        catch (const std::ifstream::failure& e)
        {
            std::cerr << "Could not open " << file_path << "! Please verify that the file exists." << std::endl;
            throw;
        }
        std::string line;
        T number;
        for( size_t i = 0; getline(file, line); ++i )
        {
            std::stringstream line_stream(line);
            for( size_t j = 0; line_stream >> number; ++j )
            {
                raw_data_[i*width_ + j] = number;
            }
        }
        file.close();
    }

    /** Copy constructor
     *  \param other Object to copy from
     */
    DataContainer(const DataContainer<T>& other)
        : height_(other.Height())
        , width_(other.Width())
        , raw_data_(nullptr)
    {
        try
        {
            raw_data_ = new T[height_*width_];
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for new array!" << std::endl;
            throw;
        }
        std::memcpy(raw_data_, other.RawData(), sizeof(T)*height_*width_);
    }

    /** Move constructor
     *  \param other Object to move from
     */
    DataContainer(DataContainer<T>&& other) noexcept
        : height_(other.Height())
        , width_(other.Width())
        , raw_data_(other.RawData())
    {
        other.RawData(nullptr);
    }

    /** Default destructor */
    virtual ~DataContainer()
    {
        if( raw_data_ != nullptr )
        {
            delete[] raw_data_;
        }
        raw_data_ = nullptr;
    }

    /** Access height_
     * \return The current value of height_
     */
    size_t Height() const noexcept
    {
        return height_;
    }
    /** Set height_
     * \param val New value to set
     */
    void Height(size_t val) noexcept
    {
        height_ = val;
    }

    /** Access width_
     * \return The current value of width_
     */
    size_t Width() const noexcept
    {
        return width_;
    }
    /** Set width_
     * \param val New value to set
     */
    void Width(size_t val) noexcept
    {
        width_ = val;
    }

    /** Access raw_data_
     * \return The current value of raw_data_
     */
    T* RawData() const noexcept
    {
        return raw_data_;
    }
    /** Set raw_data_
     * \param val New value to set
     */
    void RawData(T* const data) noexcept
    {
        raw_data_ = data;
    }

    /** Empty test operator
     *  \return True if empty.
     */
    bool IsEmpty() const noexcept
    {
        if( height_ == 0 &&
                width_ == 0 &&
                raw_data_ == nullptr)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    /** Comparison operator equal
     *  \param other Object to compare to
     */
    template <class U> bool operator==(const DataContainer<U>& other) const noexcept
    {
        if( ! std::is_same<T,U>::value ||
                height_ != other.height_ ||
                width_ != other.Width() )
        {
            return false;
        }

        if( raw_data_ == other.RawData() )
        {
            return true;
        }
        else
        {
            return std::equal(raw_data_, raw_data_+(height_*width_), other.RawData());
        }
    }

    /** Comparison operator not-equal
     *  \param other Object to compare to
     */
    template <class U> bool operator!=(const DataContainer<U>& other) const noexcept
    {
        return !(*this == other);
    }

    /** Assignment operator
     *  \param other Object to assign to current object
     */
    DataContainer<T>& operator=(const DataContainer<T>& other) = default;
    template <class U> DataContainer<T>& operator=(const DataContainer<U>& other)
    {
        if( this != &other )
        {
            height_ = other.Height();
            width_ = other.Width();

            try
            {
                raw_data_ = new T[height_*width_];
            }
            catch (const std::bad_alloc& ba)
            {
                std::cerr << "Could not allocate memory for new array!" << std::endl;
                throw;
            }

            // If data types are not the same, we need to cast every value to the new type
            if( std::is_same<T,U>::value )
            {
                std::copy(raw_data_, raw_data_ + (height_*width_), other.RawData());
            }
            else
            {
                #pragma GCC ivdep
                for( size_t iter = 0; iter < height_*width_; ++iter )
                {
                    raw_data_[iter] = (T) other.RawData()[iter];
                }
            }

        }
        return *this;
    }

    /** Move operator
     *  \param other Object to move to current object
     */
    DataContainer<T>& operator=(DataContainer<T>&& other) noexcept
    {
        if (this != &other)
        {
            height_ = other.Height();
            width_ = other.Width();
            raw_data_ = other.RawData();
            other.SetRawData(nullptr);
        }
        return *this;
    }

    /** Additive operator
     *  \param other Object to add to current object
     */
    template <class U> DataContainer<T> operator+(const DataContainer<U>& other) const
    {
        T* result_data = nullptr;

        if( height_ != 0 && height_ == other.Height() &&
            width_ != 0 && width_ == other.Width() &&
            raw_data_ != nullptr && other.RawData() != nullptr )
        {
            try
            {
                result_data = new T[height_*width_];
            }
            catch (const std::bad_alloc& ba)
            {
                std::cerr << "Could not allocate memory for resulting array!" << std::endl;
                throw;
            }

            #pragma GCC ivdep
            for(size_t i = 0; i < height_*width_; ++i)
            {
                result_data[i] = raw_data_[i] + (T) other.RawData()[i];
            }
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
        }
        return DataContainer(result_data, height_, width_);
    }

    /** Additive operator in-place
     *  \param other Object to add to current object
     */
    template <class U> DataContainer<T>& operator+=(const DataContainer<U>& other)
    {
        this = this + other;
        return *this;
    }

    /** Subtractive operator
     *  \param other Object to subtract to current object
     */
    template <class U> DataContainer<T> operator-(const DataContainer<U>& other) const
    {
        T* result_data = nullptr;

        if( height_ != 0 && height_ == other.Height() &&
            width_ != 0 && width_ == other.Width() &&
            raw_data_ != nullptr && other.RawData() != nullptr )
        {
            try
            {
                result_data = new T[height_*width_];
            }
            catch (const std::bad_alloc& ba)
            {
                std::cerr << "Could not allocate memory for resulting array!" << std::endl;
                throw;
            }

            #pragma GCC ivdep
            for(size_t i = 0; i < height_*width_; ++i)
            {
                result_data[i] = raw_data_[i] - (T) other.RawData()[i];
            }
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
        }
        return DataContainer(result_data, height_, width_);
    }

    /** Subtractive operator in-place
     *  \param other Object to subtract to current object
     */
    template <class U> DataContainer<T>& operator-=(const DataContainer<U>& other)
    {
        this = this - other;
        return *this;
    }

    /** Multiplicative operator
     *  \param other Object to multiply to current object
     */
    template <class U> DataContainer<T> operator*(const DataContainer<U>& other) const
    {
        T* result_data = nullptr;

        if( height_ != 0 && other.Height() != 0 &&
            width_ != 0 && other.Width() != 0 &&
            width_ == other.Height() &&
            raw_data_ != nullptr && other.RawData() != nullptr)
        {
            try
            {
                result_data = new T[height_*other.Width()] {}; // init to zero
            }
            catch (const std::bad_alloc& ba)
            {
                std::cerr << "Could not allocate memory for resulting array!" << std::endl;
                throw;
            }

            #pragma omp parallel for
            for(size_t i = 0; i < height_; ++i)
            {
                for(size_t k = 0; k < width_; ++k)
                {
                    for(size_t j = 0; j < other.Width(); ++j)
                    {
                        result_data[i*other.Width() + j] += raw_data_[i*width_ + k] * (T) other.RawData()[k*other.Width() + j];
                    }
                }
            }
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
        }
        return DataContainer(result_data, height_, other.Width());
    }

    /** Multiplicative operator in-place
     *  \param other Object to multiply to current object
     */
    template <class U> DataContainer<T>& operator*=(const DataContainer<U>& other)
    {
        this = this * other;
        return *this;
    }

    /** Element-wise multiply
     *  \param other Object to multiply to current object, element-wise
     */
    template <class U> DataContainer<T> operator->*(const DataContainer<U>& other) const
    {
        T* result_data = nullptr;

        if( height_ != 0 && height_ == other.Height() &&
            width_ != 0 && width_ == other.Width() &&
            raw_data_ != nullptr && other.RawData() != nullptr )
        {
            try
            {
                result_data = new T[height_*width_];
            }
            catch (const std::bad_alloc& ba)
            {
                std::cerr << "Could not allocate memory for resulting array!" << std::endl;
                throw;
            }

            #pragma GCC ivdep
            for(size_t i = 0; i < height_*width_; ++i)
            {
                result_data[i] = raw_data_[i] * (T) other.RawData()[i];
            }
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
        }
        return DataContainer(result_data, height_, width_);
    }

    /** Element-wise divide
     *  \param other Object to divide from current object, element-wise
     */
    template <class U> DataContainer<T> operator/(const DataContainer<U>& other) const
    {
        T* result_data = nullptr;

        if( height_ != 0 && height_ == other.Height() &&
            width_ != 0 && width_ == other.Width() &&
            raw_data_ != nullptr && other.RawData() != nullptr )
        {
            try
            {
                result_data = new T[height_*width_];
            }
            catch (const std::bad_alloc& ba)
            {
                std::cerr << "Could not allocate memory for resulting array!" << std::endl;
                throw;
            }

            #pragma GCC ivdep
            for(size_t i = 0; i < height_*width_; ++i)
            {
                result_data[i] = (T) ((double) raw_data_[i] / (double) other.RawData()[i]);
            }
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
        }
        return DataContainer(result_data, height_, width_);
    }

    /** Inner product
     *  \param other Object to divide from current object, element-wise
     */
    template <class U> T Inner(const DataContainer<U>& other) const noexcept
    {
        return (T) std::inner_product(  raw_data_,
                                        raw_data_ + (height_*width_),
                                        other.RawData(),
                                        (T) 0);
    }

    /** Log
     *  Applies the log function to all elements
     */
    DataContainer<T> Log() const
    {
        T* result_data = nullptr;

        try
        {
            result_data = new T[height_*width_];
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for resulting array!" << std::endl;
            throw;
        }

        #pragma GCC ivdep
        for(size_t i = 0; i < height_*width_; ++i)
        {
            result_data[i] = (T) std::log(raw_data_[i]);
        }

        return DataContainer(result_data, height_, width_);
    }

    /** Norm
     *  Returns the norm of the vector, not implemented yet for matrices.
     */
    T Norm(const NormType l_norm) const noexcept
    {
        switch(l_norm)
        {
        case one:
        {
            return std::accumulate(raw_data_,
                                   raw_data_ + (height_*width_),
                                   (T) 0,
                                   [](T acc, T next){return acc + abs(next);});
        }
        case two:
        {
            return (T) std::sqrt(std::inner_product(raw_data_,
                                                    raw_data_ + (height_*width_),
                                                    raw_data_,
                                                    (T) 0));
        }
        case two_squared:
        {
            return (T) std::inner_product(  raw_data_,
                                            raw_data_ + (height_*width_),
                                            raw_data_,
                                            (T) 0);
        }
        case inf:
        {
            return *std::max_element(raw_data_, raw_data_ + (height_*width_));
        }
        default:
        {
            return (T) 0;
        }
        }
    }

    /** Sum
     *  Returns the sum of the vector, not implemented yet for matrices.
     */
    T Sum() const noexcept
    {
        return std::accumulate(raw_data_, raw_data_ + (height_*width_), (T) 0);
    }

    /** Transpose
    *   Returns the transposed vector/matrix
    */
    DataContainer<T> Transpose() const
    {
        T* transposed_data = nullptr;

        try
        {
            transposed_data = new T[height_*width_]{};
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for new array!" << std::endl;
            throw;
        }

        if( height_ == 1 || width_ == 1 ) // vector
        {
            std::copy(raw_data_, raw_data_ + (height_*width_), transposed_data);
        }
        else // matrix
        {
            #pragma omp parallel for
            for ( size_t i = 0; i < height_; ++i )
			{
			    #pragma GCC ivdep
				for ( size_t j = 0 ; j < width_; ++j )
				{
					transposed_data[ j * height_ + i ] = raw_data_[ i * width_ + j ];
				}
			}
        }
        return DataContainer(transposed_data, width_, height_); // width <--> height
    }

    /** Print
    *   Prints the data to the console
    */
    void Print() const noexcept
    {
        for( size_t i = 0; i < height_; ++i )
        {
            std::cout << std::endl;
            for( size_t j = 0; j < width_; ++j )
            {
                std::cout << raw_data_[i*width_ + j] << " ";
            }
        }
        std::cout << std::endl;
    }

protected:

private:
    size_t height_; //!< Member variable "height_"
    size_t width_; //!< Member variable "width_"
    T* raw_data_; //!< Member variable "raw_data_"
};
} // namespace astroqut

#endif // ASTROQUT_DATACONTAINER_HPP
