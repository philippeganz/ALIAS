///
/// \file include/datacontainer.hpp
/// \brief DataContainer class header
/// \details Provide matrix container with multiple matrix operations used throughout the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-06
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_DATACONTAINER_HPP
#define ASTROQUT_DATACONTAINER_HPP

#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace astroqut{
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
    DataContainer( const T* const data, const size_t height, const size_t width)
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
        std::memset( raw_data_, number, height_*width_ );
    }

    /** File constructor
     *  \param file_path Path to the data file
     *  \param width Width of the data
     *  \param height Height of the data
     */
    DataContainer(const std::string& file_path, const size_t height, const size_t width)
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
    DataContainer(const astroqut::DataContainer<T>& other)
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
    DataContainer(astroqut::DataContainer<T>&& other) noexcept
        : height_(other.Height())
        , width_(other.Width())
        , raw_data_(other.RawData())
    {
        other.SetRawData(nullptr);
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
    size_t Height() const noexcept { return height_; }
    /** Set height_
     * \param val New value to set
     */
    void Height(size_t val) noexcept { height_ = val; }

    /** Access width_
     * \return The current value of width_
     */
    size_t Width() const noexcept { return width_; }
    /** Set width_
     * \param val New value to set
     */
    void Width(size_t val) noexcept { width_ = val; }

    /** Access raw_data_
     * \return The current value of raw_data_
     */
    T* RawData() const noexcept { return raw_data_; }
    /** Set raw_data_
     * \param val New value to set
     */
    void RawData(T* const data) noexcept { raw_data_ = data; }

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
        if( std::is_same<T,U>::value &&
            raw_data_ == other.RawData() &&
            height_ == other.height_ &&
            width_ == other.Width() )
        {
            return true;
        }
        else
        {
            return false;
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
                std::memcpy(raw_data_, other.RawData(), sizeof(T)*height_*width_);
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
        if( height_ != 0 && height_ == other.Height() &&
            width_ != 0 && width_ == other.Width() &&
            raw_data_ != nullptr && other.RawData() != nullptr )
        {
            #pragma GCC ivdep
            for(size_t i = 0; i < height_*width_; ++i)
            {
                raw_data_[i] += (T) other.RawData()[i];
            }
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
        }
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
                result_data = new T[height_*other.Width()];
            }
            catch (const std::bad_alloc& ba)
            {
                std::cerr << "Could not allocate memory for resulting array!" << std::endl;
                throw;
            }

            #pragma omp parallel for
            for(size_t i = 0; i < height_; ++i)
            {
                for(size_t j = 0; j < other.Width(); ++j)
                {
                    result_data[i*other.Width() + j] = 0;
                    for(size_t k = 0; k < width_; ++k)
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
        return DataContainer(result_data, width_, height_);
    }

    /** Multiplicative operator in-place
     *  \param other Object to multiply to current object
     */
    template <class U> DataContainer<T>& operator*=(const DataContainer<U>& other)
    {
        if( height_ != 0 && other.Height() != 0 &&
            width_ != 0 && other.Width() != 0 &&
            width_ == other.Width() && height_ == other.Height() &&
            raw_data_ != nullptr && other.RawData() != nullptr)
        {
            #pragma omp parallel for
            for(size_t i = 0; i < height_; ++i)
            {
                for(size_t j = 0; j < other.Width(); ++j)
                {
                    T temp = 0;
                    for(size_t k = 0; k < width_; ++k)
                    {
                        temp += raw_data_[i*width_ + k] * (T) other.RawData()[k*other.Width() + j];
                    }
                    raw_data_[i*other.Width() + j] = temp;
                }
            }
            return *this;
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
        }
        return *this;
    }

protected:

private:
    size_t height_; //!< Member variable "height_"
    size_t width_; //!< Member variable "width_"
    T* raw_data_; //!< Member variable "raw_data_"
};
} // namespace astroqut

#endif // ASTROQUT_DATACONTAINER_HPP
