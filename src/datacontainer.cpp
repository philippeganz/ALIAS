///
/// \file src/datacontainer.cpp
/// \brief Implementation of the DataContainer class
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-01
/// \copyright GPL-3.0
///

#include "datacontainer.hpp"

astroqut::DataContainer::DataContainer() noexcept
    : picture_size_(0)
    , raw_data_ (nullptr)
    {}

astroqut::DataContainer::DataContainer(unsigned int* data, const size_t picture_size) noexcept
    : picture_size_(picture_size)
    , raw_data_ (data)
    {}

astroqut::DataContainer::DataContainer(const std::string& file_path, const size_t picture_size)
    : picture_size_(picture_size)
    , raw_data_ (nullptr)
{
    try
    {
        raw_data_ = new unsigned int[picture_size_*picture_size_];
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
        std::string line;
        unsigned int number;
        for( size_t i = 0; getline(file, line); ++i )
        {
            std::stringstream line_stream(line);
            for( size_t j = 0; line_stream >> number; ++j )
            {
                raw_data_[i*picture_size_ + j] = number;
            }
        }
    }
    catch (const std::ifstream::failure& e)
    {
        std::cerr << "Could not open " << file_path << "! Please verify that the file exists." << std::endl;
        throw;
    }
    file.close();
}

astroqut::DataContainer::DataContainer(const astroqut::DataContainer& other) noexcept
    : picture_size_(other.PictureSize())
    , raw_data_ (other.RawData())
    {}

astroqut::DataContainer::~DataContainer()
{
    delete[] raw_data_;
    raw_data_ = nullptr;
}

astroqut::DataContainer& astroqut::DataContainer::operator=(const DataContainer& other) noexcept
{
    if (this != &other)
    {
        picture_size_ = other.PictureSize();
        raw_data_ = other.RawData();
    }
    return *this;
}

astroqut::DataContainer& astroqut::DataContainer::operator=(DataContainer&& other) noexcept
{
    if (this != &other)
    {
        picture_size_ = other.PictureSize();
        raw_data_ = other.RawData();
        other.SetRawData(nullptr);
    }
    return *this;
}

astroqut::DataContainer astroqut::DataContainer::operator+(const astroqut::DataContainer& other) const
{
    unsigned int* result_data = nullptr;

    if(picture_size_ != 0 && picture_size_ == other.PictureSize() &&
       raw_data_ != nullptr && other.RawData() != nullptr)
    {
        try
        {
            result_data = new unsigned int[picture_size_*picture_size_];
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for resulting array!" << std::endl;
            throw;
        }

        unsigned int* other_data = other.RawData();
        #pragma omp parallel for
        for(size_t i = 0; i < picture_size_; ++i)
        {
            #pragma GCC ivdep
            for(size_t j = 0; j < picture_size_; ++j)
            {
                result_data[i*picture_size_ + j] = raw_data_[i*picture_size_ + j] + other_data[i*picture_size_ + j];
            }
        }
    }
    else
    {
        throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
    }
    return astroqut::DataContainer(result_data, picture_size_);
}

astroqut::DataContainer& astroqut::DataContainer::operator+=(const DataContainer& other)
{
    if(picture_size_ != 0 && picture_size_ == other.PictureSize() &&
       raw_data_ != nullptr && other.RawData() != nullptr)
    {
        unsigned int* other_data = other.RawData();
        #pragma omp parallel for
        for(size_t i = 0; i < picture_size_; ++i)
        {
            #pragma GCC ivdep
            for(size_t j = 0; j < picture_size_; ++j)
            {
                raw_data_[i*picture_size_ + j] += other_data[i*picture_size_ + j];
            }
        }
        return *this;
    }
    else
    {
        throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
    }
}

astroqut::DataContainer astroqut::DataContainer::operator*(const astroqut::DataContainer& other) const
{
    unsigned int* result_data = nullptr;

    if(picture_size_ != 0 && picture_size_ == other.PictureSize() &&
       raw_data_ != nullptr && other.RawData() != nullptr)
    {
        try
        {
            result_data = new unsigned int[picture_size_*picture_size_];
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for resulting array!" << std::endl;
            throw;
        }

        unsigned int* other_data = other.RawData();
        #pragma omp parallel for
        for(size_t i = 0; i < picture_size_; ++i)
        {
            for(size_t j = 0; j < picture_size_; ++j)
            {
                result_data[i*picture_size_ + j] = 0;
                #pragma GCC ivdep
                for(size_t k = 0; k < picture_size_; ++k)
                {
                    result_data[i*picture_size_ + j] += raw_data_[i*picture_size_ + k] * other_data[k*picture_size_ + j];
                }
            }
        }
    }
    else
    {
        throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
    }
    return astroqut::DataContainer(result_data, picture_size_);
}

astroqut::DataContainer& astroqut::DataContainer::operator*=(const DataContainer& other)
{
    if(picture_size_ != 0 && picture_size_ == other.PictureSize() &&
       raw_data_ != nullptr && other.RawData() != nullptr)
    {
        unsigned int* other_data = other.RawData();
        #pragma omp parallel for
        for(size_t i = 0; i < picture_size_; ++i)
        {
            for(size_t j = 0; j < picture_size_; ++j)
            {
                unsigned int temp = 0;
                #pragma GCC ivdep
                for(size_t k = 0; k < picture_size_; ++k)
                {
                    temp += raw_data_[i*picture_size_ + k] * other_data[k*picture_size_ + j];
                }
                raw_data_[i*picture_size_ + j] = temp;
            }
        }
        return *this;
    }
    else
    {
        throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
    }
}
