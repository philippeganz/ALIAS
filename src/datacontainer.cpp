///
/// \file src/datacontainer.cpp
/// \brief Implementation of the DataContainer class
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-04-17
/// \copyright GPL-3.0
///

#include "datacontainer.hpp"

astroqut::DataContainer::DataContainer() noexcept
    : width_(0)
    , height_(0)
    , raw_data_ (nullptr)
    {}

astroqut::DataContainer::DataContainer(unsigned int* data, int width, int height) noexcept
    : width_(width)
    , height_(height)
    , raw_data_ (data)
    {}

astroqut::DataContainer::DataContainer(const std::string& file_path)
    : width_(0)
    , height_(0)
    , raw_data_ (nullptr)
{
    std::ifstream file;
    file.exceptions( std::ifstream::failbit | std::ifstream::badbit );
    try {
        file.open(file_path.c_str(), std::ifstream::in);

        file.close();
    }
    catch (const std::ifstream::failure& e)
    {
        std::cerr << "Could not open " << file_path << "! Please verify that the file exists." << std::endl;
        throw;
    }
}

astroqut::DataContainer::DataContainer(const astroqut::DataContainer& other) noexcept
    : width_(other.Width())
    , height_(other.Height())
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
        width_ = other.Width();
        height_ = other.Height();
        raw_data_ = other.RawData();
    }
    return *this;
}

astroqut::DataContainer& astroqut::DataContainer::operator=(DataContainer&& other) noexcept
{
    if (this != &other)
    {
        width_ = other.Width();
        height_ = other.Height();
        raw_data_ = other.RawData();
        other.SetRawData(nullptr);
    }
    return *this;
}

astroqut::DataContainer astroqut::DataContainer::operator+(const astroqut::DataContainer& other) const
{
    unsigned int* result_data = nullptr;

    if(width_ != 0 && width_ == other.Width() &&
       height_ != 0 && height_ == other.Height() &&
       raw_data_ != nullptr && other.RawData() != nullptr)
    {
        try
        {
            result_data = new unsigned int[width_*height_];
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for resulting array!" << std::endl;
            throw;
        }

        unsigned int* other_data = other.RawData();
        #pragma omp parallel for
        for(size_t i = 0; i < height_; ++i)
        {
            #pragma GCC ivdep
            for(size_t j = 0; j < width_; ++j)
            {
                result_data[i*width_ + j] = raw_data_[i*width_ + j] + other_data[i*width_ + j];
            }
        }
    }
    else
    {
        throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
    }
    return astroqut::DataContainer(result_data, width_, height_);
}

void astroqut::DataContainer::operator+=(const DataContainer& other)
{
    if(width_ != 0 && width_ == other.Width() &&
       height_ != 0 && height_ == other.Height() &&
       raw_data_ != nullptr && other.RawData() != nullptr)
    {
        unsigned int* other_data = other.RawData();
        #pragma omp parallel for
        for(size_t i = 0; i < height_; ++i)
        {
            #pragma GCC ivdep
            for(size_t j = 0; j < width_; ++j)
            {
                raw_data_[i*width_ + j] += other_data[i*width_ + j];
            }
        }
    }
    else
    {
        throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
    }
}

astroqut::DataContainer astroqut::DataContainer::operator*(const astroqut::DataContainer& other) const
{
    unsigned int* result_data = nullptr;

    if(width_ != 0 && other.Width() != 0 &&
       height_ != 0 && other.Height() != 0 &&
       width_ == other.Height() &&
       raw_data_ != nullptr && other.RawData() != nullptr)
    {
        try
        {
            result_data = new unsigned int[other.Width()*height_];
        }
        catch (const std::bad_alloc& ba)
        {
            std::cerr << "Could not allocate memory for resulting array!" << std::endl;
            throw;
        }

        unsigned int* other_data = other.RawData();
        size_t other_width = other.Width();
        #pragma omp parallel for
        for(size_t i = 0; i < height_; ++i)
        {
            for(size_t j = 0; j < other_width; ++j)
            {
                result_data[i*other_width + j] = 0;
                #pragma GCC ivdep
                for(size_t k = 0; k < width_; ++k)
                {
                    result_data[i*other_width + j] += raw_data_[i*width_ + k] * other_data[k*other_width + j];
                }
            }
        }
    }
    else
    {
        throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
    }
    return astroqut::DataContainer(result_data, other.Width(), height_);
}

void astroqut::DataContainer::operator*=(const DataContainer& other)
{
    astroqut::DataContainer result;
    try
    {
        result = *this * other;
    }
    catch(const std::exception& e)
    {
        throw;
    }
    delete[] raw_data_;
    *this = std::move(result);
}
