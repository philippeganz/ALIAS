///
/// \file include/datacontainer.hpp
/// \brief DataContainer class header
/// \details Provide matrix container with multiple matrix operations used throughout the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-01
/// \copyright GPL-3.0
///

#ifndef DATACONTAINER_HPP
#define DATACONTAINER_HPP

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace astroqut{
    class DataContainer
    {
        public:
            /** Default constructor
             *  Create an empty container
             */
            DataContainer() noexcept;
            /** Full member constructor
             *  \param data 2D array containing the pixels
             *  \param picture_size Size of the picture
             */
            DataContainer(unsigned int* data, const size_t picture_size) noexcept;
            /** File constructor
             *  \param file_path Path to the picture file
             *  \param picture_size Size of the picture
             */
            DataContainer(const std::string& file_path, const size_t picture_size);
            /** Copy constructor
             *  \param other Object to copy from
             */
            DataContainer(const DataContainer& other) noexcept;

            /** Default destructor */
            virtual ~DataContainer();

            /** Access picture_size_
             * \return The current value of picture_size_
             */
            size_t PictureSize() const noexcept { return picture_size_; }
            /** Set picture_size_
             * \param val New value to set
             */
            void SetPictureSize(size_t val) noexcept { picture_size_ = val; }

            /** Access height_
             * \return The current value of height_
             */
            unsigned int* RawData() const noexcept { return raw_data_; }
            /** Set raw_data_
             * \param val New value to set
             */
            void SetRawData(unsigned int* const val) noexcept { raw_data_ = val; }

            /** Assignment operator
             *  \param other Object to assign to current object
             */
            DataContainer& operator=(const DataContainer& other) noexcept;

            /** Move operator
             *  \param other Object to move to current object
             */
            DataContainer& operator=(DataContainer&& other) noexcept;

            /** Additive operator
             *  \param other Object to add to current object
             */
            DataContainer operator+(const DataContainer& other) const;

            /** In-place additive operator
             *  \param other Object to add to current object
             */
            DataContainer& operator+=(const DataContainer& other);

            /** Multiplicative operator
             *  \param other Object to multiply to current object
             */
            DataContainer operator*(const DataContainer& other) const;

            /** In-place multiplicative operator
             *  \param other Object to multiply to current object
             */
            DataContainer& operator*=(const DataContainer& other);

        protected:

        private:
            size_t picture_size_; //!< Member variable "picture_size_"
            unsigned int* raw_data_; //!< Member variable "raw_data_"
    };
} // namespace astroqut

#endif // DATACONTAINER_HPP
