///
/// \file include/datacontainer.hpp
/// \brief DataContainer class header
/// \details Provide matrix container with multiple matrix operations used throughout the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-04-17
/// \copyright GPL-3.0
///

#ifndef DATACONTAINER_HPP
#define DATACONTAINER_HPP

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
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
             *  \param width Width of the picture
             *  \param height Height of the picture
             */
            DataContainer(unsigned int* data, int width, int height) noexcept;
            /** File constructor
             *  \param file_path Path to the picture file
             */
            DataContainer(const std::string& file_path);
            /** Copy constructor
             *  \param other Object to copy from
             */
            DataContainer(const DataContainer& other) noexcept;

            /** Default destructor */
            virtual ~DataContainer();

            /** Access width_
             * \return The current value of width_
             */
            size_t Width() const noexcept { return width_; }
            /** Set width_
             * \param val New value to set
             */
            void SetWidth(size_t val) noexcept { width_ = val; }

            /** Access height_
             * \return The current value of height_
             */
            size_t Height() const noexcept { return height_; }
            /** Set height_
             * \param val New value to set
             */
            void SetHeight(size_t val) noexcept { height_ = val; }

            /** Access raw_data_
             * \return The current value of raw_data_
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
            void operator+=(const DataContainer& other);

            /** Multiplicative operator
             *  \param other Object to multiply to current object
             */
            DataContainer operator*(const DataContainer& other) const;

            /** In-place multiplicative operator
             *  \param other Object to multiply to current object
             */
            void operator*=(const DataContainer& other);

        protected:

        private:
            size_t width_; //!< Member variable "width_"
            size_t height_; //!< Member variable "height_"
            unsigned int* raw_data_; //!< Member variable "raw_data_"
    };
} // namespace astroqut

#endif // DATACONTAINER_HPP
