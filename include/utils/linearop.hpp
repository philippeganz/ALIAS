///
/// \file include/utils/linearop.hpp
/// \brief LinearOp class header
/// \details Provide generic linear operator base class.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.2.0
/// \date 2018-01-04
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_LINEAROP_HPP
#define ASTROQUT_UTILS_LINEAROP_HPP

#include "const.hpp"

#include <cstddef>
#include <iostream>
#include <stdexcept>

namespace astroqut
{

/** Types of argument tests
 */
enum ArgTestType {mult, add};

/** Forward declaration of the Matrix class
 */
template <class T> class Matrix;

template <class T>
class LinearOp
{
protected:
    size_t height_; //!< Member variable "height_"
    size_t width_; //!< Member variable "width_"
    size_t length_; //!< Member variable "length_"

public:

    /** Default constructor
     *  Create an empty container
     */
    LinearOp() noexcept
        : height_(0)
        , width_(0)
        , length_(0)
    {
#ifdef DEBUG
        std::cout << "LinearOp : Default constructor called" << std::endl;
#endif // DEBUG
    }

    /** Full member constructor
     *  \param height Height of the data
     *  \param width Width of the data
     */
    LinearOp(size_t height, size_t width) noexcept
        : height_(height)
        , width_(width)
        , length_(height*width)
    {
#ifdef DEBUG
        std::cout << "LinearOp : Full member constructor called" << std::endl;
#endif // DEBUG
    }

    /** Default destructor
     */
    virtual ~LinearOp()
    {
#ifdef DEBUG
        std::cout << "LinearOp : Default destructor called" << std::endl;
#endif // DEBUG
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
        length_ = height_*width_;
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
        length_ = height_*width_;
    }

    /** Access length_
     * \return The current value of length_
     */
    size_t Length() const noexcept
    {
        return length_;
    }

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    virtual bool IsValid() const = 0;

    /** Argument test
     *  \param other Other object to test
     *  \return Throws an error message if instance is not valid.
     */
    bool ArgTest(const LinearOp& other, ArgTestType type) const
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

};
} // namespace astroqut

#endif // ASTROQUT_UTILS_LINEAROP_HPP
