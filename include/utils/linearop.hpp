///
/// \file include/utils/linearop.hpp
/// \brief LinearOp class header
/// \details Provide generic linear operator interface.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.2.0
/// \date 2017-12-29
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_LINEAROP_HPP
#define ASTROQUT_UTILS_LINEAROP_HPP

#include <cstddef>

namespace astroqut
{

/** Types of norm currently implemented
 */
enum NormType {one, two, two_squared, inf};

/** Types of argument tests
 */
enum ArgTestType {mult, add};

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
    {}

    /** Full member constructor
     *  \param height Height of the data
     *  \param width Width of the data
     */
    LinearOp(const size_t height, const size_t width) noexcept
        : height_(height)
        , width_(width)
        , length_(height*width)
    {}

    /** Default destructor
     */
    virtual ~LinearOp() = default;

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
    void Height(const size_t val) noexcept
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
    void Width(const size_t val) noexcept
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

};
} // namespace astroqut

#endif // ASTROQUT_UTILS_LINEAROP_HPP
