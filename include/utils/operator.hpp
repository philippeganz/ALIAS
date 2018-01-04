///
/// \file include/utils/operator.hpp
/// \brief Operator class header
/// \details Provide operator container with operator-vector operations
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.2.0
/// \date 2018-01-04
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_HPP
#define ASTROQUT_UTILS_OPERATOR_HPP

#include "utils/matrix.hpp"

#include <functional>

namespace astroqut
{

template<class T>
class Operator : public LinearOp<T>
{
protected:
    Matrix<T> data_; //!< Member variable "data_"
    bool transposed_; //!< Member variable "transposed_"

public:
    /** Default constructor
     */
    Operator()
        : data_()
        , transposed_(false)
    {}

    /** Copy constructor
     *  \param other Object to copy from
     */
    Operator(const Operator& other)
        : LinearOp<T>(other.height_, other.width_)
        , data_(other.data_)
        , transposed_(other.transposed_)
    {}

    /** Empty constructor
     *  \param height Height of the operator
     *  \param width Width of the operator
     */
    Operator(size_t height, size_t width) noexcept
        : LinearOp<T>(height, width)
        , data_(nullptr, 0, 0)
        , transposed_(false)
    {}

    /** Full member constructor
     *  \param data Array containing data used by the operator
     *  \param height Height of the operator
     *  \param width Width of the operator
     *  \param transposed If the operator is transposed
     */
    Operator(const Matrix<T>& data, size_t height, size_t width, bool transposed) noexcept
        : LinearOp<T>(height, width)
        , data_(data)
        , transposed_(transposed)
    {}

    /** Clone function
     *  \return A copy of the current instance
     */
    virtual Operator* Clone() const
    {
        return new Operator(*this);
    }

    /** Default destructor
     */
    virtual ~Operator()
    {}

    /** Access data_
     * \return The current value of data_
     */
    Matrix<T> Data() const noexcept
    {
        return data_;
    }
    /** Set data_
     * \param data New value to set
     */
    void Data(const Matrix<T>& data) noexcept
    {
        data_ = data;
    }

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    virtual bool IsValid() const override
    {
        if( this->height_ != 0 &&
            this->width_ != 0 )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument("Operator dimensions must be non-zero and function shall not be nullptr!");
        }
    }

    /** Access transposed_
     * \return The current value of transposed_
     */
    bool IsTransposed() const noexcept
    {
        return transposed_;
    }

    /** Comparison operator equal
     *  \param other Object to compare to
     *  \return True if both object are the same element-wise, False else
     */
    bool operator==(const Operator& other) const noexcept
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
    }

    /** Comparison operator not-equal
     *  \param other Object to compare to
     *  \return False if both object are the same element-wise, True else
     */
    bool operator!=(const Operator& other) const noexcept
    {
        return !(*this == other);
    }

    /** Assignment operator
     *  \param other Object to assign to current object
     *  \return A reference to this
     */
    Operator& operator=(const Operator& other)
    {
        if( this != &other )
        {
            this->height_ = other.height_;
            this->width_ = other.width_;
            this->length_ = other.length_;
            data_ = other.data_;
            transposed_ = other.transposed_;
        }
        return *this;
    }

    /** Move operator
     *  \param other Object to move to current object
     *  \return A reference to this
     */
    Operator& operator=(Operator&& other) noexcept
    {
        if( this != &other )
        {
            this->height_ = other.height_;
            this->width_ = other.width_;
            this->length_ = other.length_;
            std::swap(data_, other.data_);
            transposed_ = other.transposed_;
        }
        return *this;
    }

    virtual Matrix<T> operator*(const Matrix<T>& other) const
    {
        std::cerr << "Do not use the Operator class directly, please use one of its derived classes." << std::endl;

        return Matrix(other);
    }

    /** Transpose in-place
     *   \return A reference to this
     */
    virtual Operator& Transpose()
    {
        std::swap(this->height_, this->width_);
        transposed_ = !transposed_;
        return *this;
    }
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_HPP