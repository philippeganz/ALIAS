///
/// \file include/utils/linearop/operator.hpp
/// \brief Operator class header
/// \details Provide operator container with operator-vector operations
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-06-02
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_HPP
#define ASTROQUT_UTILS_OPERATOR_HPP

#include "utils/linearop/matrix.hpp"

#include <functional>
#include <utility>

namespace astroqut
{

template<class T>
class Operator : public LinearOp
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
    {
#ifdef DEBUG
        std::cout << "Operator : Default constructor called" << std::endl;
#endif // DEBUG
    }

    /** Copy constructor
     *  \param other Object to copy from
     */
    Operator(const Operator& other)
        : LinearOp(other.height_, other.width_)
        , data_(other.data_)
        , transposed_(other.transposed_)
    {
#ifdef DEBUG
        std::cout << "Operator : Copy constructor called" << std::endl;
#endif // DEBUG
    }

    /** Move constructor
     *  \param other Object to copy from
     */
    Operator(Operator&& other)
        : LinearOp(other.height_, other.width_)
        , data_(std::move(other.data_))
        , transposed_(other.transposed_)
    {
#ifdef DEBUG
        std::cout << "Operator : Move constructor called" << std::endl;
#endif // DEBUG
    }

    /** Empty constructor
     *  \param height Height of the operator
     *  \param width Width of the operator
     */
    Operator(size_t height, size_t width) noexcept
        : LinearOp(height, width)
        , data_(nullptr, 0, 0)
        , transposed_(false)
    {
#ifdef DEBUG
        std::cout << "Operator : Empty constructor called" << std::endl;
#endif // DEBUG
    }

    /** Full member constructor
     *  \param data Array containing data used by the operator
     *  \param height Height of the operator
     *  \param width Width of the operator
     *  \param transposed If the operator is transposed
     */
    Operator(Matrix<T> data, size_t height, size_t width, bool transposed) noexcept
        : LinearOp(height, width)
        , data_(Matrix<T>(data))
        , transposed_(transposed)
    {
#ifdef DEBUG
        std::cout << "Operator : Full member constructor called" << std::endl;
#endif // DEBUG
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    virtual Operator* Clone() const = 0;

    /** Default destructor
     */
    virtual ~Operator()
    {
#ifdef DEBUG
        std::cout << "Operator : Default destructor called" << std::endl;
#endif // DEBUG
    }

    /** Access data_
     * \return The current value of data_
     */
    const Matrix<T>& Data() const noexcept
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
            throw std::invalid_argument("Operator dimensions must be non-zero!");
        }
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

    /** Copy assignment operator
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

    /** Move assignment operator
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

    virtual Matrix<T> operator*(const Matrix<T>& other) const = 0;
};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_HPP
