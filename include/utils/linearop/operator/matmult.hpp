///
/// \file include/utils/linearop/operator/matmult.hpp
/// \brief Matrix Multiplication class header
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.6.0
/// \date 2019-01-19
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_MATMULT_HPP
#define ASTROQUT_UTILS_OPERATOR_MATMULT_HPP

#include "utils/linearop/operator.hpp"

namespace astroqut
{

template <class T = double>
class MatMult : public Operator<T>
{
public:

    /** Default constructor
     */
    MatMult()
        : Operator<T>(0, 0)
    {
#ifdef DEBUG
        std::cout << "MatMult : Default constructor called" << std::endl;
#endif // DEBUG
    }

    /** Full member constructor
     *  \param data Array containing data used by the operator
     *  \param height Height of the operator
     *  \param width Width of the operator
     */
    MatMult(Matrix<T> data, size_t height, size_t width)
        : Operator<T>(data, height, width, false)
    {
#ifdef DEBUG
        std::cout << "MatMult : Full member constructor called" << std::endl;
#endif // DEBUG
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    virtual MatMult* Clone() const override
    {
        return new MatMult(*this);
    }

    /** Default destructor
     */
    virtual ~MatMult()
    {
#ifdef DEBUG
        std::cout << "MatMult : Default destructor called" << std::endl;
#endif // DEBUG
    }

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    bool IsValid() const override final
    {
        if( this->height_ != 0 &&
            this->width_ != 0 &&
            !this->data_.IsEmpty() )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument("Operator dimensions must be non-zero and function shall not be nullptr!");
        }
    }

    Matrix<T> operator*(const Matrix<T>& other) const override final
    {
#ifdef DO_ARGCHECKS
        try
        {
            this->ArgTest(other, mult);
        }
        catch (const std::exception&)
        {
            throw;
        }
#endif // DO_ARGCHECKS

        return std::move(this->data_ * other);
    }

    /** Transpose in-place
     *   \return A reference to this
     */
    virtual MatMult& Transpose() override
    {
        std::swap(this->height_, this->width_);
        this->transposed_ = !this->transposed_;
        this->data_ = std::move(this->data_.Transpose());
        return *this;
    }

};

} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_MATMULT_HPP
