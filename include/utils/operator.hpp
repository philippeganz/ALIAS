///
/// \file include/utils/operator.hpp
/// \brief Operator class header
/// \details Provide operator container with operator-vector operations
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017
/// \version 0.2.0
/// \date 2017-12-28
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_HPP
#define ASTROQUT_UTILS_OPERATOR_HPP

#include "utils/container.hpp"

#include "functional"

namespace astroqut
{

template<class T>
class Operator : public LinearOp<T>
{
private:
    std::function<Matrix(Matrix)> fun_;

public:
    /** Default constructor
     */
    Operator()
    {}

    /** Default destructor
     */
    virtual ~Operator()
    {}

    /** Copy constructor
     *  \param other Object to copy from
     */
    Operator(const Operator& other)
        : LinearOp<T>(other.Height(), other.Width())
        , fun_(other.Fun())
    {}

    /** Full member constructor
     *  \param fun Function describing the matrix-like operator
     *  \param height Height of the operator
     *  \param width Width of the operator
     */
    Matrix(std::function<Matrix(Matrix)> const fun, const size_t height, const size_t width) noexcept
        : LinearOp<T>(height, width)
        , fun_(fun)
    {}

    /** Access fun_
     * \return The current value of fun_
     */
    std::function<Matrix(Matrix)> Fun() const noexcept
    {
        return fun_
    }
    /** Set fun_
     * \param val New value to set
     */
    void Fun(std::function<Matrix(Matrix)> val) noexcept
    {
        fun_ = val;
    }
};
} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_HPP
