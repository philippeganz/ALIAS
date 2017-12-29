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

#include "utils/matrix.hpp"

#include <functional>

namespace astroqut
{

template<class T>
class Operator : public LinearOp<T>
{
private:
    std::function<Matrix<T> (Matrix<T>&, Matrix<T>&)> fun_;

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
    Operator(std::function<Matrix<T> (Matrix<T>&, Matrix<T>&)> fun, const size_t height, const size_t width) noexcept
        : LinearOp<T>(height, width)
        , fun_(fun)
    {}

    /** Access fun_
     * \return The current value of fun_
     */
    std::function<Matrix<T> (Matrix<T>&, Matrix<T>&)> Fun() const noexcept
    {
        return fun_;
    }
    /** Set fun_
     * \param val New value to set
     */
    void Fun(std::function<Matrix<T> (Matrix<T>&, Matrix<T>&)> fun) noexcept
    {
        fun_ = fun;
    }

//    Matrix<T> operator*(const Matrix<T>& other) const override final
//    {
//        Matrix<T> result()
//
//        if( ArgTest(other, mult) )
//        {
//            try
//            {
//                result_data = new T[this->height_*other.width_] {}; // init to zero
//            }
//            catch (const std::bad_alloc& ba)
//            {
//                std::cerr << "Could not allocate memory for resulting array!" << std::endl;
//                throw;
//            }
//
//        }
//        return fun_(other, result_data);
//    }
};
} // namespace astroqut

#endif // ASTROQUT_UTILS_OPERATOR_HPP
