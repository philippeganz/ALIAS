///
/// \file include/fista/poisson.hpp
/// \brief FISTA (Fast Iterative Shrinkage Tresholding Algorithm) solver for Poisson distributed noise.
/// \author Hatef Monajemi <monajemi@stanford.edu> 2012-2014
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.4.0
/// \date 2018-02-25
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_FISTA_POISSON_HPP
#define ASTROQUT_FISTA_POISSON_HPP

#include "utils/linearop/matrix.hpp"
#include "utils/linearop/operator/matmult.hpp"
#include "utils/linearop/operator/convolution.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>

namespace astroqut{
namespace fista{
namespace poisson{

template<class T>
struct Parameters
{
    /** Default constructor
     *  Parameters with default values
     */
    Parameters() noexcept
        : tol(1e-6)
        , max_iter(2000)
        , init_value{}
        , log(true)
        , log_period(10)
    {}

    T tol; //!< Member variable "tol"
    size_t max_iter = 2000; //!< Member variable "max_iter"
    Matrix<T> init_value; //!< Member variable "init_value"
    bool log; //!< Member variable "log"
    unsigned int log_period; //!< Member variable "log_period"
};

/** Poisson distributed noise solver
 *  \param A Explicit regression matrix
 *  \param u Background shift
 *  \param b Response data to the regression matrix
 *  \param lambda Regularization parameter
 *  \param options Parameters that defines various value for FISTA to work
 */
template<class T>
T Func(const Matrix<T>& Axu,
       const Matrix<T>& b )
{
    // sum(A*x+u - b.*log(A*x+u))
    return (Axu - (b & (Axu.Log()))).Sum();
}
template<class T>
Matrix<T> FuncGrad(const Matrix<T>& Axu,
                   const Operator<T>& At,
                   const Matrix<T>& b )
{
    // A' * ((A*x+u - b) ./ (A*x+u))
    return At * ((Axu - b) / Axu);
}
template<class T>
T FLasso(const Matrix<T>& Axu,
         const Matrix<T>& x_woi,
         const Matrix<T>& b,
         T lambda )
{
    // sum(A*x+u - b.*log(A*x+u)) + lambda * norm(x[-0],1)
    return Func(Axu, b) + lambda*x_woi.Norm(one);
}
template<class T>
T FLassoApprox(const Matrix<T>& Ayu,
               const Operator<T>& At,
               const Matrix<T>& x,
               const Matrix<T>& x_woi,
               const Matrix<T>& y,
               const Matrix<T>& b,
               T lambda,
               T L )
{
    // sum(A*x+u - b.*log(A*x+u)) + <(x-y), gradfunc(A,y,u,b,w)> + 0.5*L*||x-y||^2  + lambda * norm(x[-0],1)
    return Func(Ayu, b) + Inner( x-y, FuncGrad(Ayu, At, b) ) + 0.5*L*(x-y).Norm(two_squared) + lambda*x_woi.Norm(one);
}
template<class T>
Matrix<T> Solve(const Operator<T>& A,
                const Matrix<T>& u,
                const Matrix<T>& b,
                T lambda,
                const Parameters<T>& options )
{
    std::cout << std::defaultfloat;
    std::cout << std::string(37, '*') << " FISTA " << std::string(36, '*') << std::endl;
    std::cout << "A : " << A.Height() << "x" << A.Width() << " matrix";
    std::cout << ", u : " << u.Height() << " vector";
    std::cout << ", b : " << b.Height() << " vector" << std::endl;
    std::cout << "lambda :" << lambda << ", tol :" << options.tol << std::endl;
    std::cout << std::string(80, '*') << std::endl << std::endl;
    std::cout << " iter" << " | " << "         tol        " << " | " << "       FLasso       " << " | " << "     Lf      " << " | " << " lambda " << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::scientific;

    // x and y variables
    Matrix<T> x(0.0L, A.Width(), 1);
    if( !options.init_value.IsEmpty() )
        x = options.init_value;
    Matrix<T> x_next(x);
    Matrix<T> x_next_woi(x_next.Data()+1, x_next.Height()-1, 1); // points to second element of x_next;
    Matrix<T> y(x);

    // intermediate results
    Matrix<T> Axu = A*x+u;
    Matrix<T> Ax_nextu;
    Matrix<T> Ayu = Axu;
    Operator<T> &At = A.Clone()->Transpose();
    T f_lasso_next = 0.0L;
    T f_lasso_previous[10]{};
    f_lasso_previous[0] = std::abs(FLasso(Axu, x_next_woi, b, lambda));
    Matrix<T> grad_current = FuncGrad(Axu, At, b);

    // FISTA variables
    T tol = std::numeric_limits<T>::infinity();
    T Lf = 1.0L;
    T eta = 2.0L;
    T L_bar = 0.0L;
//    T t = 1.0L;
//    T t_next = (1.0L + std::sqrt(1.0L + 4.0L * t * t)) / 2.0L;
    size_t k = 0;

    // main loop
    while( std::abs(tol) > options.tol && k < options.max_iter )
    {
        // backtracking loop
        T beta = std::numeric_limits<T>::infinity();
        for( int ik = 0; beta > 0; ++ik )
        {
            L_bar = std::pow(eta, ik) * Lf;
            x_next = y - (grad_current/L_bar);
            if(x_next[0] < 0.0)
                x_next[0] = 0;
            x_next_woi.Data(x_next.Data()+1); // points to second element of new x_next
            std::move(x_next_woi).Shrink(lambda/L_bar); //cast to an rvalue to allow in-place shrinkage
            Ax_nextu = (A*x_next)+u;
            if( Ax_nextu.ContainsNeg() ) // skip function evaluation if we have negative values
                continue;
            f_lasso_next = FLasso(Ax_nextu, x_next_woi, b, lambda);
            beta = f_lasso_next - FLassoApprox(Ayu, At, x_next, x_next_woi, y, b, lambda, L_bar);
        }

        // FISTA step
        y = x_next;
//        y = x_next + (x_next - x) * ((t - 1.0)/t_next);

        // compute tol from previous function value
        T f_lasso_previous_sum = std::accumulate(f_lasso_previous, f_lasso_previous+10, 0.0L) / std::min((T) k+1, 10.0L);
        tol = std::abs( f_lasso_next - f_lasso_previous_sum ) / f_lasso_previous_sum;

        // actualize values for next iteration
        ++k;
        x = x_next;
        Axu = Ax_nextu;
        Ayu = A*y+u;
        f_lasso_previous[k % 10] = f_lasso_next;
        grad_current = FuncGrad(Axu, At, b);
        Lf = (k % 100 == 0 ? 1.0L : L_bar / 2.0L);
//        t = t_next;
//        t_next = (1.0L + std::sqrt(1.0L + 4.0L * t * t)) / 2.0L;

        if( options.log )
        {
            if( k % (options.log_period * 50) == 0 )
            {
                std::cout << std::endl << " iter" << " | " << "         tol        " << " | " << "       FLasso       " << " | " << "     Lf      " << " | " << " lambda " << std::endl;
                std::cout << std::string(80, '-') << std::endl;
            }
            if( k % options.log_period == 0 )
            {
                std::cout << std::setw(5) << k << " | " << std::scientific << std::setprecision(10) << std::setw(20) << std::abs(tol) << " | " << std::setw(20) << f_lasso_next << " | " << std::defaultfloat << std::setw(13) << Lf << " | " << std::setw(8) << lambda << std::endl;
            }
        }
    }

    std::cout << std::setw(5) << k << " | " << std::scientific << std::setprecision(10) << std::setw(20) << std::abs(tol) << " | " << std::setw(20) << f_lasso_next << " | " << std::defaultfloat << std::setw(13) << Lf << " | " << std::setw(8) << lambda << std::endl;

    std::cout << "FISTA: converged in " << k << " iterations" << std::endl;
    std::cout << "FISTA: Relative error: " << std::abs(tol) << std::endl;

    x_next_woi.Data(nullptr); // release pointer

    return x;
}

} // namespace poisson
} // namespace fista
} // namespace astroqut

#endif // ASTROQUT_FISTA_POISSON_HPP

