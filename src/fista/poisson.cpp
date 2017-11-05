///
/// \file src/fista.cpp
/// \brief FISTA implementation.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-08-13
/// \copyright GPL-3.0
///

#include "fista/poisson.hpp"

namespace astroqut{
namespace fista{
namespace poisson{

inline double Func( const DataContainer<double>& Axu,
                    const DataContainer<double>& b )
{
    // sum(A*x+u - b.*log(A*x+u))
    return (Axu - (b ->* (Axu.Log()))).Sum();
}

inline DataContainer<double> FuncGrad( const DataContainer<double>& Axu,
                                       const DataContainer<double>& At,
                                       const DataContainer<double>& b )
{
    // A' * ((A*x+u - b) ./ (A*x+u))
    return At * ((Axu - b) / Axu);
}

inline double FLasso( const DataContainer<double>& Axu,
                      const DataContainer<double>& x_woi,
                      const DataContainer<double>& b,
                      const double lambda )
{
    // sum(A*x+u - b.*log(A*x+u)) + lambda * norm(x[-0],1)
    return Func(Axu, b) + lambda*x_woi.Norm(one);
}

inline double FLassoApprox( const DataContainer<double>& Ayu,
                            const DataContainer<double>& At,
                            const DataContainer<double>& x,
                            const DataContainer<double>& x_woi,
                            const DataContainer<double>& y,
                            const DataContainer<double>& b,
                            double lambda,
                            double L )
{
    // sum(A*x+u - b.*log(A*x+u)) + <(x-y), gradfunc(A,y,u,b,w)> + 0.5*L*||x-y||^2  + lambda * norm(x[-0],1)
    return Func(Ayu, b) + (x-y).Inner(FuncGrad(Ayu, At, b)) + 0.5*L*(x-y).Norm(two_squared) + lambda*x_woi.Norm(one);
}

DataContainer<double> solve( const DataContainer<double>& A,
                             const DataContainer<double>& u,
                             const DataContainer<double>& b,
                             const double lambda,
                             const Parameters& options )
{
    std::cout << std::defaultfloat;
    std::cout << std::string(37, '*') << " FISTA " << std::string(36, '*') << std::endl;
    std::cout << "A : " << A.Height() << "x" << A.Width() << " matrix";
    std::cout << ", u : " << u.Height() << " vector";
    std::cout << ", b : " << b.Height() << " vector" << std::endl;
    std::cout << "lambda :" << lambda << ", tol :" << options.tol << std::endl;
    std::cout << std::string(80, '*') << std::endl << std::endl;
    std::cout << " iter" << " | " << "    tol    " << " | " << "   FLasso   " << " | " << "     Lf    " << " | " << "   lambda  " << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::scientific;

    // x and y variables
    DataContainer<double> x(options.init_value);
    if( x.IsEmpty() )
    {
        x.Height(A.Width());
        x.Width(1);
        x.RawData(new double[A.Width()]{});
    }
    DataContainer<double> x_next(x);
    DataContainer<double> x_next_woi(x_next.RawData()+1, x_next.Height()-1, 1); // points to second element of x_next;
    DataContainer<double> y(x);

    // intermediate results
    DataContainer<double> Axu = A*x+u;
    DataContainer<double> Ax_nextu;
    DataContainer<double> Ayu = Axu;
    DataContainer<double> At = A.Transpose();
    double f_lasso_next = 0.0;
    double f_lasso_previous[10]{};
    f_lasso_previous[0] = std::abs(FLasso(Axu, x_next_woi, b, lambda));
    DataContainer<double> grad_current = FuncGrad(Axu, At, b);

    // FISTA variables
    double tol = std::numeric_limits<double>::infinity();
    double Lf = 1.0;
    double eta = 2.0;
    double L_bar = 0.0;
    double t = 1.0;
    double t_next = 1.0; // (1.0 + std::sqrt(1.0 + 4.0 * t * t)) / 2.0
    size_t k = 0;

    // main loop
    for( k = 0; std::abs(tol) > options.tol && k < options.max_iter; ++k )
    {
        // backtracking loop
        double beta = std::numeric_limits<double>::infinity();
        for( int ik = 0; beta > 0; ++ik )
        {
            L_bar = std::pow(eta, ik) * Lf;
            x_next = y - (grad_current/L_bar);
            x_next_woi.RawData(x_next.RawData()+1); // points to second element of new x_next
            x_next_woi.Shrink(lambda/L_bar);
            Ax_nextu = (A*x_next)+u;
            if( Ax_nextu.ContainsNeg() ) // skip function evaluation if we have negative values
            {
                continue;
            }
            f_lasso_next = FLasso(Ax_nextu, x_next_woi, b, lambda);
            beta = f_lasso_next - FLassoApprox(Ayu, At, x_next, x_next_woi, y, b, lambda, L_bar);
        }

        // FISTA step
        y = x_next + (x_next - x) * ((t - 1.0)/t_next);

        // compute tol from previous function value
        double f_lasso_previous_sum = std::accumulate(f_lasso_previous, f_lasso_previous+10, 0.0) / std::min((double) k+1, 10.0);
        tol = std::abs( f_lasso_next - f_lasso_previous_sum ) / f_lasso_previous_sum;

        // actualize values for next iteration
        x = x_next;
        Axu = Ax_nextu;
        Ayu = A*y+u;
        f_lasso_previous[(k+1) % 10] = f_lasso_next;
        grad_current = FuncGrad(Axu, At, b);
        Lf = (k % 100 == 0 ? 1.0 : L_bar / 2.0);
        t = t_next;
        t_next = 1.0; // (1.0 + std::sqrt(1.0 + 4.0 * t * t)) / 2.0

        if( options.log )
        {
            if( (k+1) % (options.log_period * 50) == 0 )
            {
                std::cout << std::endl << " iter" << " | " << "    tol    " << " | " << "   FLasso   " << " | " << "     Lf    " << " | " << "   lambda  " << std::endl;
                std::cout << std::string(80, '-') << std::endl;
            }
            if( (k+1) % options.log_period == 0 )
            {
                std::cout << std::setw(5) << k+1 << " | " << std::setprecision(4) << std::abs(tol) << " | " << std::setw(12) << f_lasso_next << " | " << Lf << " | " << lambda << std::endl;
            }
        }
    }

    std::cout << std::setw(5) << k << " | " << std::setprecision(4) << std::abs(tol) << " | " << std::setw(12) << f_lasso_next << " | " << Lf << " | " << lambda << std::endl << std::endl;

    x_next_woi.RawData(nullptr); // release pointer

    return x;
}

} // namespace poisson
} // namespace fista
} // namespace astroqut
