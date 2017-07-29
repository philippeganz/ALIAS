///
/// \file src/fista.cpp
/// \brief FISTA implementation.
/// \author Philippe Ganz <philippe.ganz@gmail.com>
/// \version 0.1.0
/// \date 2017-07-29
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

inline DataContainer<double> FuncGrad(  const DataContainer<double>& Axu,
                                        const DataContainer<double>& At,
                                        const DataContainer<double>& b )
{
    // A' * ((A*x+u - b) ./ (A*x+u))
    return At * ((Axu - b) / Axu);
}

inline double FLasso(   const DataContainer<double>& Axu,
                        const DataContainer<double>& x,
                        const DataContainer<double>& b,
                        const double lambda )
{
    // sum(A*x+u - b.*log(A*x+u)) + lambda * norm(x[-0],1)
    return Func(Axu, b) + lambda*(x.Norm(one) - std::abs(x.RawData()[0]));
}

inline double FLassoApprox( const DataContainer<double>& Ayu,
                            const DataContainer<double>& At,
                            const DataContainer<double>& x,
                            const DataContainer<double>& y,
                            const DataContainer<double>& b,
                            double lambda,
                            double L )
{
    DataContainer<double> x_y = x-y;
    // sum(A*x+u - b.*log(A*x+u)) + <(x-y), gradfunc(A,y,u,b,w)> + 0.5*L*||x-y||^2  + lambda * norm(x[-0],1)
    return Func(Ayu, b) + x_y.Inner(FuncGrad(Ayu, At, b)) + 0.5*L*x_y.Norm(two_squared) + lambda*(x.Norm(one) - std::abs(x.RawData()[0]));
}

DataContainer<double> solve(const DataContainer<double>& model,
                            const DataContainer<double>& background,
                            const DataContainer<double>& response,
                            const double lambda,
                            const double lambda_max,
                            const Parameters& options )
{
    DataContainer<double> x(options.InitValue());
    if( x.IsEmpty() )
    {
        x.Height(model.Width());
        x.Width(1);
        x.RawData(new double[model.Width()]{});
    }

    double tol = std::numeric_limits<double>::infinity();
    double Lf = 1.0;
    double eta = 2.0;
    DataContainer<double> y(x);

    for( size_t k = 0; k < options.MaxIter(), std::abs(tol) > options.Tol(); ++k )
    {

    }

}

} // namespace poisson
} // namespace fista
} // namespace astroqut
