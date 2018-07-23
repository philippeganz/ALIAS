///
/// \file include/WS/astroQUT.hpp
/// \brief AstroQUT solver header
/// \author Jairo Diaz <jairo.diaz@unige.ch> MATLAB version 2016-2017
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.5.0
/// \date 2018-07-21
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_WS_ASTROQUT_HPP
#define ASTROQUT_WS_ASTROQUT_HPP

#include "fista/poisson.hpp"
#include "utils/linearop/matrix.hpp"

namespace astroqut
{
namespace WS
{

template<class T>
struct Parameters
{
    /** Default constructor
     *  Parameters with default values
     */
    Parameters()
        : pic_size(0)
        , model_size(0)
        , blur_thresh(0.01)
        , blur_alpha(1.449)
        , blur_R0(2.2364)
        , wavelet{3,8}
        , MC_max(100)
        , MC_quantile_PF(80)
        , MC_quantile_PS(99)
        , beta0(std::numeric_limits<double>::infinity())
        , lambda(1.0)
        , standardize{}
        , x0{}
        , fista_params{}
    {
#ifdef DEBUG
        std::cerr << "WS::Parameters Default constructor called." << std::endl;
#endif // DEBUG
    }

    /** Default destructor
     */
    ~Parameters()
    {
#ifdef DEBUG
        std::cerr << "WS::Parameters Default destructor called." << std::endl;
#endif // DEBUG
    }

    size_t pic_size;
    size_t model_size;
    T blur_thresh; //!< Member variable "blur_thresh" threshold for bluring mask size
    T blur_alpha; //!< Member variable "blur_alpha" alpha in psf
    T blur_R0; //!< Member variable "blur_R0" r0 in psf
    int wavelet[2]; //!< Member variable "wavelet" wavelet type
    size_t MC_max;
    size_t MC_quantile_PF;
    size_t MC_quantile_PS;
    T beta0; //!< Member variable "beta0" intercept
    T lambda;
    Matrix<T> standardize;
    Matrix<T> x0;
    fista::poisson::Parameters<T> fista_params;
};

void Prepare(const Matrix<double>& image,
             const Matrix<double>& sensitivity,
             const Matrix<double>& background,
             Parameters<double>& options);

Matrix<double> Solve(const Matrix<double>& image,
                     const Matrix<double>& sensitivity,
                     const Matrix<double>& background,
                     Parameters<double>& options );

} // namespace WS
} // namespace astroqut

#endif // ASTROQUT_WS_ASTROQUT_HPP

