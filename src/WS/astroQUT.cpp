///
/// \file src/WS/astroQUT.cpp
/// \brief AstroQUT solver implementation.
/// \author Jairo Diaz <jairo.diaz@unige.ch> 2016-2017
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.5.0
/// \date 2018-07-07
/// \copyright GPL-3.0
///

#include "fista/poisson.hpp"
#include "utils/linearop/operator/astrooperator.hpp"
#include "WS/astroQUT.hpp"

#include <algorithm>
#include <chrono>
#include <random>

#include <omp.h>

namespace astroqut
{
namespace WS
{


static void MC_compute(const Matrix<double>& mu_hat,
                       Matrix<double>& MC_results,
                       const AstroOperator<double>& astro,
                       bool standardize,
                       bool ps,
                       const Parameters<double>& options)
{
    #pragma omp parallel
    {
        Matrix<double> mu_hat_rnd(mu_hat.Height(), mu_hat.Width());
        std::random_device rnd;
        std::default_random_engine generator(rnd() + omp_get_thread_num());

        #pragma omp for schedule(dynamic)
        for(size_t MC_id = 0; MC_id < options.MC_max; ++MC_id)
        {
            if(MC_id % (options.MC_max/100) == 0)
                std::cout << "\r" + std::to_string(1 + 100 * MC_id/options.MC_max) + "/100";

            #pragma omp simd
            for(size_t i = 0; i < mu_hat.Length(); ++i)
            {
                std::poisson_distribution<int> poisson_gen(mu_hat[i]);
                mu_hat_rnd[i] = poisson_gen(generator);
            }

            mu_hat_rnd -= mu_hat;
            mu_hat_rnd /= mu_hat;

            Matrix<double> rdn_result = astro.WtAtBt(mu_hat_rnd, standardize, true, true, ps).Abs();

            #pragma omp simd
            for(size_t i = 0; i < MC_results.Height(); ++i )
                MC_results[i*MC_results.Width() + MC_id] = rdn_result[i];
        }
    }
    std::cout << std::endl;
}

static void BetaZero(const Matrix<double>& image,
                     const AstroOperator<double>& astro,
                     Parameters<double>& options)
{
    std::cout << "Computing beta0..." << std::endl;
    Matrix<double> null_model = astro.BAW(options.x0, false, true, true, false);

    Matrix<double> non_zero_values(image.Height(), image.Width());
    size_t non_zero_values_amount = 0;
    for(size_t i = 0; i < image.Length(); ++i)
        if( ! IsEqual(image[i], 0.0) )
            non_zero_values[non_zero_values_amount++] = image[i];

    // partial sort
    std::nth_element(&non_zero_values[0],
                     &non_zero_values[non_zero_values_amount/2],
                     &non_zero_values[non_zero_values_amount]);
    // get median
    double non_zero_values_median = non_zero_values[non_zero_values_amount/2];

    double null_model_sum = null_model.Sum();

    options.beta0 = non_zero_values_amount*non_zero_values_median/null_model_sum;
    std::cout << "beta0 = " << std::scientific << options.beta0 << std::endl;
}

static void Standardize(const Matrix<double>& mu_hat,
                        const Matrix<double>& background,
                        const AstroOperator<double>& astro,
                        Parameters<double>& options)
{
    std::cout << "Computing standardization matrix..." << std::endl;
    Matrix<double> MC_results(options.pic_size*2, options.MC_max);

    MC_compute(mu_hat, MC_results, astro, false, false, options);

    options.standardize = Matrix<double>(1.0, options.model_size, 1);

    #pragma omp parallel for simd
    for(size_t i = 0; i < options.pic_size*2; ++i)
    {
        // partial sort
        std::nth_element(&MC_results[i*options.MC_max],
                         &MC_results[i*options.MC_max+options.MC_quantile_PF],
                         &MC_results[(i+1)*options.MC_max]);
        // get quantile
        options.standardize[i] = MC_results[i*options.MC_max+options.MC_quantile_PF];
    }

    // max rescale of wavelets and remove high frequency wavelets
    double wavelet_max_value = *std::max_element(&options.standardize[1], &options.standardize[options.pic_size]);
    for(size_t i = 0; i < options.pic_size/2; ++i)
        options.standardize[i] = wavelet_max_value;
    for(size_t i = options.pic_size/2; i < options.pic_size; ++i)
        options.standardize[i] = std::numeric_limits<double>::infinity();

    // spline ok

    // remove point sources from centre of the picture
    for(size_t i = 3*options.pic_size/8; i < 5*options.pic_size/8; ++i)
        for(size_t j = 3*options.pic_size/8; j < 5*options.pic_size/8; ++j)
            options.standardize[options.pic_size*2 + i*options.pic_size + j] = std::numeric_limits<double>::infinity();
    // remove point sources from empty background
    for(size_t i = 0; i < background.Length(); ++i)
        if( IsEqual(background[i], 0.0) )
        {
            options.standardize[options.pic_size*2 + i] = std::numeric_limits<double>::infinity();
            std::cout << i << " put to zero." << std::endl;
        }
}

static void Lambda(const Matrix<double>& mu_hat,
                   AstroOperator<double>& astro,
                   Parameters<double>& options)
{
    std::cout << "Computing lambda and lambdaI..." << std::endl;
    Matrix<double> MC_results(options.model_size, options.MC_max);

    Matrix<double> lambda_standardize(options.standardize);

    lambda_standardize[0] = std::numeric_limits<double>::infinity();

    astro.Standardize(lambda_standardize);

    MC_compute(mu_hat, MC_results, astro, true, true, options);

    Matrix<double> max_values(options.MC_max, 1);
    #pragma omp parallel for simd
    for(size_t MC_id = 0; MC_id < options.MC_max; ++MC_id)
    {
        // compute max value for each MC simulation in wavelet and spline results
        double max_val = -std::numeric_limits<double>::infinity();
        for(size_t i = 0; i < options.pic_size*2; ++i)
            max_val = std::max(max_val, MC_results[i*options.MC_max + MC_id]);
        max_values[MC_id] = max_val;
    }

    // determine the value of lambda
    std::nth_element(&max_values[0],
                     &max_values[options.MC_quantile_PF],
                     &max_values[options.MC_max]);
    double lambda = max_values[options.MC_quantile_PF];
    std::cout << "lambda = " << std::scientific << lambda << std::endl;

    #pragma omp parallel for simd
    for(size_t MC_id = 0; MC_id < options.MC_max; ++MC_id)
    {
        // compute max value for each MC simulation in point source results
        double max_val = -std::numeric_limits<double>::infinity();
        for(size_t i = options.pic_size*2; i < options.model_size; ++i)
            max_val = std::max(max_val, MC_results[i*options.MC_max + MC_id]);
        max_values[MC_id] = max_val;
    }

    // determine the value of lambdaI
    std::nth_element(&max_values[0],
                     &max_values[options.MC_quantile_PS],
                     &max_values[options.MC_max]);
    double lambdaI = max_values[options.MC_quantile_PS];
    std::cout << "lambdaI = " << std::scientific << lambdaI << std::endl;

    double PS_standardize_ratio = lambdaI/lambda;

    for(size_t i = options.pic_size*2; i < options.model_size; ++i)
        options.standardize[i] *= PS_standardize_ratio;

}

Matrix<double> Solve(const Matrix<double>& image,
                     const Matrix<double>& sensitivity,
                     const Matrix<double>& background,
                     Parameters<double>& options )
{
    options.pic_size = (size_t)std::sqrt(image.Length());
    options.model_size = (options.pic_size + 2) * options.pic_size;
    options.x0 = Matrix<double>(0.0, options.model_size, 1);
    options.x0[0] = 1;

    AstroOperator<double> astro(options.pic_size, options.pic_size, options.pic_size/2, sensitivity, Matrix<double>(1, options.model_size, 1), false, options);

    BetaZero(image, astro, options);

    options.x0[0] = options.beta0;
    Matrix<double> I(0.0, options.pic_size*2, 1);
    I[0] = 1.0;
    Matrix<double> u = astro.BAW(I, false, true, true, false);
    Matrix<double> mu_hat = background + u*options.beta0;
    astro.Transpose();

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    Standardize(mu_hat, background, astro, options);

    Lambda(mu_hat, astro, options);

    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time_MC = end-start;

    "data/512_chandra/computeddivx.data" << options.standardize;

    astro.Transpose();
    astro.Standardize(options.standardize);

    fista::poisson::Parameters<double> params;
    params.log_period = 10;
    params.tol = 1.0e-10;
    params.iter_max = options.iter_max;
    params.init_value = Matrix<double>(0.0, options.model_size, 1);
    params.init_value[0] = options.beta0 * options.standardize[0];
    params.indices = Matrix<size_t>(0,1+options.pic_size+options.pic_size*options.pic_size,1);
    for(size_t i = 1; i < 1+options.pic_size+options.pic_size*options.pic_size; ++i)
        params.indices[i] = i + options.pic_size - 1;

    start = std::chrono::high_resolution_clock::now();

    Matrix<double> solution = fista::poisson::Solve(astro, background, image, options.lambda, params);

    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time_FISTA = end-start;

    std::cout << std::defaultfloat << std::endl;
    std::cout << "Time for MC simulations: " << elapsed_time_MC.count()    << " seconds" << std::endl;
    std::cout << "Time for FISTA solver: "   << elapsed_time_FISTA.count() << " seconds" << std::endl;
    std::cout << std::endl;

    return solution;
}

} // namespace WS
} // namespace astroqut
