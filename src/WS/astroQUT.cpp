///
/// \file src/WS/astroQUT.cpp
/// \brief AstroQUT solver implementation.
/// \author Jairo Diaz <jairo.diaz@unige.ch> 2016-2017
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.5.0
/// \date 2018-08-05
/// \copyright GPL-3.0
///

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


static Matrix<double> MCCompute(const Matrix<double>& mu_hat,
                                std::default_random_engine generator)
{
    Matrix<double> mu_hat_rnd(mu_hat.Height(), mu_hat.Width());

    #pragma omp simd
    for(size_t i = 0; i < mu_hat.Length(); ++i)
    {
        std::poisson_distribution<int> poisson_gen(mu_hat[i]);
        mu_hat_rnd[i] = poisson_gen(generator);
    }

    mu_hat_rnd -= mu_hat;
    mu_hat_rnd /= mu_hat;

    return mu_hat_rnd;
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

static void Standardize(const Matrix<double> mu_hat,
                        const Matrix<double>& background,
                        const AstroOperator<double>& astro,
                        Parameters<double>& options)
{
    std::cout << "Computing standardization matrix..." << std::endl;
    Matrix<double> MC_astro(options.pic_size*2, options.MC_max);

    std::random_device rnd;
    #pragma omp parallel for schedule(dynamic)
    for(size_t MC_id = 0; MC_id < options.MC_max; ++MC_id)
    {
        std::default_random_engine generator(rnd() + omp_get_team_num());
        if(MC_id % (options.MC_max/100) == 0)
            std::cout << "\r" + std::to_string(1 + 100 * MC_id/options.MC_max) + "/100";

        Matrix<double> rnd_result = astro.WtAtBt(MCCompute(mu_hat, generator), false, true, true, false).Abs();

        #pragma omp simd
        for(size_t i = 0; i < rnd_result.Height(); ++i )
            MC_astro[i*MC_astro.Width() + MC_id] = rnd_result[i];
    }
    std::cout << std::endl;

    options.standardize = Matrix<double>(1.0, options.model_size, 1);

    #pragma omp parallel for simd
    for(size_t i = 0; i < options.pic_size*2; ++i)
    {
        // partial sort
        std::nth_element(&MC_astro[i*options.MC_max],
                         &MC_astro[i*options.MC_max+options.MC_quantile_PF],
                         &MC_astro[(i+1)*options.MC_max]);
        // get quantile
        options.standardize[i] = MC_astro[i*options.MC_max+options.MC_quantile_PF];
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
#ifdef VERBOSE
            std::cout << i << " put to zero." << std::endl;
#endif // VERBOSE
        }
}

static void Lambda(const Matrix<double> mu_hat,
                   AstroOperator<double>& astro,
                   Parameters<double>& options)
{
    std::cout << "Computing lambda and lambdaI..." << std::endl;
    Matrix<double> lambda_standardize(options.standardize);
    lambda_standardize[0] = std::numeric_limits<double>::infinity();
    astro.Standardize(lambda_standardize);

    Matrix<double> WS_max_values(-std::numeric_limits<double>::infinity(), options.MC_max, 1);
    Matrix<double> PS_max_values(-std::numeric_limits<double>::infinity(), options.MC_max, 1);

    std::random_device rnd;
    #pragma omp parallel for schedule(dynamic)
    for(size_t MC_id = 0; MC_id < options.MC_max; ++MC_id)
    {

        std::default_random_engine generator(rnd() + omp_get_thread_num());

        if(MC_id % (options.MC_max/100) == 0)
            std::cout << "\r" + std::to_string(1 + 100 * MC_id/options.MC_max) + "/100";

        Matrix<double> rnd_result = astro.WtAtBt(MCCompute(mu_hat, generator)).Abs();

        // compute max value for each MC simulation in wavelet and spline results
        WS_max_values[MC_id] = *std::max_element(&rnd_result[0], &rnd_result[options.pic_size*2]);

        // compute max value for each MC simulation in point source results
        PS_max_values[MC_id] = *std::max_element(&rnd_result[options.pic_size*2], &rnd_result[options.model_size]);
    }
    std::cout << std::endl;

    // determine the value of lambda
    std::nth_element(&WS_max_values[0],
                     &WS_max_values[options.MC_quantile_PF],
                     &WS_max_values[options.MC_max]);
    options.lambda = WS_max_values[options.MC_quantile_PF];
    std::cout << "lambda = " << std::scientific << options.lambda << std::endl;

    // determine the value of lambdaI
    std::nth_element(&PS_max_values[0],
                     &PS_max_values[options.MC_quantile_PS],
                     &PS_max_values[options.MC_max]);
    double lambdaI = PS_max_values[options.MC_quantile_PS];
    std::cout << "lambdaI = " << std::scientific << lambdaI << std::endl;

    double PS_standardize_ratio = lambdaI/options.lambda;
    for(size_t i = options.pic_size*2; i < options.model_size; ++i)
        options.standardize[i] *= PS_standardize_ratio;

}

static void StandardizeAndRegularize(const Matrix<double>& background,
                                     AstroOperator<double>& astro,
                                     Parameters<double>& options)
{
    std::cout << "Computing standardization and regularization values with ";
    std::cout << options.MC_max << " MC simulations..." << std::endl;
    Matrix<double> I(0.0, options.pic_size*2, 1);
    I[0] = 1.0;
    Matrix<double> u = astro.BAW(I, false, true, true, false);
    Matrix<double> mu_hat = background + u*options.beta0;

    astro.Transpose();

    Standardize(mu_hat, background, astro, options);

    Lambda(mu_hat, astro, options);

    astro.Transpose();
    astro.Standardize(options.standardize);
    std::cout << std::endl;
}

static Matrix<double> Estimate(const Matrix<double>& image,
                               const Matrix<double>& background,
                               const AstroOperator<double>& astro,
                               Parameters<double>& options)
{
    std::cout << "Computing static estimate..." << std::endl;
    options.fista_params.iter_max = 1000;
    options.fista_params.init_value = Matrix<double>(0.0, options.model_size, 1);
    options.fista_params.init_value[0] = options.beta0 * options.standardize[0];
    options.fista_params.indices = Matrix<size_t>(0,1+options.pic_size+options.pic_size*options.pic_size,1);
    for(size_t i = 1; i < 1+options.pic_size+options.pic_size*options.pic_size; ++i)
        options.fista_params.indices[i] = i + options.pic_size - 1;

    Matrix<double> result = fista::poisson::Solve(astro, background, image, options.lambda, options.fista_params);

    result /= options.standardize;
    result.RemoveNeg(options.pic_size*2, options.model_size);

    std::cout << std::endl;

    return result;
}

static Matrix<double> EstimateNonZero(const Matrix<double>& image,
                                      const Matrix<double>& background,
                                      const Matrix<double>& solution_static,
                                      const AstroOperator<double>& astro,
                                      Parameters<double>& options)
{

    std::cout << "Getting non zero elements..." << std::endl;
    options.fista_params.iter_max = 1000;
    Matrix<size_t> non_zero_elements_indices = solution_static.NonZeroIndices();
    Matrix<double> non_zero_elements_operator(image.Length(),
                                              non_zero_elements_indices.Length());

    std::cout << "Generating non zero elements operator..." << std::endl;
    #pragma omp parallel
    {
        Matrix<double> identity(0.0, options.model_size, 1);
        #pragma omp for schedule(dynamic)
        for(size_t i = 0; i < non_zero_elements_indices.Length(); ++i)
        {
            if(i % (non_zero_elements_indices.Length()/100) == 0)
                std::cout << "\r" + std::to_string(1 + 100 * i/non_zero_elements_indices.Length()) + "/100";
            identity[non_zero_elements_indices[i]] = 1.0;
            Matrix<double> local_result = astro * identity;
            identity[non_zero_elements_indices[i]] = 0.0;

            for(size_t j = 0; j < local_result.Length(); ++j)
                non_zero_elements_operator[j*non_zero_elements_indices.Length() + i] = local_result[j];
        }
    }
    std::cout << std::endl;

    options.fista_params.init_value = Matrix<double>(non_zero_elements_indices.Length(), 1);

    #pragma omp parallel for simd
    for(size_t i = 0; i < non_zero_elements_indices.Length(); ++i)
        options.fista_params.init_value[i] = solution_static[non_zero_elements_indices[i]] * options.standardize[non_zero_elements_indices[i]];

    std::vector<size_t> remove_neg_indices = {0};
    for(size_t i = 1; i < non_zero_elements_indices.Length(); ++i)
        if(non_zero_elements_indices[i] > options.pic_size)
            remove_neg_indices.push_back(i);
    options.fista_params.indices = Matrix<size_t>(&remove_neg_indices[0], remove_neg_indices.size(), remove_neg_indices.size(), 1);

    std::cout << "Generating the new beta estimates..." << std::endl;
    Matrix<double> beta_new = fista::poisson::Solve(MatMult<double>(non_zero_elements_operator, non_zero_elements_operator.Height(), non_zero_elements_operator.Width()),
                                                    background,
                                                    image,
                                                    0.0,
                                                    options.fista_params);

    Matrix<double> result(solution_static);
    for(size_t i = 0; i < non_zero_elements_indices.Length(); ++i)
        result[non_zero_elements_indices[i]] = beta_new[i];

    return result;
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
    options.fista_params.log_period = 20;

    AstroOperator<double> astro(options.pic_size, options.pic_size, options.pic_size/2, sensitivity, Matrix<double>(1, options.model_size, 1), false, options);

    BetaZero(image, astro, options);

    options.x0[0] = options.beta0;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

    start = std::chrono::high_resolution_clock::now();
    StandardizeAndRegularize(background, astro, options);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time_MC = end-start;

    start = std::chrono::high_resolution_clock::now();
    Matrix<double> solution_static = Estimate(image, background, astro, options);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time_FISTA = end-start;

    start = std::chrono::high_resolution_clock::now();
    Matrix<double> solution = EstimateNonZero(image, background, solution_static, astro, options);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time_NZ = end-start;
    options.x0[0] = solution[0]/options.standardize[0];
    std::cout << "New beta0 = " << options.x0[0] << std::endl;

    start = std::chrono::high_resolution_clock::now();
    StandardizeAndRegularize(background, astro, options);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time_MC2 = end-start;

    for(size_t i = 0; i < options.pic_size; ++i)
        solution[i] += 1e-10;
    Matrix<size_t> zero_elements = solution.ZeroIndices();
    for(size_t i = 0; i < zero_elements.Length(); ++i)
        options.standardize[zero_elements[i]] = std::numeric_limits<double>::infinity();
    astro.Standardize(options.standardize);
    "data/512_chandra/computed_divx.data" << options.standardize;

    start = std::chrono::high_resolution_clock::now();
    solution_static = Estimate(image, background, astro, options);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time_FISTA2 = end-start;
    "data/512_chandra/computed_sol_static.data" << solution_static;

    start = std::chrono::high_resolution_clock::now();
    solution = EstimateNonZero(image, background, solution_static, astro, options);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time_NZ2 = end-start;

    Matrix<double> solution_normalized = solution / options.standardize;
    "data/512_chandra/computed_sol.data" << solution_normalized;

    options.x0[0] = solution_normalized[0];
    std::cout << "New beta0 = " << options.x0[0] << std::endl;

    std::cout << std::defaultfloat << std::endl;
    std::cout << "Time for MC simulations: " << elapsed_time_MC.count()    << " seconds" << std::endl;
    std::cout << "Time for FISTA solver: "   << elapsed_time_FISTA.count() << " seconds" << std::endl;
    std::cout << "Time for GLM fit on non zero elements: "   << elapsed_time_NZ.count() << " seconds" << std::endl;
    std::cout << "Time for MC simulations 2: " << elapsed_time_MC2.count()    << " seconds" << std::endl;
    std::cout << "Time for FISTA solver 2: "   << elapsed_time_FISTA2.count() << " seconds" << std::endl;
    std::cout << "Time for GLM fit on non zero elements 2: "   << elapsed_time_NZ2.count() << " seconds" << std::endl;
    std::cout << "Total time: "   << elapsed_time_MC.count() + elapsed_time_FISTA.count() + elapsed_time_NZ.count() + elapsed_time_MC2.count() + elapsed_time_FISTA2.count()+ elapsed_time_NZ2.count() << " seconds" << std::endl;
    std::cout << std::endl;

    return solution_normalized;
}

} // namespace WS
} // namespace astroqut
