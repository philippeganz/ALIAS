///
/// \file src/utils/linearop/operator/wavelet.cpp
/// \brief Wavelet transform class implementation
/// \details Provide various functions and tools for the Wavelet class
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2018
/// \version 0.3.0
/// \date 2018-04-08
/// \copyright GPL-3.0
///

#include "utils/linearop/operator/wavelet.hpp"

namespace astroqut
{
namespace wavelet
{

/** Make orthonormal  QMF filter
 *  \brief Generate orthonormal QMF filter for Wavelet Transform
 *  \param type Wavelet type, can be one of haar, beylkin, coiflet, daubechies, symmlet, vaidyanathan, battle
 *  \param parameter Integer parameter specific to each wavelet type
 *  \param f_type Low or high pass filter
 *  \author Jonathan Buckheit and David Donoho, MATLAB version in Wavelab 850, 1993-1995
 *  \author Philippe Ganz <philippe.ganz@gmail.com> 2018
 */
Matrix<double> MakeONFilter(wavelet_type type, int parameter, filter_type f_type)
{
    size_t data_size = 0;
    double data[59] = {0.0};

    switch(type)
    {
    case haar:
        {
            data_size = 2;
            double data_local[2] = {1/std::sqrt(2), 1/std::sqrt(2)};
            std::copy(data_local, data_local+data_size, data);
            break;
        }

    case beylkin:
        {
            data_size = 18;
            double data_local[18] = {0.099305765374, 0.424215360813, 0.699825214057, 0.449718251149, -0.110927598348, -0.264497231446, 0.026900308804, 0.155538731877, -0.017520746267, -0.088543630623, 0.019679866044, 0.042916387274, -0.017460408696, -0.014365807969, 0.010040411845, 0.001484234782, -0.002736031626, 0.000640485329};
            std::copy(data_local, data_local+data_size, data);
            break;
        }
    case coiflet:
        {
            switch(parameter)
            {
            case 1:
                {
                    data_size = 6;
                    double data_local[6] = {0.038580777748, -0.126969125396, -0.077161555496, 0.607491641386, 0.745687558934, 0.226584265197};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 2:
                {
                    data_size = 12;
                    double data_local[12] = {0.016387336463, -0.041464936782, -0.067372554722, 0.386110066823, 0.812723635450, 0.417005184424, -0.076488599078, -0.059434418646, 0.023680171947, 0.005611434819, -0.001823208871, -0.000720549445};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 3:
                {
                    data_size = 18;
                    double data_local[18] = {-0.003793512864, 0.007782596426, 0.023452696142, -0.065771911281, -0.061123390003, 0.405176902410, 0.793777222626, 0.428483476378, -0.071799821619, -0.082301927106, 0.034555027573, 0.015880544864, -0.009007976137, -0.002574517688, 0.001117518771, 0.000466216960, -0.000070983303, -0.000034599773};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 4:
                {
                    data_size = 24;
                    double data_local[24] = {0.000892313668, -0.001629492013, -0.007346166328, 0.016068943964, 0.026682300156, -0.081266699680, -0.056077313316, 0.415308407030, 0.782238930920,0.434386056491, -0.066627474263, -0.096220442034, 0.039334427123, 0.025082261845, -0.015211731527, -0.005658286686, 0.003751436157, 0.001266561929, -0.000589020757, -0.000259974552, 0.000062339034, 0.000031229876, -0.000003259680, -0.000001784985};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 5:
                {
                    data_size = 30;
                    double data_local[30] = {-0.000212080863, 0.000358589677, 0.002178236305, -0.004159358782, -0.010131117538, 0.023408156762, 0.028168029062, -0.091920010549, -0.052043163216, 0.421566206729, 0.774289603740, 0.437991626228, -0.062035963906, -0.105574208706, 0.041289208741, 0.032683574283, -0.019761779012, -0.009164231153, 0.006764185419, 0.002433373209, -0.001662863769, -0.000638131296, 0.000302259520, 0.000140541149, -0.000041340484, -0.000021315014, 0.000003734597, 0.000002063806, -0.000000167408, -0.000000095158};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            default:
                {
                    break;
                }
            }
            break;
        }
    case daubechies:
        {
            switch(parameter)
            {
            case 4:
                {
                    data_size = 4;
                    double data_local[4] = {0.482962913145, 0.836516303738, 0.224143868042, -0.129409522551};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 6:
                {
                    data_size = 6;
                    double data_local[6] = {0.332670552950, 0.806891509311, 0.459877502118, -0.135011020010, -0.085441273882, 0.035226291882};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 8:
                {
                    data_size = 8;
                    double data_local[8] = {0.230377813309, 0.714846570553, 0.630880767930, -0.027983769417, -0.187034811719, 0.030841381836, 0.032883011667, -0.010597401785};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 10:
                {
                    data_size = 10;
                    double data_local[10] = {0.160102397974, 0.603829269797, 0.724308528438, 0.138428145901, -0.242294887066, -0.032244869585, 0.077571493840, -0.006241490213, -0.012580751999, 0.003335725285};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 12:
                {
                    data_size = 12;
                    double data_local[12] = {0.111540743350, 0.494623890398, 0.751133908021, 0.315250351709, -0.226264693965, -0.129766867567, 0.097501605587, 0.027522865530, -0.031582039317, 0.000553842201, 0.004777257511, -0.001077301085};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 14:
                {
                    data_size = 14;
                    double data_local[14] = {0.077852054085, 0.396539319482, 0.729132090846, 0.469782287405, -0.143906003929, -0.224036184994, 0.071309219267, 0.080612609151, -0.038029936935, -0.016574541631, 0.012550998556, 0.000429577973, -0.001801640704, 0.000353713800};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 16:
                {
                    data_size = 16;
                    double data_local[16] = {0.054415842243, 0.312871590914, 0.675630736297, 0.585354683654, -0.015829105256, -0.284015542962, 0.000472484574, 0.128747426620, -0.017369301002, -0.044088253931, 0.013981027917, 0.008746094047, -0.004870352993, -0.000391740373, 0.000675449406, -0.000117476784};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 18:
                {
                    data_size = 18;
                    double data_local[18] = {0.038077947364, 0.243834674613, 0.604823123690, 0.657288078051, 0.133197385825, -0.293273783279, -0.096840783223, 0.148540749338, 0.030725681479, -0.067632829061, 0.000250947115, 0.022361662124, -0.004723204758, -0.004281503682, 0.001847646883, 0.000230385764, -0.000251963189, 0.000039347320};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 20:
                {
                    data_size = 20;
                    double data_local[20] = {0.026670057901, 0.188176800078, 0.527201188932, 0.688459039454, 0.281172343661, -0.249846424327, -0.195946274377, 0.127369340336, 0.093057364604, -0.071394147166, -0.029457536822, 0.033212674059, 0.003606553567, -0.010733175483, 0.001395351747, 0.001992405295, -0.000685856695, -0.000116466855, 0.000093588670, -0.000013264203};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            default:
                {
                    break;
                }
            }
            break;
        }
    case symmlet:
        {
            switch(parameter)
            {
            case 4:
                {
                    data_size = 8;
                    double data_local[8] = {-0.107148901418, -0.041910965125, 0.703739068656, 1.136658243408, 0.421234534204, -0.140317624179, -0.017824701442, 0.045570345896};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 5:
                {
                    data_size = 10;
                    double data_local[10] = {0.038654795955, 0.041746864422, -0.055344186117, 0.281990696854, 1.023052966894, 0.896581648380, 0.023478923136, -0.247951362613, -0.029842499869, 0.027632152958};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 6:
                {
                    data_size = 12;
                    double data_local[12] = {0.021784700327, 0.004936612372, -0.166863215412, -0.068323121587, 0.694457972958, 1.113892783926, 0.477904371333, -0.102724969862, -0.029783751299, 0.063250562660, 0.002499922093, -0.011031867509};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 7:
                {
                    data_size = 14;
                    double data_local[14] = {0.003792658534, -0.001481225915, -0.017870431651, 0.043155452582, 0.096014767936, -0.070078291222, 0.024665659489, 0.758162601964, 1.085782709814, 0.408183939725, -0.198056706807, -0.152463871896, 0.005671342686, 0.014521394762};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 8:
                {
                    data_size = 16;
                    double data_local[16] = {0.002672793393, -0.000428394300, -0.021145686528, 0.005386388754, 0.069490465911, -0.038493521263, -0.073462508761, 0.515398670374, 1.099106630537, 0.680745347190, -0.086653615406, -0.202648655286, 0.010758611751, 0.044823623042, -0.000766690896, -0.004783458512};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 9:
                {
                    data_size = 18;
                    double data_local[18] = {0.001512487309, -0.000669141509, -0.014515578553, 0.012528896242, 0.087791251554, -0.025786445930, -0.270893783503, 0.049882830959, 0.873048407349, 1.015259790832, 0.337658923602, -0.077172161097, 0.000825140929, 0.042744433602, -0.016303351226, -0.018769396836, 0.000876502539, 0.001981193736};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 10:
                {
                    data_size = 20;
                    double data_local[20] = {0.001089170447, 0.000135245020, -0.012220642630, -0.002072363923, 0.064950924579, 0.016418869426, -0.225558972234, -0.100240215031, 0.667071338154, 1.088251530500, 0.542813011213, -0.050256540092, -0.045240772218, 0.070703567550, 0.008152816799, -0.028786231926, -0.001137535314, 0.006495728375, 0.000080661204, -0.000649589896};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            default:
                {
                    break;
                }
            }
            break;
        }
    case vaidyanathan:
        {
            data_size = 24;
            double data_local[24] = {-0.000062906118, 0.000343631905, -0.000453956620, -0.000944897136, 0.002843834547, 0.000708137504, -0.008839103409, 0.003153847056, 0.019687215010, -0.014853448005, -0.035470398607, 0.038742619293, 0.055892523691, -0.077709750902, -0.083928884366, 0.131971661417, 0.135084227129, -0.194450471766, -0.263494802488, 0.201612161775, 0.635601059872, 0.572797793211, 0.250184129505, 0.045799334111};
            std::copy(data_local, data_local+data_size, data);
            break;
        }
    case battle:
        {
            switch(parameter)
            {
            case 1:
                {
                    data_size = 23;
                    double data_local[23] = {-0.000086752300000, -0.000158601000000, 0.000361781000000, 0.000652922000000, -0.001557010000000, -0.002745880000000, 0.007064420000000, 0.012003000000000, -0.036730900000000, -0.048861800000000, 0.280931000000000, 0.578163000000000, 0.280931000000000, -0.048861800000000, -0.036730900000000, 0.012003000000000, 0.007064420000000, -0.002745880000000, -0.001557010000000, 0.000652922000000, 0.000361781000000, -0.000158601000000, -0.000086752300000};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 3:
                {
                    data_size = 41;
                    double data_local[41] = {0.000103307000000, -0.000164264000000, -0.000201818000000, 0.000326749000000, 0.000395946000000, -0.000655620000000, -0.000780468000000, 0.001330860000000, 0.001546240000000, -0.002745290000000, -0.003078630000000, 0.005799320000000, 0.006141430000000, -0.012715400000000, -0.012145500000000, 0.029746800000000, 0.022684600000000, -0.077807900000000, -0.035498000000000, 0.306830000000000, 0.541736000000000, 0.306830000000000, -0.035498000000000, -0.077807900000000, 0.022684600000000, 0.029746800000000, -0.012145500000000, -0.012715400000000, 0.006141430000000, 0.005799320000000, -0.003078630000000, -0.002745290000000, 0.001546240000000, 0.001330860000000, -0.000780468000000, -0.000655620000000, 0.000395946000000, 0.000326749000000, -0.000201818000000, -0.000164264000000, 0.000103307000000};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            case 5:
                {
                    data_size = 59;
                    double data_local[59] = {0.000101113000000, 0.000110709000000, -0.000159168000000, -0.000172685000000, 0.000251419000000, 0.000269842000000, -0.000398759000000, -0.000422485000000, 0.000635563000000, 0.000662836000000, -0.001019120000000, -0.001042070000000, 0.001646590000000, 0.001641320000000, -0.002686460000000, -0.002588160000000, 0.004440020000000, 0.004078820000000, -0.007468480000000, -0.006398860000000, 0.012875400000000, 0.009906350000000, -0.022995100000000, -0.014853700000000, 0.043354400000000, 0.020841400000000, -0.091406800000000, -0.026177100000000, 0.312869000000000, 0.528374000000000, 0.312869000000000, -0.026177100000000, -0.091406800000000, 0.020841400000000, 0.043354400000000, -0.014853700000000, -0.022995100000000, 0.009906350000000, 0.012875400000000, -0.006398860000000, -0.007468480000000, 0.004078820000000, 0.004440020000000, -0.002588160000000, -0.002686460000000, 0.001641320000000, 0.001646590000000, -0.001042070000000, -0.001019120000000, 0.000662836000000, 0.000635563000000, -0.000422485000000, -0.000398759000000, 0.000269842000000, 0.000251419000000, -0.000172685000000, -0.000159168000000, 0.000110709000000, 0.000101113000000};
                    std::copy(data_local, data_local+data_size, data);
                    break;
                }
            default:
                {
                    break;
                }
            }
            break;
        }
    default:
        {
            break;
        }
    }

    if( f_type == high )
    {
        for( size_t i = 1; i < data_size; i += 2 )
        {
            data[i] = -data[i];
        }
    }

    Matrix<double> result(data, data_size, data_size, 1);
    return result/result.Norm(two);
}

/** Forward Wavelet Transform (periodized, orthogonal)
 *  \brief Applies a periodized and orthogonal discrete wavelet transform.
 *  \param signal Signal to transform, must be length a power of 2
 *  \param wcoef Result array, must be the same size as signal
 *  \param column Column to transform
 *  \param coarsest_level Coarsest level of the wavelet transform
 *  \param low_pass_filter Quadrature mirror filter for low pass filtering
 *  \param high_pass_filter Mirrored quadrature mirror filter for high pass filtering
 *  \param intermediate Temporary array of size 1 x Height of signal
 *  \param intermediate_temp Temporary array of size 1 x Height of signal
 *  \author David Donoho <donoho@stat.stanford.edu> 1993
 *  \author Philippe Ganz <philippe.ganz@gmail.com> 2018
 */
void FWT_PO(const Matrix<double>& signal,
            Matrix<double>& wcoef,
            unsigned int column,
            unsigned int coarsest_level,
            const Matrix<double>& low_pass_filter,
            const Matrix<double>& high_pass_filter,
            double* intermediate,
            double* intermediate_temp )
{
    size_t level_max = (size_t) std::ceil(std::log2(signal.Height()));
    size_t level_offset = signal.Height();

#ifdef DO_ARGCHECKS
    if( (size_t) std::pow(2,level_max) != signal.Length() )
    {
        std::cerr << "Signal height must be length a power of two." << std::endl;
        throw;
    }

    if( coarsest_level >= level_max )
    {
        std::cerr << "The coarsest level must be in the [0, " << level_max << ") range." << std::endl;
        throw;
    }

    if( column >= signal.Width() )
    {
        std::cerr << "The column must be in the [0, " << signal.Width() << ") range." << std::endl;
        throw;
    }
#endif // DO_ARGCHECKS

    for( size_t i = 0; i < signal.Height(); ++i )
    {
        intermediate[i] = signal[i*wcoef.Width() + column];
    }

    for( size_t level = level_max; level > coarsest_level; --level )
    {
        for( size_t pass_index = 0; pass_index < level_offset/2; ++pass_index )
        {
            double low_pass_local_coef = 0.0;
            size_t low_pass_offset = 2*pass_index;
            double high_pass_local_coef = 0.0;
            int high_pass_offset = 2*pass_index+1;

            for( size_t filter_index = 0; filter_index < low_pass_filter.Length(); ++filter_index )
            {
                low_pass_local_coef += low_pass_filter[filter_index] * intermediate[low_pass_offset];

                ++low_pass_offset;
                if( low_pass_offset >= level_offset )
                {
                    low_pass_offset -= level_offset;
                }

                high_pass_local_coef += high_pass_filter[filter_index] * intermediate[high_pass_offset];

                --high_pass_offset;
                if( high_pass_offset < 0 )
                {
                    high_pass_offset += level_offset;
                }
            }

            intermediate_temp[pass_index] = low_pass_local_coef;
            intermediate_temp[pass_index + level_offset/2] = high_pass_local_coef;
        }

        for( size_t i = 0; i < level_offset; ++i )
        {
            intermediate[i] = intermediate_temp[i];
        }

        level_offset /= 2;
    }

    for( size_t i = 0; i < signal.Height(); ++i )
    {
        wcoef[i*wcoef.Width() + column] = intermediate[i];
    }
}

/** Inverse Wavelet Transform (periodized, orthogonal)
 *  \brief Applies a periodized and orthogonal inverse discrete wavelet transform.
 *  \param wcoef Wavelet coefficients to transform back, must be length a power of 2
 *  \param signal Result array, must be the same size as wcoef
 *  \param column Column to transform, -1 to transform all
 *  \param coarsest_level Coarsest level of the wavelet transform
 *  \param low_pass_filter Quadrature mirror filter for low pass filtering
 *  \param high_pass_filter Mirrored quadrature mirror filter for high pass filtering
 *  \param intermediate Temporary array of size 1 x Height of signal
 *  \param intermediate_temp Temporary array of size 1 x Height of signal
 *  \author David Donoho <donoho@stat.stanford.edu> 1993
 *  \author Philippe Ganz <philippe.ganz@gmail.com> 2018
 */
void IWT_PO(const Matrix<double>& wcoef,
            Matrix<double>& signal,
            unsigned int column,
            unsigned int coarsest_level,
            const Matrix<double>& low_pass_filter,
            const Matrix<double>& high_pass_filter,
            double* intermediate,
            double* intermediate_temp )
{
    size_t level_max = (size_t) std::ceil(std::log2(signal.Height()));
    size_t level_offset = 1;
    size_t filter_length = low_pass_filter.Length();
    size_t filter_length_half_even = (filter_length + 1) / 2;
    size_t filter_length_half_odd = filter_length / 2;

#ifdef DO_ARGCHECKS
    if( (size_t) std::pow(2,level_max) != signal.Length() )
    {
        std::cerr << "Signal height must be length a power of two." << std::endl;
        throw;
    }

    if( coarsest_level >= level_max )
    {
        std::cerr << "The coarsest level must be in the [0, " << level_max << ") range." << std::endl;
        throw;
    }

    if( column >= signal.Width() )
    {
        std::cerr << "The column must be in the [0, " << signal.Width() << ") range." << std::endl;
        throw;
    }
#endif // DO_ARGCHECKS

    for( size_t i = 0; i < (size_t) std::pow(2, coarsest_level); ++i )
    {
        intermediate[i] = wcoef[i*wcoef.Width() + column];
    }

    for( size_t level = (size_t) std::pow(2, coarsest_level); level <= level_max; ++level )
    {
        for( size_t pass_index = 0; pass_index < level_offset; ++pass_index )
        {
            double even_local_coef = 0.0;
            int low_pass_offset = pass_index;
            double odd_local_coef = 0.0;
            size_t high_pass_offset = pass_index;

            for( size_t filter_index = 0; filter_index < filter_length_half_even; ++filter_index )
            {
                even_local_coef += low_pass_filter[2*filter_index] * intermediate[low_pass_offset];

                --low_pass_offset;
                if( low_pass_offset < 0 )
                {
                    low_pass_offset += level_offset;
                }

                odd_local_coef += high_pass_filter[2*filter_index] * wcoef[(level_offset + high_pass_offset)*wcoef.Width() + column];

                ++high_pass_offset;
                if( high_pass_offset >= level_offset )
                {
                    high_pass_offset -= level_offset;
                }
            }

            low_pass_offset = pass_index;
            high_pass_offset = pass_index;
            for( size_t filter_index = 0; filter_index < filter_length_half_odd; ++filter_index )
            {
                odd_local_coef += low_pass_filter[2*filter_index+1] * intermediate[low_pass_offset];

                --low_pass_offset;
                if( low_pass_offset < 0 )
                {
                    low_pass_offset += level_offset;
                }

                even_local_coef += high_pass_filter[2*filter_index+1] * wcoef[(level_offset + high_pass_offset)*wcoef.Width() + column];

                ++high_pass_offset;
                if( high_pass_offset >= level_offset )
                {
                    high_pass_offset -= level_offset;
                }
            }

            intermediate_temp[2*pass_index] = even_local_coef;
            intermediate_temp[2*pass_index + 1] = odd_local_coef;
        }

        for( size_t i = 0; i < 2*level_offset; ++i )
        {
            intermediate[i] = intermediate_temp[i];
        }

        level_offset *= 2;
    }

    for( size_t i = 0; i < signal.Height(); ++i )
    {
        signal[i*signal.Width() + column] = intermediate[i];
    }

}



} // namespace wavelet
} // namespace astroqut
