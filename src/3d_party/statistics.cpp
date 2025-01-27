//----------------------------------------------------------------------------------------------------------------------
///
/// \file       statistics.cpp
/// \brief      Различные функции математической статистики
/// \date       08.02.21 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#include <statistics.h>

namespace SPML /// Специальная библиотека программных модулей (СБПМ)
{
namespace Statistics
{
const int FactNumMaxInt32 = 12; ///< Максимальное целое число, факториал которого укладывается в 4 байта
//----------------------------------------------------------------------------------------------------------------------
unsigned int Factorial( unsigned int n )
{
    assert( n <= FactNumMaxInt32 );
    if( n == 0 ) {
        return 1;
    }
    return ( n * Factorial( n - 1 ) );
}
//----------------------------------------------------------------------------------------------------------------------
unsigned int C_k_n( unsigned int k, unsigned int n )
{
    if( k > n ) {
        return 0;
    } else if( k == 0 || k == n ) {
        return 1;
    } else if( k == 1 || k == n - 1 ) {
        return n;
    } else {
        return C_k_n( k - 1, n - 1 ) * n / k; // recursive
    }
}

//----------------------------------------------------------------------------------------------------------------------
bool IsEven( int value )
{
    if( value % 2 == 0 ) {
        return true; // чётное
    }
    return false; // нечётное
}
//----------------------------------------------------------------------------------------------------------------------
double VolumeOfSphere_formula( unsigned int n, double r )
{    
    assert( r >= 0 );

    if( n == 0 ) {
        return 1.0;
    }
    if( n == 1 ) {
        return ( 2.0 * r );
    }
    int k = 0; //    int k = n;
    double result = 0;

    if( IsEven( n ) ) { // чётное n
        k = static_cast<int>( std::floor( n / 2 ) );
        result = ( ( std::pow( Consts::PI, k ) / Factorial( k ) ) * std::pow( r, 2 * k ) );
    } else {
        k = static_cast<int>( std::floor( ( n - 1 ) / 2 ) );
        result = ( ( ( 2 * Factorial( k ) * std::pow( ( 4 * Consts::PI ), k ) ) / Factorial( ( 2 * k ) + 1 ) )
            * std::pow( r, ( 2 * k ) + 1 ) );
    }
    return result;
}

double VolumeOfSphere_hardcoded( unsigned int n, double r )
{    
    assert( r >= 0 );

    if( n == 0 ) {
        return 1.0;
    }
    if( n == 1 ) {
        return ( 2.0 * r );
    }
    double result = 0;
    switch( n ) { // в зависимости об размерности пространства измерений
        case 2:
        {
            result = Consts::PI;
            break;
        }
        case 3:
        {
            result = 4.0 * Consts::PI / 3.0;
            break;
        }
        case 4:
        {
            result = Consts::PI * Consts::PI / 2.0;
            break;
        }
        case 5:
        {
            result = 8.0 * Consts::PI * Consts::PI / 15.0;
            break;
        }
        case 6:
        {
            result = Consts::PI * Consts::PI * Consts::PI / 6.0;
            break;
        }
        case 7:
        {
            result = 16.0 * Consts::PI * Consts::PI * Consts::PI / 105.0;
            break;
        }
        case 8:
        {
            result = Consts::PI * Consts::PI * Consts::PI * Consts::PI / 24.0;
            break;
        }
        case 9:
        {
            result = 32.0 * Consts::PI * Consts::PI * Consts::PI * Consts::PI / 945.0;
            break;
        }
        case 10:
        {
            result = Consts::PI * Consts::PI * Consts::PI * Consts::PI * Consts::PI / 120.0;
            break;
        }
        default:
            assert( false );
    }
    return ( result * std::pow( r, n ) );
}
//----------------------------------------------------------------------------------------------------------------------
double CalcRVSTR( double det, double strobe, double volume, int sizeY )
{
    assert( !Compare::IsZeroAbs( volume ) ); // Объем зоны не должен быть равен нулю
    assert( det > 0.0 );
    double Vcoef = VolumeOfSphere_formula( sizeY ); // Коэффициент объема эллипсоида рассеяния
//    double result = ( ( Vcoef * ( std::sqrt( std::abs( det ) ) ) * std::pow( strobe, sizeY ) ) / volume );
    double result = ( ( Vcoef * ( std::sqrt( det ) ) * std::pow( strobe, sizeY ) ) / volume );
    return result;
}
//----------------------------------------------------------------------------------------------------------------------
double ErfInv( double x )
{
    if( ( x < -1.0 ) || ( x > 1.0 ) ) {
        return std::numeric_limits<double>::quiet_NaN();//        return NAN;
    } else if( x == 1.0 ) {
        return std::numeric_limits<double>::infinity();//        return INFINITY;
    } else if( x == -1.0 ) {
        return -std::numeric_limits<double>::infinity();//        return -INFINITY;
    }

    const double LN2 = 6.931471805599453094172321214581e-1L;

    const double A0_ = 1.1975323115670912564578e0L;
    const double A1_ = 4.7072688112383978012285e1L;
    const double A2_ = 6.9706266534389598238465e2L;
    const double A3_ = 4.8548868893843886794648e3L;
    const double A4_ = 1.6235862515167575384252e4L;
    const double A5_ = 2.3782041382114385731252e4L;
    const double A6_ = 1.1819493347062294404278e4L;
    const double A7_ = 8.8709406962545514830200e2L;

    const double B0_ = 1.0000000000000000000e0L;
    const double B1_ = 4.2313330701600911252e1L;
    const double B2_ = 6.8718700749205790830e2L;
    const double B3_ = 5.3941960214247511077e3L;
    const double B4_ = 2.1213794301586595867e4L;
    const double B5_ = 3.9307895800092710610e4L;
    const double B6_ = 2.8729085735721942674e4L;
    const double B7_ = 5.2264952788528545610e3L;

    const double C0_ = 1.42343711074968357734e0L;
    const double C1_ = 4.63033784615654529590e0L;
    const double C2_ = 5.76949722146069140550e0L;
    const double C3_ = 3.64784832476320460504e0L;
    const double C4_ = 1.27045825245236838258e0L;
    const double C5_ = 2.41780725177450611770e-1L;
    const double C6_ = 2.27238449892691845833e-2L;
    const double C7_ = 7.74545014278341407640e-4L;

    const double D0_ = 1.4142135623730950488016887e0L;
    const double D1_ = 2.9036514445419946173133295e0L;
    const double D2_ = 2.3707661626024532365971225e0L;
    const double D3_ = 9.7547832001787427186894837e-1L;
    const double D4_ = 2.0945065210512749128288442e-1L;
    const double D5_ = 2.1494160384252876777097297e-2L;
    const double D6_ = 7.7441459065157709165577218e-4L;
    const double D7_ = 1.4859850019840355905497876e-9L;

    const double E0_ = 6.65790464350110377720e0L;
    const double E1_ = 5.46378491116411436990e0L;
    const double E2_ = 1.78482653991729133580e0L;
    const double E3_ = 2.96560571828504891230e-1L;
    const double E4_ = 2.65321895265761230930e-2L;
    const double E5_ = 1.24266094738807843860e-3L;
    const double E6_ = 2.71155556874348757815e-5L;
    const double E7_ = 2.01033439929228813265e-7L;

    const double F0_ = 1.414213562373095048801689e0L;
    const double F1_ = 8.482908416595164588112026e-1L;
    const double F2_ = 1.936480946950659106176712e-1L;
    const double F3_ = 2.103693768272068968719679e-2L;
    const double F4_ = 1.112800997078859844711555e-3L;
    const double F5_ = 2.611088405080593625138020e-5L;
    const double F6_ = 2.010321207683943062279931e-7L;
    const double F7_ = 2.891024605872965461538222e-15L;

    double abs_x = std::abs(x);

    if( abs_x <= 0.85L ) {
        double r = 0.180625L - 0.25L * x * x;
        double num = (((((((A7_ * r + A6_) * r + A5_) * r + A4_) * r + A3_) * r + A2_) * r + A1_) * r + A0_);
        double den = (((((((B7_ * r + B6_) * r + B5_) * r + B4_) * r + B3_) * r + B2_) * r + B1_) * r + B0_);
        return ( x * num / den );
    }

    double r = std::sqrt( LN2 - logl( 1.0L - abs_x ) );

    double num, den;
    if( r <= 5.0L ) {
        r = r - 1.6L;
        num = (((((((C7_ * r + C6_) * r + C5_) * r + C4_) * r + C3_) * r + C2_) * r + C1_) * r + C0_);
        den = (((((((D7_ * r + D6_) * r + D5_) * r + D4_) * r + D3_) * r + D2_) * r + D1_) * r + D0_);
    } else {
        r = r - 5.0L;
        num = (((((((E7_ * r + E6_) * r + E5_) * r + E4_) * r + E3_) * r + E2_) * r + E1_) * r + E0_);
        den = (((((((F7_ * r + F6_) * r + F5_) * r + F4_) * r + F3_) * r + F2_) * r + F1_) * r + F0_);
    }
    return std::copysignl( num / den, x );

//    if( ( x < -1.0 ) || ( x > 1.0 ) ) {
//        return std::numeric_limits<double>::quiet_NaN();//        return NAN;
//    } else if( x == 1.0 ) {
//        return std::numeric_limits<double>::infinity();//        return INFINITY;
//    } else if( x == -1.0 ) {
//        return -std::numeric_limits<double>::infinity();//        return -INFINITY;
//    }

//    const long double LN2 = 6.931471805599453094172321214581e-1L;

//    const long double A0 = 1.1975323115670912564578e0L;
//    const long double A1 = 4.7072688112383978012285e1L;
//    const long double A2 = 6.9706266534389598238465e2L;
//    const long double A3 = 4.8548868893843886794648e3L;
//    const long double A4 = 1.6235862515167575384252e4L;
//    const long double A5 = 2.3782041382114385731252e4L;
//    const long double A6 = 1.1819493347062294404278e4L;
//    const long double A7 = 8.8709406962545514830200e2L;

//    const long double B0 = 1.0000000000000000000e0L;
//    const long double B1 = 4.2313330701600911252e1L;
//    const long double B2 = 6.8718700749205790830e2L;
//    const long double B3 = 5.3941960214247511077e3L;
//    const long double B4 = 2.1213794301586595867e4L;
//    const long double B5 = 3.9307895800092710610e4L;
//    const long double B6 = 2.8729085735721942674e4L;
//    const long double B7 = 5.2264952788528545610e3L;

//    const long double C0 = 1.42343711074968357734e0L;
//    const long double C1 = 4.63033784615654529590e0L;
//    const long double C2 = 5.76949722146069140550e0L;
//    const long double C3 = 3.64784832476320460504e0L;
//    const long double C4 = 1.27045825245236838258e0L;
//    const long double C5 = 2.41780725177450611770e-1L;
//    const long double C6 = 2.27238449892691845833e-2L;
//    const long double C7 = 7.74545014278341407640e-4L;

//    const long double D0 = 1.4142135623730950488016887e0L;
//    const long double D1 = 2.9036514445419946173133295e0L;
//    const long double D2 = 2.3707661626024532365971225e0L;
//    const long double D3 = 9.7547832001787427186894837e-1L;
//    const long double D4 = 2.0945065210512749128288442e-1L;
//    const long double D5 = 2.1494160384252876777097297e-2L;
//    const long double D6 = 7.7441459065157709165577218e-4L;
//    const long double D7 = 1.4859850019840355905497876e-9L;

//    const long double E0 = 6.65790464350110377720e0L;
//    const long double E1 = 5.46378491116411436990e0L;
//    const long double E2 = 1.78482653991729133580e0L;
//    const long double E3 = 2.96560571828504891230e-1L;
//    const long double E4 = 2.65321895265761230930e-2L;
//    const long double E5 = 1.24266094738807843860e-3L;
//    const long double E6 = 2.71155556874348757815e-5L;
//    const long double E7 = 2.01033439929228813265e-7L;

//    const long double F0 = 1.414213562373095048801689e0L;
//    const long double F1 = 8.482908416595164588112026e-1L;
//    const long double F2 = 1.936480946950659106176712e-1L;
//    const long double F3 = 2.103693768272068968719679e-2L;
//    const long double F4 = 1.112800997078859844711555e-3L;
//    const long double F5 = 2.611088405080593625138020e-5L;
//    const long double F6 = 2.010321207683943062279931e-7L;
//    const long double F7 = 2.891024605872965461538222e-15L;

//    long double abs_x = std::abs(x);

//    if( abs_x <= 0.85L ) {
//        long double r = 0.180625L - 0.25L * x * x;
//        long double num = (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0);
//        long double den = (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + B0);
//        return ( x * num / den );
//    }

//    long double r = std::sqrt( LN2 - logl( 1.0L - abs_x ) );

//    long double num, den;
//    if( r <= 5.0L ) {
//        r = r - 1.6L;
//        num = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0);
//        den = (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + D0);
//    } else {
//        r = r - 5.0L;
//        num = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0);
//        den = (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + F0);
//    }
//    return std::copysignl( num / den, x );
}
//----------------------------------------------------------------------------------------------------------------------
double ChiSquaredCDF_PeizerPratt( double arg, unsigned int degree )
{
    assert( degree > 1 ); // Ограничение на число степеней свободы для Пейзера-Пратта
    double t = arg * arg;
    double x = 0.0;
    double degree_ = static_cast<double>( degree );
    if( Compare::AreEqualAbs( t, degree_ - 1.0 ) ) { // if( t == v-1 )
        x = ( ( ( -1.0 / 3.0 ) + ( 0.08 / degree_ ) ) / std::sqrt( ( 2.0 * degree_ ) - 2.0 ) );
    } else {
        x = ( ( t - degree_ + ( 2.0 / 3.0 ) - ( 0.08 / degree_ ) ) / std::abs( t - ( degree_ - 1.0 ) ) ) *
            std::sqrt( ( degree_ - 1.0 ) * std::log( ( degree_ - 1.0 ) / t ) + t - ( degree_ - 1.0 ) );
    }
    double result = GaussCDF( x, 0, 1 );
    return result;
}
//----------------------------------------------------------------------------------------------------------------------
double ChiSquaredCDF_EvenDegree( double arg, unsigned int degree )
{
    assert( IsEven( degree ) );
    assert( degree > 0 && degree < FactNumMaxInt32 ); // Ограничение на число степеней свободы
    double sum = 0;
    double arg_ = arg * arg;
    int degree_ = degree / 2;
    for( int k = 0; k < degree_; k++ ) {
        sum += ( ( 1.0 / Factorial( k ) ) * std::pow( ( arg_ * 0.5 ), k ) );
    }
    double result = ( 1.0 - ( std::exp( -arg_ * 0.5 ) * sum ) );
    return result;
}
//----------------------------------------------------------------------------------------------------------------------
/*
double Freqs( double arg, int degree )
{
    assert( degree > 2 && degree < 10 );        // Ограничение на число степеней свободы
    double F = ChiSquaredCDF_PeizerPratt( arg, degree ); // Используем аппрокисимацию
//    boost::math::chi_squared chi_squared_degree( degree );
//    double F = chi_squared_degree( arg );
    F = 1.0 - F;
    // Ограничим сверху и снизу
    if( F < 0.00001 ) {
        F = 0.00001;
    }
    if( F > 0.99999 ) {
        F = 0.99999;
    }
    return F;
}

double FreqsNew( double arg, int degree )
{
    assert( degree > 2 && degree < 10 );        // Ограничение на число степеней свободы
//    double F = ChiSqPeizerPratt( arg, degree ); // Используем аппрокисимацию
    boost::math::chi_squared chi_squared_degree( degree );
    double F = boost::math::cdf( chi_squared_degree, arg * arg );
    F = 1.0 - F;
    // Ограничим сверху и снизу
    if( F < 0.00001 ) {
        F = 0.00001;
    }
    if( F > 0.99999 ) {
        F = 0.99999;
    }
    return F;
}

double FreqsNew2( double arg, int degree )
{
    assert( degree > 2 && degree < 10 );        // Ограничение на число степеней свободы

    double F = ChiSquaredCDF_EvenDegree( arg, degree );
    F = 1.0 - F;
    // Ограничим сверху и снизу
    if( F < 0.00001 ) {
        F = 0.00001;
    }
    if( F > 0.99999 ) {
        F = 0.99999;
    }
    return F;
}
*/
//----------------------------------------------------------------------------------------------------------------------
double ChiSquaredCDF( double arg, unsigned int degree, const TChiSquaredMethod method )
{
    double result = 0;
    switch( method ) {
        case( TChiSquaredMethod::CSM_Boost ):
        {
            boost::math::chi_squared chi_squared_degree( degree );
            result = boost::math::cdf( chi_squared_degree, arg * arg );
            break;
        }
        case( TChiSquaredMethod::CSM_EvenDegree ):
        {
            if( IsEven( degree ) ) {
                result = ChiSquaredCDF_EvenDegree( arg, degree );
                break;
            }
            // Иначе провалиться ниже в расчет по аппроксимации
        }
        case( TChiSquaredMethod::CSM_Approx ):
        {
            result = ChiSquaredCDF_PeizerPratt( arg, degree );
            break;
        }
        default:
            assert( false );
    }
    // Ограничим результат сверху и снизу
//    double min = 1.0e-5;
    double min = 1.0e-7;
    double max = 1.0 - min;
    if( result < min ) {
        result = min;
    }
    if( result > max ) {
        result = max;
    }
    return result;
}

double InvChiSquaredCDF( double arg, unsigned int degree, const TChiSquaredMethod method )
{
    double result = 1.0 - ChiSquaredCDF( arg, degree, method );
    return result;
}

//----------------------------------------------------------------------------------------------------------------------
bool WaldCriterion( double pl, double &pc, double d, double &dcorr, double f, double ktsol_max, double &ktsol_corr,
    double &ad, double &ar, double &b )
{
    // Начальные значения порогов
    ad = -1.0e6;
    ar = 1.0e6;
    b = 1.0;

    int kt = static_cast<int>( ktsol_max ); // Скорректированное число тактов на принятие решения - его надо выдать!

    // Ограничим pl
    double tmp = 1e-6;
    double pli = 0.0;
    pli = std::max( pl, tmp );
    pli = std::min( pl, 1.0 - tmp );

    double tmp0, tmp1, tmp2, tmp3;
    int ntakts;
    dcorr = d; // Начальное значение dcorr

    const double pc_start = pc;
    const double pc_end = 1.0;
    double pc_step = 0.01;
    do {        
        for( pc = pc_start; pc <= pc_end; pc = pc + pc_step ) {
            if( pc <= pli ) continue;
            for( int j = 1; j < kt; j++ ) {
                dcorr = d / ( 1.0 - std::pow( ( 1.0 - pc ), j ) );
                if( dcorr > 1.0 ) continue;

                tmp0 = std::log( f / dcorr );
                tmp1 = std::log( ( 1.0 - f ) / ( 1.0 - dcorr ) );
                tmp2 = std::log( pli / pc );
                tmp3 = std::log( ( 1.0 - pli ) / ( 1.0 - pc ) );

                ad = - tmp0 / ( tmp2 - tmp3 );
                ar = - tmp1 / ( tmp2 - tmp3 );
                b = tmp2 / ( tmp2 - tmp3 );

                double ntakts_double = std::ceil( ( tmp0 * tmp1 ) / ( tmp2 * tmp3 ) ); // Округление вверх всегда
                ntakts = static_cast<int>( ntakts_double );

                if( ( kt - ntakts ) >= j ) {
                    kt = kt - j;
                    goto LoopEnd; // нужно выйти из двух вложенных циклов - goto оправдан
                }
            }
        }
LoopEnd:;
        if( ( std::abs( ad ) + std::abs( ar ) ) > 1.0e6 ) {
            dcorr *= 0.9; // correct d
        } else {
            break; // break do-while loop - we 've finished
        }
    } while( dcorr >= 0.7 ); // 0.7 - нижний порог понижения d
    ktsol_corr = static_cast<double>( kt );

    if( dcorr <= 0.7 ) {
        return false; // Fail!
    }
    return true; // OK
}

int KfromNcriterion( double f, int n, double pl )
{
    int k;
    for( k = 1; k <= n; k++ ) {
        double f_calculated = SPML::Statistics::C_k_n( k, n ) * std::pow( pl, k ) * std::pow( ( 1.0 - pl ), ( n - k ) );
        if( f_calculated < f ) {
            return k;
        }
    }
    return n;
}

}
}
/// \}
