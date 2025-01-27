//----------------------------------------------------------------------------------------------------------------------
///
/// \file       geodesy.cpp
/// \brief      Земные эллипсоиды, геодезические задачи (персчеты координат)
/// \date       06.11.19 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#include <geodesy.h>

namespace SPML /// Специальная библиотека программных модулей (СБПМ)
{
namespace Geodesy /// Геодезические функции и функции перевода координат
{

//----------------------------------------------------------------------------------------------------------------------
std::string CEllipsoid::Name() const { return name; }

double CEllipsoid::A() const { return a; }

double CEllipsoid::B() const { return b; }

double CEllipsoid::F() const { return f; }

double CEllipsoid::Invf() const { return invf; }

double CEllipsoid::EccentricityFirst() const
{
    return ( std::sqrt( ( a * a ) - ( b * b ) ) / a );
}

double CEllipsoid::EccentricityFirstSquared() const
{
    return ( 1.0 - ( ( b * b ) / ( a * a ) ) );
}

double CEllipsoid::EccentricitySecond() const
{
    return ( std::sqrt( ( a * a ) - ( b * b ) ) / b );
}

double CEllipsoid::EccentricitySecondSquared() const
{
    return ( ( ( a * a ) / ( b * b ) ) - 1.0 );
}

CEllipsoid::CEllipsoid()
{
    a = 0.0;
    b = 0.0;
    f = 0.0;
    invf = 0.0;
}

CEllipsoid::CEllipsoid( std::string ellipsoidName, double semiMajorAxis, double semiMinorAxis, double inverseFlattening, bool isInvfDef )
{
    name = ellipsoidName;
    a = semiMajorAxis;
    invf = inverseFlattening;
    if( isInvfDef && ( Compare::IsZeroAbs( inverseFlattening ) || std::isinf( inverseFlattening ) ) ) {
        b = semiMajorAxis;
        f = 0.0;
    } else if ( isInvfDef ) {
        b = ( 1.0 - ( 1.0 / inverseFlattening ) ) * semiMajorAxis;
        f = 1.0 / inverseFlattening;
    } else {
        b = semiMinorAxis;
        f = 1.0 / inverseFlattening;
    }
}

//CEllipsoid::CEllipsoid( std::string ellipsoidName, double semiMajorAxis, double semiMinorAxis, double inverseFlattening )
//{
//    name = ellipsoidName;
//    a = semiMajorAxis;
//    invf = inverseFlattening;
//    if( Compare::IsZero( semiMinorAxis ) && !Compare::IsZero( inverseFlattening ) ) { // b = 0, invf != 0
//        b = ( 1.0 - ( 1.0 / inverseFlattening ) ) * semiMajorAxis;
//        f = 1.0 / inverseFlattening;
//    } else if( !Compare::IsZero( semiMinorAxis ) && Compare::IsZero( inverseFlattening ) ) {
//        b = semiMinorAxis;
//        f = 0.0;
//    } else {
//        b = semiMinorAxis;
//        f = 1.0 / inverseFlattening;
//    }
//}
//----------------------------------------------------------------------------------------------------------------------
Geographic::Geographic() : Lat( 0.0 ), Lon( 0.0 )
{}

Geographic::Geographic( double lat, double lon ) : Lat( lat ), Lon( lon )
{}
//----------------------------------------------------------------------------------------------------------------------
Geodetic::Geodetic() : Geographic( 0.0, 0.0 ), Height( 0.0 )
{}

Geodetic::Geodetic( double lat, double lon, double h ) : Geographic( lat, lon ), Height( h )
{}
//----------------------------------------------------------------------------------------------------------------------
RAD::RAD() : R( 0.0 ), Az( 0.0 ), AzEnd( 0.0 )
{}

RAD::RAD( double r, double az, double azEnd ) : R( r ), Az( az ), AzEnd( azEnd )
{}
//----------------------------------------------------------------------------------------------------------------------
XYZ::XYZ() : X( 0.0 ), Y( 0.0 ), Z( 0.0 )
{}

XYZ::XYZ( double x, double y, double z ) : X( x ), Y( y ), Z( z )
{}
//----------------------------------------------------------------------------------------------------------------------
ENU::ENU() : E( 0.0 ), N( 0.0 ), U( 0.0 )
{}

ENU::ENU( double e, double n, double u ) : E( e ), N( n ), U( u )
{}
//----------------------------------------------------------------------------------------------------------------------
AER::AER() : A( 0.0 ), E( 0.0 ), R( 0.0 )
{}

AER::AER( double a, double e, double r ) : A( a ), E( e ), R( r )
{}
//----------------------------------------------------------------------------------------------------------------------
void GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double latEnd, double lonEnd, double &d, double &az, double &azEnd )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    double f = ellipsoid.F();

    // По умолчанию Радианы:
    double _latStart = latStart;
    double _lonStart = lonStart;
    double _latEnd = latEnd;
    double _lonEnd = lonEnd;

    // При необходимости переведем входные данные в Радианы:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _latStart *= Convert::DgToRd;
            _lonStart *= Convert::DgToRd;
            _latEnd *= Convert::DgToRd;
            _lonEnd *= Convert::DgToRd;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    if( Compare::AreEqualAbs( a, b ) ) { // При расчете на сфере используем упрощенные формулы
        // Azimuth
        double fact1, fact2, fact3;
        fact1 = std::cos( _latEnd ) * std::sin( _lonEnd - _lonStart );
        fact2 = std::cos( _latStart ) * std::sin( _latEnd );
        fact3 = std::sin( _latStart ) * std::cos( _latEnd ) * std::cos( _lonEnd - _lonStart );
        az = Convert::AngleTo360( std::atan2( fact1, fact2 - fact3 ), Units::AU_Radian ); // [рад] - Прямой азимут в начальной точке

        // ReverseAzimuth
        fact1 = std::cos( _latStart ) * std::sin( _lonEnd - _lonStart );
        fact2 = std::cos( _latStart ) * std::sin( _latEnd ) * std::cos( _lonEnd - _lonStart );
        fact3 = std::sin( _latStart ) * std::cos( _latEnd );
        azEnd = Convert::AngleTo360( ( std::atan2( fact1, fact2 - fact3 ) ), Units::AU_Radian ); // [рад] - Прямой азимут в конечной точке

        // Distance
        double temp1, temp2, temp3;
        temp1 = std::sin( _latStart ) * std::sin( _latEnd );
        temp2 = std::cos( _latStart ) * std::cos( _latEnd ) * std::cos( _lonEnd - _lonStart );
        temp3 = temp1 + temp2;        
        d = std::acos( temp3 ) * a ; // [м]
    } else { // Для эллипсоида используем формулы Винсента
        double L = _lonEnd - _lonStart;

        double U1 = std::atan( ( 1.0 - f ) * std::tan( _latStart ) );
        double U2 = std::atan( ( 1.0 - f ) * std::tan( _latEnd ) );

        double sinU1 = std::sin( U1 );
        double cosU1 = std::cos( U1 );
        double sinU2 = std::sin( U2 );
        double cosU2 = std::cos( U2 );

        // eq. 13
        double lambda = L;
        double lambda_new = 0.0;
        int iterLimit = 100;

        double sinSigma = 0.0;
        double cosSigma = 0.0;
        double sigma = 0.0;
        double sinAlpha = 0.0;
        double cosSqAlpha = 0.0;
        double cos2SigmaM = 0.0;
        double c = 0.0;
        double sinLambda = 0.0;
        double cosLambda = 0.0;

        do {
            sinLambda = std::sin( lambda );
            cosLambda = std::cos( lambda );

            // eq. 14
            sinSigma = std::sqrt( ( ( cosU2 * sinLambda ) * ( cosU2 * sinLambda ) +
                ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) * ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) ) );
            if( Compare::IsZeroAbs( sinSigma ) ) { // co-incident points
                d = 0.0;
                az = 0.0;
                azEnd = 0.0;
                return;
            }

            // eq. 15
            cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

            // eq. 16
            sigma = std::atan2( sinSigma, cosSigma );

            // eq. 17    Careful!  sin2sigma might be almost 0!
            sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
            cosSqAlpha = 1 - sinAlpha * sinAlpha;

            // eq. 18    Careful!  cos2alpha might be almost 0!
            cos2SigmaM = cosSigma - 2.0 * sinU1 * sinU2 / cosSqAlpha;

            if( std::isnan( cos2SigmaM ) ) {
                cos2SigmaM = 0; // equatorial line: cosSqAlpha = 0
            }

            // eq. 10
            c = ( f / 16.0 ) * cosSqAlpha * ( 4.0 + f * ( 4.0 - 3.0 * cosSqAlpha ) );

            lambda_new = lambda;

            // eq. 11 (modified)
            lambda = L + ( 1.0 - c ) * f * sinAlpha *
                ( sigma + c * sinSigma * ( cos2SigmaM + c * cosSigma * ( -1.0 + 2.0 * cos2SigmaM * cos2SigmaM ) ) );

        } while( std::abs( ( lambda - lambda_new ) / lambda ) > 1.0e-15 && --iterLimit > 0 ); // see how much improvement we got

        double uSq = cosSqAlpha * ( a * a - b * b ) / ( b * b );

        // eq. 3
        double A = 1 + uSq / 16384.0 * ( 4096.0 + uSq * ( -768.0 + uSq * ( 320.0 - 175.0 * uSq ) ) );

        // eq. 4
        double B = uSq / 1024.0 * ( 256.0 + uSq * ( -128.0 + uSq * ( 74.0 - 47.0 * uSq ) ) );

        // eq. 6
        double deltaSigma = B * sinSigma *
            ( cos2SigmaM + ( B / 4.0 ) * ( cosSigma * ( -1.0 + 2.0 * cos2SigmaM * cos2SigmaM ) -
            ( B / 6.0 ) * cos2SigmaM * ( -3.0 + 4.0 * sinSigma * sinSigma ) * ( -3.0 + 4.0 * cos2SigmaM * cos2SigmaM ) ) );

        // eq. 19        
        d = b * A * ( sigma - deltaSigma ); // [m]

        // eq. 20
        az = Convert::AngleTo360( std::atan2( ( cosU2 * sinLambda ),
            ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) ), Units::TAngleUnit::AU_Radian ); // Прямой азимут в начальной точке, [рад]

        // eq. 21
        azEnd = Convert::AngleTo360( std::atan2( ( cosU1 * sinLambda ), ( -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda ) ), Units::TAngleUnit::AU_Radian ); // Прямой азимут в конечной точке, [рад]
    }
    // az, azEnd, d сейчас в радианах и метрах соответственно

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            az *= Convert::RdToDg;
            azEnd *= Convert::RdToDg;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            d *= 0.001;
            break;
        }
        default:
            assert( false );
    }
    return;
}

RAD GEOtoRAD(const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geographic &start, const Geographic &end )
{
    double d, az, azEnd;
    GEOtoRAD( ellipsoid, rangeUnit, angleUnit, start.Lat, start.Lon, end.Lat, end.Lon, d, az, azEnd );
    return RAD( d, az, azEnd );
}

//----------------------------------------------------------------------------------------------------------------------
void RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double d, double az, double &latEnd, double &lonEnd, double &azEnd )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    double f = ellipsoid.F();

    // по умолчанию Метры-Радианы:
    double _latStart = latStart;    // [рад]
    double _lonStart = lonStart;    // [рад]
    double _d = d;                  // [м]
    double _az = az;                // [рад]

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _latStart *= Convert::DgToRd;
            _lonStart *= Convert::DgToRd;
            _az *= Convert::DgToRd;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _d *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    if( Compare::AreEqualAbs(a, b) ) { // При расчете на сфере используем упрощенные формулы
        _d = _d / a; // Нормирование
        // latitude
        double temp1, temp2, temp3;
        temp1 = std::sin( _latStart ) * std::cos( _d );
        temp2 = std::cos( _latStart ) * std::sin( _d ) * std::cos( _az );
        latEnd = std::asin( temp1 + temp2 ); // [рад]

        // longitude
        temp1 = std::sin( _d ) * std::sin( _az );
        temp2 = std::cos( _latStart ) * std::cos( _d );
        temp3 = std::sin( _latStart ) * std::sin( _d ) * std::cos( _az );
        lonEnd = _lonStart + std::atan2( temp1, temp2 - temp3 ); // [рад]

        // final bearing
        temp1 = std::cos( _latStart ) * std::sin( _az );
        temp2 = std::cos( _latStart ) * std::cos( _d ) * std::cos( _az );
        temp3 = std::sin( _latStart ) * std::sin( _d );
        azEnd = Convert::AngleTo360( std::atan2( temp1, temp2 - temp3 ), Units::TAngleUnit::AU_Radian ); // [рад] - Прямой азимут в конечной точке
    } else { // Для эллипсоида используем формулы Винсента
        double cosAlpha1 = std::cos( _az );
        double sinAlpha1 = std::sin( _az );
        double s = _d; // distance [m]
        double tanU1 = ( 1.0 - f ) * std::tan( _latStart );
        double cosU1 = 1.0 / std::sqrt( ( 1.0 + tanU1 * tanU1 ) );
        double sinU1 = tanU1 * cosU1;

        // eq. 1
        double sigma1 = std::atan2( tanU1, cosAlpha1 );

        // eq. 2
        double sinAlpha = cosU1 * sinAlpha1;
        double cosSqAlpha = 1 - sinAlpha * sinAlpha;
        double uSq = cosSqAlpha * ( a * a - b * b ) / ( b * b );

        // eq. 3
        double A = 1.0 + ( uSq / 16384.0 ) * ( 4096.0 + uSq * ( -768.0 + uSq * ( 320.0 - 175.0 * uSq ) ) );

        // eq. 4
        double B = ( uSq / 1024.0 ) * ( 256.0 + uSq * ( -128.0 + uSq * ( 74.0 - 47.0 * uSq ) ) );

        // iterate until there is a negligible change in sigma
        double sOverbA = s / ( b * A );
        double sigma = sOverbA;
        double prevSigma = sOverbA;
        double cos2SigmaM = 0.0;
        double sinSigma = 0.0;
        double cosSigma = 0.0;
        double deltaSigma = 0.0;

        int iterations = 0;

        while( true ) {
            // eq. 5
            cos2SigmaM = std::cos( 2.0 * sigma1 + sigma );
            sinSigma = std::sin( sigma );
            cosSigma = std::cos( sigma );

            // eq. 6
            deltaSigma = B * sinSigma * ( cos2SigmaM +
                ( B / 4.0 ) * ( cosSigma * ( -1.0 + 2.0 * cos2SigmaM * cos2SigmaM ) -
                ( B / 6.0 ) * cos2SigmaM * ( -3.0 + 4.0 * sinSigma * sinSigma ) * ( -3.0 + 4.0 * cos2SigmaM * cos2SigmaM ) ) );

            // eq. 7
            sigma = sOverbA + deltaSigma;

            // break after converging to tolerance
            if( std::abs( sigma - prevSigma ) < 1.0e-15 || std::isnan( std::abs( sigma - prevSigma ) ) ) {
                break;
            }
            prevSigma = sigma;

            iterations++;
            if( iterations > 1000 ) {
                break;
            }
        }
        cos2SigmaM = std::cos( 2.0 * sigma1 + sigma );
        sinSigma = std::sin( sigma );
        cosSigma = std::cos( sigma );

        double tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1;

        // eq. 8
        latEnd = std::atan2( sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
            ( 1.0 - f ) * std::sqrt( ( sinAlpha * sinAlpha + tmp * tmp) ) ); // [рад]

        // eq. 9
        double lambda = std::atan2( ( sinSigma * sinAlpha1 ), ( cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1 ) );

        // eq. 10
        double c = ( f / 16.0 ) * cosSqAlpha * ( 4.0 + f * ( 4.0 - 3.0 * cosSqAlpha ) );

        // eq. 11
        double L = lambda - ( 1.0 - c ) * f * sinAlpha * ( sigma + c * sinSigma *
            ( cos2SigmaM + c * cosSigma * ( -1.0 + 2.0 * cos2SigmaM * cos2SigmaM ) ) );

        //double phi = ( _lonStart + L + 3 * PI ) % ( 2 * PI ) - PI;   // to -180.. 180 original!
        //lonEnd = ( _lonStart + L ) * RdToDgD;   // [град] my
        lonEnd = _lonStart + L;   // [рад] my

        // eq. 12
        double alpha2 = std::atan2( sinAlpha, -tmp ); // final bearing, if required
        azEnd = Convert::AngleTo360( alpha2, Units::TAngleUnit::AU_Radian ); // Прямой азимут в конечной точке, [рад]
    }
    // latEnd, lonEnd, azEnd сейчас в радианах

    // Проверим, нужен ли перевод:    
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            latEnd *= Convert::RdToDg;
            lonEnd *= Convert::RdToDg;
            azEnd *= Convert::RdToDg;
            break;
        }
        default:
            assert( false );
    }
    return;
}

Geographic RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geographic &start, const RAD &rad, double &azEnd )
{
    double latEnd, lonEnd;
    RADtoGEO( ellipsoid, rangeUnit, angleUnit,
        start.Lat, start.Lon, rad.R, rad.Az, latEnd, lonEnd, azEnd );
    return Geographic( latEnd, lonEnd );
}
//----------------------------------------------------------------------------------------------------------------------
void GEOtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat, double lon, double h, double &x, double &y, double &z )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
//    double b = ellipsoid.B();

    // по умолчанию Метры-Радианы:
    double _lat = lat; // [рад]
    double _lon = lon; // [рад]
    double _h = h;     // [м]

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat *= Convert::DgToRd;
            _lon *= Convert::DgToRd;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце
    assert( !Compare::IsZeroAbs( a * a ) );
//    double es = 1.0 - ( ( b * b ) / ( a * a ) ); // e^2
    double es = ellipsoid.EccentricityFirstSquared();
    assert( ( 1.0 - ( es * ( std::sin( _lat ) * std::sin( _lat ) ) ) ) > 0 );
    double v = a / std::sqrt( 1.0 - ( es * ( std::sin( _lat ) * std::sin( _lat ) ) ) );

    x = ( v + _h ) * std::cos( _lat ) * std::cos( _lon );         // [м]
    y = ( v + _h ) * std::cos( _lat ) * std::sin( _lon );         // [м]
    z = ( ( ( 1.0 - es ) * v ) + _h ) * std::sin( _lat );    // [м]

    // x, y, x сейчас в метрах

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            x /= 1000.0;
            y /= 1000.0;
            z /= 1000.0;
            break;
        }
        default:
            assert( false );
    }
}

XYZ GEOtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic point )
{
    double x, y, z;
    GEOtoECEF( ellipsoid, rangeUnit, angleUnit, point.Lat, point.Lon, point.Height, x, y, z );
    return XYZ( x, y, z );
}
//----------------------------------------------------------------------------------------------------------------------
void ECEFtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double x, double y, double z, double &lat, double &lon, double &h )
{
    static const double AD_C = 1.0026000;   // Toms region 1 constant. ( h_min = -1e5 [m], h_max = 2e6 [m] ) - MOST CASES!
//    static const double AD_C = 1.00092592;  // Toms region 2 constant. ( h_min = 2e6 [m], h_max = 6e6 [m] )
//    static const double AD_C = 0.999250297; // Toms region 3 constant. ( h_min = 6e6 [m], h_max = 18e6 [m] )
//    static const double AD_C = 0.997523508; // Toms region 4 constant. ( h_min = 18e5 [m], h_max = 1e9 [m] )
    static const double COS_67P5 = 0.38268343236508977; // Cosine of 67.5 degrees

    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();

//    double es = 1.0 - ( b * b ) / ( a * a );    // Eccentricity squared : (a^2 - b^2)/a^2
//    double ses = ( a * a ) / ( b * b ) - 1.0;   // Second eccentricity squared : (a^2 - b^2)/b^2
    double es = ellipsoid.EccentricityFirstSquared();   // Eccentricity squared : (a^2 - b^2)/a^2
    double ses = ellipsoid.EccentricitySecondSquared(); // Second eccentricity squared : (a^2 - b^2)/b^2

    bool At_Pole = false; // indicates whether location is in polar region */

    double _x = x; // [м]
    double _y = y; // [м]
    double _z = z; // [м]

    // При необходимости переведем в Метры:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _x *= 1000.0;
            _y *= 1000.0;
            _z *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце
    double _lon = 0;
    double _lat = 0;
    double _h = 0;

    if( !Compare::IsZeroAbs( x ) ) { //if (x != 0.0)
        _lon = std::atan2( _y, _x );
    } else {
        if( _y > 0 ) {
            _lon = Consts::PI / 2;
        } else if ( _y < 0 ) {
            _lon = -Consts::PI * 0.5;
        } else {
            At_Pole = true;
            _lon = 0.0;
            if( _z > 0.0 ) {            // north pole
                _lat = Consts::PI_05;
            } else if ( _z < 0.0 ) {    // south pole
                _lat = -Consts::PI_05;
            } else {                    // center of earth
                //lat = PI * 0.5 * RdToDgD; // [град]
                lat = Consts::PI_05; // [рад] TODO: Как тут улучшить?
                lon = 0;
                h = -b;
                return;
            }
        }
    }
    double W2 = _x * _x + _y * _y; // Square of distance from Z axis
    assert( W2 > 0 );
    double W = std::sqrt( W2 ); // distance from Z axis
    double T0 = _z * AD_C; // initial estimate of vertical component
    assert( ( T0 * T0 + W2 ) > 0 );
    double S0 = std::sqrt( T0 * T0 + W2 ); //initial estimate of horizontal component
    double Sin_B0 = T0 / S0; //std::sin(B0), B0 is estimate of Bowring aux variable
    double Cos_B0 = W / S0; //std::cos(B0)
    double Sin3_B0 = Sin_B0 * Sin_B0 * Sin_B0;//Math.Pow(Sin_B0, 3);
    double T1 = _z + b * ses * Sin3_B0; //corrected estimate of vertical component
    double Sum = W - a * es * Cos_B0 * Cos_B0 * Cos_B0; //numerator of std::cos(phi1)
    assert( ( T1 * T1 + Sum * Sum ) > 0 );
    double S1 = std::sqrt( T1 * T1 + Sum * Sum ); //corrected estimate of horizontal component
    double Sin_p1 = T1 / S1; //std::sin(phi1), phi1 is estimated latitude
    double Cos_p1 = Sum / S1; //std::cos(phi1)
    assert( ( 1.0 - es * Sin_p1 * Sin_p1 ) > 0 );
    double Rn = a / std::sqrt( 1.0 - es * Sin_p1 * Sin_p1 ); //Earth radius at location
    if( Cos_p1 >= COS_67P5 ) {
        _h = W / Cos_p1 - Rn;
    } else if ( Cos_p1 <= -COS_67P5 ) {
        assert( !Compare::IsZeroAbs( Cos_p1 ) );
        _h = W / ( -Cos_p1 ) - Rn;
    } else {
        assert( !Compare::IsZeroAbs( Sin_p1 ) );
        _h = ( _z / Sin_p1 ) + Rn * ( es - 1.0 );
    }
    if( !At_Pole ) {
        assert( !Compare::IsZeroAbs( Cos_p1 ) );
        _lat = std::atan2( Sin_p1, Cos_p1 );
    }
    lat = _lat; // [рад]
    lon = _lon; // [рад]
    h = _h;     // [м]

    // LLH сейчас в радианах и метрах соответственно

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            lat *= Convert::RdToDg;
            lon *= Convert::RdToDg;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            h /= 1000.0;
            break;
        }
        default:
            assert( false );
    }
}

Geodetic ECEFtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    XYZ &point )
{
    double lat, lon, h;
    ECEFtoGEO( ellipsoid, rangeUnit, angleUnit, point.X, point.Y, point.Z, lat, lon, h );
    return Geodetic( lat, lon, h );
}
//----------------------------------------------------------------------------------------------------------------------
double XYZtoDistance( double x1, double y1, double z1, double x2, double y2, double z2 )
{
    double res = ( ( ( x1 - x2 ) * ( x1 - x2 ) ) +
                   ( ( y1 - y2 ) * ( y1 - y2 ) ) +
                   ( ( z1 - z2 ) * ( z1 - z2 ) ) );
    assert( res >= 0 );
    res = std::sqrt( res );
    return res;
}

double XYZtoDistance( const XYZ &point1, const XYZ &point2 )
{
    return XYZtoDistance( point1.X, point1.Y, point1.Z, point2.X, point2.Y, point2.Z );
}
//----------------------------------------------------------------------------------------------------------------------
void ECEF_offset( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat1, double lon1, double h1, double lat2, double lon2, double h2, double &dX, double &dY, double &dZ )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    //double f = el.F();

    // по умолчанию Метры-Радианы:
    double _lat1 = lat1;
    double _lon1 = lon1;
    double _h1 = h1;
    double _lat2 = lat2;
    double _lon2 = lon2;
    double _h2 = h2;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat1 *= Convert::DgToRd;
            _lon1 *= Convert::DgToRd;
            _lat2 *= Convert::DgToRd;
            _lon2 *= Convert::DgToRd;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h1 *= 1000.0;
            _h2 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double s1 = std::sin( lat1 );
    double c1 = std::cos( lat1 );

    double s2 = std::sin( lat2 );
    double c2 = std::cos( lat2 );

    double p1 = c1 * std::cos( lon1 );
    double p2 = c2 * std::cos( lon2 );

    double q1 = c1 * std::sin( lon1 );
    double q2 = c2 * std::sin( lon2 );

    if( Compare::AreEqualAbs( a, b ) ) { // Сфера
        dX = a * ( p2 - p1 ) + ( h2 * p2 - h1 * p1 );
        dY = a * ( q2 - q1 ) + ( h2 * q2 - h1 * q1 );
        dZ = a * ( s2 - s1 ) + ( h2 * s2 - h1 * s1 );
    } else { // Эллипсоид
        double e2 = std::pow( ellipsoid.EccentricityFirst(), 2 ); // Квадрат 1-го эксцентриситета эллипсоида

        double w1 = 1.0 / std::sqrt( 1.0 - e2 * s1 * s1 );
        double w2 = 1.0 / std::sqrt( 1.0 - e2 * s2 * s2 );

        dX = a * ( p2 * w2 - p1 * w1 ) + ( h2 * p2 - h1 * p1 );
        dY = a * ( q2 * w2 - q1 * w1 ) + ( h2 * q2 - h1 * q1 );
        dZ = ( 1.0 - e2 ) * a * ( s2 * w2 - s1 * w1 ) + ( h2 * s2 - h1 * s1 );
    }
    // dX dY dZ сейчас в метрах

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            dX *= 0.001;
            dY *= 0.001;
            dZ *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

XYZ ECEF_offset( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic &point1, const Geodetic &point2 )
{
    double x, y, z;
    ECEF_offset( ellipsoid, rangeUnit, angleUnit,
        point1.Lat, point1.Lon, point1.Height, point2.Lat, point2.Lon, point2.Height, x, y, z );
    return XYZ( x, y, z );
}
//----------------------------------------------------------------------------------------------------------------------
void ECEFtoENU( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double dX, double dY, double dZ, double lat, double lon, double &xEast, double &yNorth, double &zUp )
{
    // по умолчанию Метры-Радианы:
    double _lat = lat;
    double _lon = lon;
    double _dX = dX;
    double _dY = dY;
    double _dZ = dZ;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat *= Convert::DgToRd;
            _lon *= Convert::DgToRd;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _dX *= 1000.0;
            _dY *= 1000.0;
            _dZ *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double cosPhi = std::cos( _lat );
    double sinPhi = std::sin( _lat );
    double cosLambda = std::cos( _lon );
    double sinLambda = std::sin( _lon );

    double t = ( cosLambda * dX ) + ( sinLambda * dY );
    xEast = ( -sinLambda * dX ) + ( cosLambda * dY );

    zUp    =  ( cosPhi * t ) + ( sinPhi * dZ );
    yNorth =  ( -sinPhi * t ) + ( cosPhi * dZ );

    // xEast yNorth zUp сейчас в метрах

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            xEast *= 0.001;
            yNorth *= 0.001;
            zUp *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

ENU ECEFtoENU( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const XYZ &shift, const Geographic &point )
{
    double e, n, u;
    ECEFtoENU( rangeUnit, angleUnit, shift.X, shift.Y, shift.Z, point.Lat, point.Lon, e, n, u );
    return ENU( e, n, u );
}
//----------------------------------------------------------------------------------------------------------------------
void ENUtoAER( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double xEast, double yNorth, double zUp, double &az, double &elev, double &slantRange )
{
    // по умолчанию Метры-Радианы:
    double _xEast = xEast;
    double _yNorth = yNorth;
    double _zUp = zUp;

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            _xEast *=  1000.0;
            _yNorth *= 1000.0;
            _zUp *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

//    r = std::sqrt( ( xEast * xEast ) + ( yNorth * yNorth ) ); // dangerous
    double r = std::hypot( xEast, yNorth ); // C++11 style

//    slantRange = sqrt( ( r * r ) + ( zUp * zUp ) ); // dangerous
    slantRange = std::hypot( r, zUp ); // C++11 style
    elev = std::atan2( zUp, r );
    az = Convert::AngleTo360( std::atan2( xEast, yNorth ), Units::TAngleUnit::AU_Radian );

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            slantRange *= 0.001;
            break;
        }
        default:
            assert( false );
    }
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            az *= Convert::RdToDg;
            elev *= Convert::RdToDg;
        }
    }
}

AER ENUtoAER( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit, const ENU &point )
{
    double a, e, r;
    ENUtoAER( rangeUnit, angleUnit, point.E, point.N, point.U, a, e, r );
    return AER( a, e, r );
}
//----------------------------------------------------------------------------------------------------------------------
void GEOtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat, double lon, double h, double lat0, double lon0, double h0, double &xEast, double &yNorth, double &zUp )
{
    double dX, dY, dZ;
    ECEF_offset( ellipsoid, rangeUnit, angleUnit, lat0, lon0, h0, lat, lon, h, dX, dY, dZ );
    ECEFtoENU( rangeUnit, angleUnit, dX, dY, dZ, lat0, lon0, xEast, yNorth, zUp );
}

ENU GEOtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic &point, const Geodetic &anchor )
{
    double e, n, u;
    GEOtoENU( ellipsoid, rangeUnit, angleUnit,
        point.Lat, point.Lon, point.Height, anchor.Lat, anchor.Lon, anchor.Height, e, n, u );
    return ENU( e, n, u );
}
//----------------------------------------------------------------------------------------------------------------------
void GEOtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat1, double lon1, double h1, double lat2, double lon2, double h2, double &az, double &elev, double &slantRange )
{
    // по умолчанию Метры-Радианы:
    double _lat1 = lat1;
    double _lon1 = lon1;
    double _h1 = h1;
    double _lat2 = lat2;
    double _lon2 = lon2;
    double _h2 = h2;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat1 *= Convert::DgToRd;
            _lon1 *= Convert::DgToRd;
            _lat2 *= Convert::DgToRd;
            _lon2 *= Convert::DgToRd;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h1 *= 1000.0;
            _h2 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double xEast, yNorth, zUp;
    GEOtoENU( ellipsoid, rangeUnit, angleUnit, lat1, lon1, h1, lat2, lon2, h2, xEast, yNorth, zUp );
    ENUtoAER( rangeUnit, angleUnit, xEast, yNorth, zUp, az, elev, slantRange );
}

AER GEOtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     const Geodetic &point1, const Geodetic &point2 )
{
    double a, e, r;
    GEOtoAER( ellipsoid, rangeUnit, angleUnit,
        point1.Lat, point1.Lon, point1.Height, point2.Lat, point2.Lon, point2.Height, a, e, r );
    return AER( a, e, r );
}
//----------------------------------------------------------------------------------------------------------------------
double CosAngleBetweenVectors( double x1, double y1, double z1, double x2, double y2, double z2 )
{
    // Исходя из формулы косинуса угла между векторами:
    double a1 =  x1 * x2;
    double a2 =  y1 * y2;
    double a3 =  z1 * z2;
    double b1 = std::sqrt( ( x1 * x1 ) + ( y1 * y1 ) + ( z1 * z1 ) );
    double b2 = std::sqrt( ( x2 * x2 ) + ( y2 * y2 ) + ( z2 * z2 ) );

    double res = ( a1 + a2 + a3 ) / ( b1 * b2 );
    return res;
}

double CosAngleBetweenVectors( const XYZ &point1, const XYZ &point2 )
{
    return CosAngleBetweenVectors( point1.X, point1.Y, point1.Z, point2.X, point2.Y, point2.Z );
}

//----------------------------------------------------------------------------------------------------------------------
double AngleBetweenVectors( double x1, double y1, double z1, double x2, double y2, double z2 )
{
    return std::acos( CosAngleBetweenVectors( x1, y1, z1, x2, y2, z2 ) );
}

double AngleBetweenVectors( const XYZ &vec1, const XYZ &vec2 )
{
    return std::acos( CosAngleBetweenVectors( vec1, vec2 ) );
}
//----------------------------------------------------------------------------------------------------------------------
void VectorFromTwoPoints( double x1, double y1, double z1, double x2, double y2, double z2, double &xV, double &yV, double &zV )
{
    xV = x2 - x1;
    yV = y2 - y1;
    zV = z2 - z1;
}

XYZ VectorFromTwoPoints( const XYZ &point1, const XYZ &point2 )
{
    XYZ result;
    VectorFromTwoPoints( point1.X, point1.Y, point1.Z, point2.X, point2.Y, point2.Z, result.X, result.Y, result.Z );
    return result;
}

}
}
/// \}
