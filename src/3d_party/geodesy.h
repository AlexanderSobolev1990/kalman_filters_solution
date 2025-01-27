//----------------------------------------------------------------------------------------------------------------------
///
/// \file       geodesy.h
/// \brief      Земные эллипсоиды, геодезические задачи (персчеты координат)
/// \date       06.11.19 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_GEODESY_H
#define SPML_GEODESY_H

// System includes:
#include <cassert>
#include <string>
#include <vector>

// SPML includes:
#include <compare.h>
#include <convert.h>
#include <units.h>

namespace SPML /// Специальная библиотека программных модулей (СБПМ)
{
namespace Geodesy /// Геодезические функции и функции перевода координат
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Земной эллипсоид
///
class CEllipsoid
{
public:
    std::string Name() const; ///< Название эллипсоида
    double A() const; ///< Большая полуось эллипсоида (экваториальный радиус), [м]
    double B() const; ///< Малая полуось эллипсоида (полярный радиус), [м]
    double F() const; ///< Cжатие f = ( a - b ) / a
    double Invf() const; ///< Обратное сжатие Invf = a / ( a - b )
    double EccentricityFirst() const; ///< Первый эксцентриситет эллипсоида \f$ e1 = \sqrt{ a^2 - b^2 } / a \f$
    double EccentricityFirstSquared() const; ///< Квадрат первого эксцентриситета эллипсоида \f$ es1 = 1 - ( b^2 / a^2 ) \f$
    double EccentricitySecond() const; ///< Второй эксцентриситет эллипсоида \f$ e2 = \sqrt{ a^2 - b^2 } / b \f$
    double EccentricitySecondSquared() const; ///< Квадрат второго эксцентриситета эллипсоида \f$ es2 = ( a^2 / b^2 ) - 1 \f$

    ///
    /// \brief Конструктор по умолчанию
    ///
    CEllipsoid();

    ///
    /// \brief Параметрический конструктор эллипсоида
    /// \param[in] ellipsoidName - название эллипсоида
    /// \param[in] semiMajorAxis - большая полуось (экваториальный радиус)
    /// \param[in] semiMinorAxis - малая полуось (полярный радиус)
    /// \param[in] inverseFlattening - обратное сжатие invf = a / ( a - b )
    /// \param[in] isInvfDef - обратное сжатие задано (малая полуось расчитана из большой и обратного сжатия)
    ///
    CEllipsoid( std::string ellipsoidName, double semiMajorAxis, double semiMinorAxis, double inverseFlattening, bool isInvfDef );

private:
    std::string name; ///< Название эллипсоида
    double a; ///< Большая полуось (экваториальный радиус), [м]
    double b; ///< Малая полуось (полярный радиус), [м]
    double invf; ///< Обратное сжатие invf = a / ( a - b )
    double f; ///< Сжатие f = ( a - b ) / a
};

namespace Ellipsoids /// Земные эллипсоиды
{
//------------------------------------------------------------------------------------------------------------------
//
//                                              Земные эллипсоиды:
//
//  1) Сфера радиусом 6371000.0 [м], https://epsg.io/7035-ellipsoid
//  2) Сфера радиусом 6378000.0 [м]
//  3) Эллипсоид WGS84, https://epsg.io/7030-ellipsoid
//  4) Эллипсоид GRS80, https://epsg.io/7019-ellipsoid
//  5) Эллипсоид ПЗ-90, https://epsg.io/7054-ellipsoid
//  6) Эллипсоид Крассовского, https://epsg.io/7024-ellipsoid
//  7) Сфера радиусом большой полуоси эллипсоида Красовского 1940 (6378245.0 [м])
//

///
/// \brief Сфера радиусом 6371000.0 [м] (EPSG:7035)
/// \details Обратное сжатие - бесконечность
///
static CEllipsoid Sphere6371()
{
    return CEllipsoid( "Сфера радиусом 6371000.0 [м] (EPSG:7035)", 6371000.0, 6371000.0, 0.0, false );
}

///
/// \brief Сфера радиусом 6378000.0 [м]
/// \details Обратное сжатие - бесконечность
///

static CEllipsoid Sphere6378()
{
    return CEllipsoid( "Сфера радиусом 6378000.0 [м]", 6378000.0, 6378000.0, 0.0, false );
}

///
/// \brief Эллипсоид WGS84 (EPSG:7030)
/// \details Главная полуось 6378137.0, обратное сжатие 298.257223563
///
static CEllipsoid WGS84()
{
    return CEllipsoid( "Эллипсоид WGS84 (EPSG:7030)", 6378137.0, 0.0, 298.257223563, true );
}

///
/// \brief Эллипсоид GRS80 (EPSG:7019)
/// \details Главная полуось 6378137.0, обратное сжатие 298.257222101
///
static CEllipsoid GRS80()
{
    return CEllipsoid( "Эллипсоид GRS80 (EPSG:7019)", 6378137.0, 0.0, 298.257222101, true );
}

///
/// \brief Эллипсоид ПЗ-90 (EPSG:7054)
/// \details Главная полуось 6378136.0, обратное сжатие 298.257839303
///
static CEllipsoid PZ90()
{
    return CEllipsoid( "Эллипсоид ПЗ-90 (EPSG:7054)", 6378136.0, 0.0, 298.257839303, true );
}

///
/// \brief Эллипсоид Красовского 1940 (EPSG:7024)
/// \details Главная полуось 6378245.0, обратное сжатие 298.3
///
static CEllipsoid Krassowsky1940()
{
    return CEllipsoid( "Эллипсоид Красовского 1940 (EPSG:7024)", 6378245.0, 0.0, 298.3, true );
}

///
/// \brief Сфера радиусом большой полуоси эллипсоида Красовского 1940 (EPSG:7024)
/// \details Обратное сжатие - бесконечность
///
static CEllipsoid SphereKrassowsky1940()
{
    return CEllipsoid( "Сфера радиусом большой полуоси эллипсоида Красовского 1940 (EPSG:7024)", 6378245.0, 6378245.0, 0.0, false );
}

///
/// \brief Возвращает доступные предопределенные эллипсоидоы
/// \return Вектор предопределенных эллипсоидов
///
[[maybe_unused]]
static const __attribute__ ((unused)) std::vector<CEllipsoid> GetPredefinedEllipsoids()
{
    return std::vector<CEllipsoid>{
        Sphere6371(),
        Sphere6378(),
        WGS84(),
        GRS80(),
        PZ90(),
        Krassowsky1940(),
        SphereKrassowsky1940()
    };
}

} // end namespace Ellipsoids

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Географические координаты (широта, долгота)
///
struct Geographic
{
    double Lat; ///< Широта
    double Lon; ///< Долгота

    Geographic(); ///< Конструктор по умолчанию

    ///
    /// \brief Параметрический конструктор
    /// \param lat - широта
    /// \param lon - долгота
    ///
    Geographic( double lat, double lon );
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Геодезические координаты (широта, долгота, высота)
///
struct Geodetic : public Geographic
{
    double Height; ///< Высота

    Geodetic(); ///< Конструктор по умолчанию

    ///
    /// \brief Параметрический конструктор
    /// \param lat - широта
    /// \param lon - долгота
    /// \param h - высота
    ///
    Geodetic( double lat, double lon, double h );
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Радиолокационные координаты (расстояние по ортодроме, азимут, конечный азимут)
/// \details Имеют два азимута: Az - это начальный азимут в точке наблюдения, AzEnd - конечный азимут (по ортодроме) на дальности R
///
struct RAD
{
    double R;       /// Дальность
    double Az;      /// Азимут в точке наблюдения
    double AzEnd;   /// Азимут в точке объекта

    RAD(); ///< Конструктор по умолчанию

    ///
    /// \brief Параметрический конструктор
    /// \param r - дальность
    /// \param az - азимут
    /// \param azEnd - конечный азимут
    ///
    RAD( double r, double az, double azEnd );
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief 3D декартовы ортогональные координаты (X, Y, Z)
///
struct XYZ
{
    double X; ///< X координата
    double Y; ///< Y координата
    double Z; ///< Z координата

    XYZ(); ///< Конструктор по умолчанию

    ///
    /// \brief Параметрический конструктор
    /// \param x - X координата
    /// \param y - Y координата
    /// \param z - Z координата
    ///
    XYZ( double x, double y, double z );
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Координаты ENU (East-North-Up)
///
struct ENU
{
    double E; ///< East координата
    double N; ///< North координата
    double U; ///< Up координата

    ENU(); ///< Конструктор по умолчанию

    ///
    /// \brief Параметрический конструктор
    /// \param e - East координата
    /// \param n - North координата
    /// \param u - Up координата
    ///
    ENU( double e, double n, double u );
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Координаты AER (Azimuth-Elevation-Range, Азимут-Угол места-Дальность)
///
struct AER
{
    double A; ///< Азимут
    double E; ///< Угол места
    double R; ///< Дальность

    AER(); ///< Конструктор по умолчанию

    ///
    /// \brief Параметрический конструктор
    /// \param a - азимут
    /// \param e - угол места
    /// \param r - дальность
    ///
    AER( double a, double e, double r );
};

//----------------------------------------------------------------------------------------------------------------------
//
//                                          Функции пересчета координат
//
[[maybe_unused]]
static double dummy_double; // Заглушка для списка параметров функций без перегрузки

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет географических координат в радиолокационные (Обратная геодезическая задача)
/// \details Расчет на эллипсоиде по формулам Винсента:
/// \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
/// application of nested equations". Survey Review. XXIII (176): 88–93.
/// \n Расчет на сфере:
/// \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] latStart - широта начальной точки
/// \param[in] lonStart - долгота начальной точки
/// \param[in] latEnd - широта конечной точки
/// \param[in] lonEnd - долгота конечной точки
/// \param[out] d - расстояние между начальной и конечной точками по ортодроме
/// \param[out] az - азимут из начальной точки на конечную
/// \param[out] azEnd - азимут в конечной точке
///
void GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double latEnd, double lonEnd, double &d, double &az, double &azEnd = dummy_double );

///
/// \brief Пересчет географических координат в радиолокационные (Обратная геодезическая задача)
/// \details Расчет на эллипсоиде по формулам Винсента:
/// \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
/// application of nested equations". Survey Review. XXIII (176): 88–93.
/// \n Расчет на сфере:
/// \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] start - начальная точка
/// \param[in] end - конечная точка
/// \return Радиолокационные координаты (расстояние между начальной и конечной точками по ортодроме,
/// азимут из начальной точки на конечную, азимут в конечной точке)
///
RAD GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geographic &start, const Geographic &end );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет радиолокационных координат в географические (Прямая геодезическая задача)
/// \details Расчет на эллипсоиде по формулам Винсента:
/// \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
/// application of nested equations". Survey Review. XXIII (176): 88–93.
/// \n Расчет на сфере:
/// \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] latStart - широта начальной точки
/// \param[in] lonStart - долгота начальной точки
/// \param[in] d - расстояние между начальной и конечной точками по ортодроме
/// \param[in] az - азимут из начальной точки на конечную
/// \param[out] latEnd - широта конечной точки
/// \param[out] lonEnd - долгота конечной точки
/// \param[out] azEnd - прямой азимут в конечной точке
///
void RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double d, double az, double &latEnd, double &lonEnd, double &azEnd = dummy_double );

///
/// \brief Пересчет радиолокационных координат в географические (Прямая геодезическая задача)
/// \details Расчет на эллипсоиде по формулам Винсента:
/// \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
/// application of nested equations". Survey Review. XXIII (176): 88–93.
/// \n Расчет на сфере:
/// \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] start - географические координаты начальной точки
/// \param[in] rad - радиолокационные координаты пути (расстояние между начальной и конечной точками по ортодроме,
/// азимут из начальной точки на конечную)
/// \param[out] azEnd - прямой азимут в конечной точке
/// \return Географические координаты конечной точки
///
Geographic RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geographic &start, const RAD &rad, double &azEnd = dummy_double );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет широты, долготы, высоты в декартовые геоцентрические координаты
/// \details При отсутствии высоты или расположении точки на поверхности эллипсоида, задать координату высоты h = 0.
/// Декартовые геоцентрические координаты (ECEF):
/// ось X - через пересечение гринвичского меридиана и экватора,
/// ось Y - через пересечение меридиана 90 [град] восточной долготы и экватора,
/// ось Z - через северный полюс.
/// Источник - Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 191
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] lat - широта точки
/// \param[in] lon - долгота точки
/// \param[in] h - высота точки над поверхностью эллипсоида
/// \param[out] x - координата по оси Х
/// \param[out] y - координата по оси Y
/// \param[out] z - координата по оси Z
///
void GEOtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat, double lon, double h, double &x, double &y, double &z );

///
/// \brief Пересчет широты, долготы, высоты в декартовые геоцентрические координаты
/// \details При отсутствии высоты или расположении точки на поверхности эллипсоида, задать координату высоты h = 0.
/// Декартовые геоцентрические координаты (ECEF):
/// ось X - через пересечение гринвичского меридиана и экватора,
/// ось Y - через пересечение меридиана 90 [град] восточной долготы и экватора,
/// ось Z - через северный полюс.
/// Источник - Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 191
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point - геодезические координаты точки
/// \return ECEF координаты точки point
///
XYZ GEOtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет декартовых геоцентрических координат в широту, долготу, высоту
/// \details Декартовые геоцентрические координаты (ECEF):
/// ось X - через пересечение гринвичского меридиана и экватора,
/// ось Y - через пересечение меридиана 90 [град] восточной долготы и экватора,
/// ось Z - через северный полюс.
/// Источник - An Improved Algorithm for Geocentric to Geodetic Coordinate Conversion, Ralph Toms, Feb 1996.
/// UCRL-JC-123138
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] x - координата по оси Х
/// \param[in] y - координата по оси Y
/// \param[in] z - координата по оси Z
/// \param[out] lat - широта точки
/// \param[out] lon - долгота точки
/// \param[out] h - высота точки над поверхностью эллипсоида
///
void ECEFtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double x, double y, double z, double &lat, double &lon, double &h );

///
/// \brief Пересчет декартовых геоцентрических координат в широту, долготу, высоту
/// \details Декартовые геоцентрические координаты (ECEF):
/// ось X - через пересечение гринвичского меридиана и экватора,
/// ось Y - через пересечение меридиана 90 [град] восточной долготы и экватора,
/// ось Z - через северный полюс.
/// Источник - An Improved Algorithm for Geocentric to Geodetic Coordinate Conversion, Ralph Toms, Feb 1996.
/// UCRL-JC-123138
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point - точка в координатах ECEF
/// \return Геодезические (широта, долгота, высота) координаты точки point
///
Geodetic ECEFtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    XYZ &point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вычисление расстояния между точками в декартовых координатах
/// \details Вычисляется как d = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
/// \attention Единицы измерения выхода соответствуют единицам измерения входа
/// \param[in] x1 - координата первой точки по оси Х
/// \param[in] y1 - координата первой точки по оси Y
/// \param[in] z1 - координата первой точки по оси Z
/// \param[in] x2 - координата второй точки по оси Х
/// \param[in] y2 - координата второй точки по оси Y
/// \param[in] z2 - координата второй точки по оси Z
/// \return Расстояние между двумя точками в декартовых координатах
///
double XYZtoDistance( double x1, double y1, double z1, double x2, double y2, double z2 );

///
/// \brief Вычисление расстояния между точками в декартовых координатах
/// \details Вычисляется как d = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
/// \param[in] point1 - 1 точка
/// \param[in] point2 - 2 точка
/// \return Расстояние между двумя точками в декартовых координатах
///
double XYZtoDistance( const XYZ &point1, const XYZ &point2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief ECEF смещение ( разница в декартовых ECEF координатах двух точек )
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] lat1 - широта 1 точки
/// \param[in] lon1 - долгота 1 точки
/// \param[in] h1 - высота 1 точки
/// \param[in] lat2 - широта 2 точки
/// \param[in] lon2 - долгота 2 точки
/// \param[in] h2 - высота 2 точки
/// \param[out] dX - смещение по оси X
/// \param[out] dY - смещение по оси Y
/// \param[out] dZ - смещение по оси Z
///
void ECEF_offset( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat1, double lon1, double h1, double lat2, double lon2, double h2, double &dX, double &dY, double &dZ );

///
/// \brief ECEF смещение ( разница в декартовых ECEF координатах двух точек )
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point1 - точка 1
/// \param[in] point2 - точка 2
/// \return ECEF смещение
///
XYZ ECEF_offset( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic &point1, const Geodetic &point2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод ECEF координат точки в ENU относительно географических координат (lat, lon)
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] dX - смещение по оси X
/// \param[in] dY - смещение по оси Y
/// \param[in] dZ - смещение по оси Z
/// \param[in] lat - широта точки
/// \param[in] lon - долгота точки
/// \param[out] xEast - ENU координата X (East)
/// \param[out] yNorth - ENU координата Y (North)
/// \param[out] zUp - ENU координата X (Up)
///
void ECEFtoENU( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double dX, double dY, double dZ, double lat, double lon, double &xEast, double &yNorth, double &zUp );

///
/// \brief Перевод ECEF координат точки в ENU относительно географических координат point
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] shift - смещение по декартовым осям
/// \param[in] point - точка
/// \return Координаты ENU точки point
///
ENU ECEFtoENU( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const XYZ &shift, const Geographic &point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод ENU коордиат точки в AER координаты
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] xEast - ENU координата X (East)
/// \param[in] yNorth - ENU координата Y (North)
/// \param[in] zUp - ENU координата X (Up)
/// \param[out] az - азимут
/// \param[out] elev - угол места
/// \param[out] slantRange - наклонная дальность
///
void ENUtoAER( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double xEast, double yNorth, double zUp, double &az, double &elev, double &slantRange );

///
/// \brief Перевод ENU коордиат точки в AER координаты
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point - точка в координатах ENU
/// \return Координаты AER точки point
///
AER ENUtoAER( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit, const ENU &point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод геодезических координат GEO точки point в координаты ENU относительно опорной точки
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] lat - широта точки
/// \param[in] lon - долгота точки
/// \param[in] h - высота точки
/// \param[in] lat0 - широта опорной точки
/// \param[in] lon0 - долгота опорной точки
/// \param[in] h0 - высота опорной точки
/// \param[out] xEast - ENU координата X (East)
/// \param[out] yNorth - ENU координата Y (North)
/// \param[out] zUp - ENU координата X (Up)
///
void GEOtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat, double lon, double h, double lat0, double lon0, double h0, double &xEast, double &yNorth, double &zUp );

///
/// \brief Перевод геодезических координат GEO точки point в координаты ENU относительно опорной точки
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point - геодезические координаты точки
/// \param[in] anchor - геодезические координаты опорной точки, относительно которой переводим
/// \return Координаты ENU точки point
///
ENU GEOtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic &point, const Geodetic &anchor );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вычисление AER координат между двумя геодезическими точками
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] lat1 - широта 1 точки
/// \param[in] lon1 - долгота 1 точки
/// \param[in] h1 - высота 1 точки
/// \param[in] lat2 - широта 2 точки
/// \param[in] lon2 - долгота 2 точки
/// \param[in] h2 - высота 2 точки
/// \param[out] az - азимут из 1 точки на 2 точку
/// \param[out] elev - угол места из 1 точки на 2 точку
/// \param[out] slantRange - наклонная дальность
///
void GEOtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat1, double lon1, double h1, double lat2, double lon2, double h2, double &az, double &elev, double &slantRange );

///
/// \brief Вычисление AER координат между двумя геодезическими точками
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point1 - геодезические координаты 1 точки
/// \param[in] point2 - геодезические координаты 2 точки
/// \return Координаты AER между точками point1 и point2
///
AER GEOtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     const Geodetic &point1, const Geodetic &point2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Косинус угла между векторами в евклидовом пространстве
/// \details Предполагается, что оба вектора начинаются в точке (0, 0, 0)
/// \param[in] x1 - X координата 1 точки
/// \param[in] y1 - Y координата 1 точки
/// \param[in] z1 - Z координата 1 точки
/// \param[in] x2 - X координата 2 точки
/// \param[in] y2 - Y координата 2 точки
/// \param[in] z2 - Z координата 2 точки
/// \return Косинус угла между векторами, [радиан]
///
double CosAngleBetweenVectors( double x1, double y1, double z1, double x2, double y2, double z2 );

///
/// \brief Косинус угла между векторами в евклидовом пространстве
/// \details Предполагается, что оба вектора начинаются в точке (0, 0, 0)
/// \param[in] point1 - 1 точка
/// \param[in] point2 - 2 точка
/// \return Косинус угла между векторами, [радиан]
///
double CosAngleBetweenVectors( const XYZ &point1, const XYZ &point2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Угол между векторами в евклидовом пространстве
/// \details Предполагается, что оба вектора начинаются в точке (0, 0, 0)
/// \param[in] x1 - X координата 1 точки
/// \param[in] y1 - Y координата 1 точки
/// \param[in] z1 - Z координата 1 точки
/// \param[in] x2 - X координата 2 точки
/// \param[in] y2 - Y координата 2 точки
/// \param[in] z2 - Z координата 2 точки
/// \return Угол между векторами, [радиан]
///
double AngleBetweenVectors( double x1, double y1, double z1, double x2, double y2, double z2 );

///
/// \brief Угол между векторами в евклидовом пространстве
/// \details Предполагается, что оба вектора начинаются в точке (0, 0, 0)
/// \param[in] vec1 - вектор 1
/// \param[in] vec2 - вектор 2
/// \return Угол между векторами, [радиан]
///
double AngleBetweenVectors( const XYZ &vec1, const XYZ &vec2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вектор из координат двух точек
/// \details Предполагается, что результирующий вектор начинается в точке (0, 0, 0)
/// \param[in] x1 - X координата 1 точки
/// \param[in] y1 - Y координата 1 точки
/// \param[in] z1 - Z координата 1 точки
/// \param[in] x2 - X координата 2 точки
/// \param[in] y2 - Y координата 2 точки
/// \param[in] z2 - Z координата 2 точки
/// \param[out] xV - X координата вектора с началом в точке 1 и концом в точке 2
/// \param[out] yV - Y координата вектора с началом в точке 1 и концом в точке 2
/// \param[out] zV - Z координата вектора с началом в точке 1 и концом в точке 2
///
void VectorFromTwoPoints( double x1, double y1, double z1, double x2, double y2, double z2, double &xV, double &yV, double &zV );

///
/// \brief Вектор, полученный из координат двух точек
/// \details Предполагается, что результирующий вектор начинается в точке (0, 0, 0)
/// \param[in] point1 - 1 точка
/// \param[in] point2 - 2 точка
/// \return Вектор, полученный из координат двух точек
///
XYZ VectorFromTwoPoints( const XYZ &point1, const XYZ &point2 );

} // end namespace SPML
} // end namespace Geodesy
#endif // SPML_GEODESY_H
/// \}
