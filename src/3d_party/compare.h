//----------------------------------------------------------------------------------------------------------------------
///
/// \file       compare.h
/// \brief      Функции сравнения чисел, массивов
/// \date       27.07.20 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_COMPARE_H
#define SPML_COMPARE_H

// System includes:
#include <cmath>

// SPML includes:
#include <consts.h>

namespace SPML /// Специальная библиотека программных модулей (СБПМ)
{
namespace Compare /// Сравнение чисел
{
//----------------------------------------------------------------------------------------------------------------------
static const double EPS = 1.0e-8; ///< Абсолютная точность по умолчанию при сравнениях чисел типа double (1.0e-8)
static const float EPS_REL = 0.01; ///< Относительная точность по умолчанию

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Сравнение двух действительных чисел (по абсолютной разнице)
/// \details Возвращает результат: abs( first - second ) < eps
/// \param[in] first - первое число
/// \param[in] second - второе число
/// \param[in] eps - абсолютная точность сравнения
/// \return true - если разница меньше точности, иначе false
///
inline bool AreEqualAbs( double first, double second, const double &eps = EPS )
{
    return ( std::abs( first - second ) <= eps );
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Сравнение двух действительных чисел (по относительной разнице)
/// \details Возарвщает результат: ( abs( ( first - second ) / first ) < eps ) && ( abs( ( first - second ) / second ) < eps )
/// \param[in] first - первое число
/// \param[in] second - второе число
/// \param[in] eps - относительная точность сравнения
/// \return true - если разница меньше точности, иначе false
///
inline bool AreEqualRel( double first, double second, const double &eps = EPS_REL )
{
    return ( ( std::abs( first - second ) <= ( eps * std::abs( first ) ) ) &&
        ( std::abs( first - second ) <= ( eps * std::abs( second ) ) ) );
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Проверка действительного числа на равенство нулю  (по абсолютной разнице)
/// \details Возвращает результат: abs( value ) < eps
/// \param[in] value - проверяемое число
/// \param[in] eps - абсолютная точность сравнения
/// \return true - если разница меньше точности, иначе false
///
inline bool IsZeroAbs( double value, const double &eps = EPS )
{
    return ( std::abs( value ) <= eps );
}

} // end namespace Compare
} // end namespace SPML
#endif // SPML_COMPARE_H
/// \}
