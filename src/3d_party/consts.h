//----------------------------------------------------------------------------------------------------------------------
///
/// \file       consts.h
/// \brief      Константы библиотеки СБПМ
/// \date       27.07.20 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_CONSTS_H
#define SPML_CONSTS_H

// System includes:
#include <cmath>

namespace SPML /// Специальная библиотека программных модулей (СБПМ)
{
namespace Consts /// Константы
{
//
// Скорость света
const double C = 3.0e8; ///< Скорость света, [м/с] в двойной точности (double)
//
// Число ПИ и его части
const double PI = std::acos( -1.0 ); ///< Число PI = 3.14... в радианах в двойной точности (double)
const double PI_2 = 2.0 * PI; ///< Число 2*PI = 6.28... в радианах в двойной точности (double)
const double PI_05 = 0.5 * PI; ///< Число PI/2 = 1.57... в радианах в двойной точности (double)
const double PI_025 = 0.25 * PI; ///< Число PI/4 = 0.785... в радианах в двойной точности (double)

} // end namespace Consts
} // end namespace SPML
#endif // SPML_CONSTS_H
/// \}
