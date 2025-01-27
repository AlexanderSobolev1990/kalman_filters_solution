//----------------------------------------------------------------------------------------------------------------------
///
/// \file       numericalmethods.h
/// \brief      Численные методы
/// \date       21.10.21 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_NUMERICALMETHODS_H
#define SPML_NUMERICALMETHODS_H

// System includes:
#include <cassert>
#include <functional>

// SPML includes:
#include <compare.h>

namespace SPML /// Специальная библиотека программных модулей (СБПМ)
{
namespace NumericalMethods /// Численные методы
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Формула метода секущих
/// \param[in] x1 - текущее значение аргумента x(k)
/// \param[in] x0 - значение аргумента на прошлом шаге x(k-1)
/// \param[in] f_x1 - текущее значение функции f(x(k))
/// \param[in] f_x0 - значение функции на прошлом шаге f(x(k-1))
/// \return Следующее значение аргумента: x(k+1)
///
inline double SecantMethodFormula( double x1, double x0, double f_x1, double f_x0 )
{
    double x = x1 - ( ( f_x1 / ( f_x1 - f_x0 ) ) * ( x1 - x0 ) );
    return x;
}

//template <typename Func, typename... Args>
//int SomeMethod( Func&& func, Args&&... args ) //
//{
//    int res = func(args...) + 777;
//    return res;
//}
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Метод секущих
/// \param[in] f - функция, корень которой ищем
/// \param[in] x_start - начальное приближение корня
/// \param[in] step - шаг поиска корня (может быть как положительным, как и отрицательным)
/// \param[in] tol - абсолютная точность нахождения корня
/// \param[in] iterMax - макисмальное число итераций
/// \param[in] plus - true: искать положительный корень
/// \return Возвращает пару значений: корень x0 и значение f(x0)
///
template <class F>
inline std::pair<double, double> SecantMethod( F f, double x_start, double step, double tol, int iterMax, bool plus )
{
    double x0 = x_start;
    double x1 = x_start + step;
    double f_x0 = f( x0 );
    double f_x1 = f( x1 );

    int iter = 0;
    double x2 = 0.0;

    while( ( ( std::abs( x1 - x0 ) > tol ) || !( ( std::abs( x1 ) > tol ) && plus ) ) && ( iter < iterMax ) )
    {
        x2 = SecantMethodFormula( x1, x0, f_x1, f_x0 );
        x0 = x1;
        x1 = x2;
        f_x0 = f( x0 );
        f_x1 = f( x1 );
        iter++;
    }
    return std::make_pair( x1, f_x1 );
}

} // end namespace AlgebraicEquations
} // end namespace SPML
#endif // SPML_NUMERICALMETHODS_H
/// \}
