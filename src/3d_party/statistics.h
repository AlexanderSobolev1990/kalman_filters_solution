//----------------------------------------------------------------------------------------------------------------------
///
/// \file       statistics.h
/// \brief      Различные функции математической статистики
/// \date       08.02.21 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_STATISTICS_H
#define SPML_STATISTICS_H

// System includes:
#include <cmath>
#include <cassert>
#include <algorithm>
#include <armadillo>
#include <boost/math/distributions.hpp>

// SPML includes:
#include <compare.h>

namespace SPML /// Специальная библиотека программных модулей (СБПМ)
{
namespace Statistics /// Статистика
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вычисление факториала
/// \details Рекурентный способ, ограничение на аргумент x <= 12, иначе результат не поместится в unsigned unt (4 байта)
/// \param x - аргумент
/// \return Факториал аргумента x
///
unsigned int Factorial( unsigned int x );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Число сочетаний из n по k (без повторений), также известно как биномиальный коэффициент
/// \details Рекурентный способ
/// \param k - набор из k элементов
/// \param n - n-элементарное множество
/// \return Число сочетаний из n по k
///
unsigned int C_k_n( unsigned int k, unsigned int n );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Проверка на четность целого числа
/// \param value - проверяемое целое число
/// \return true - если число четное, false - если нечетное
///
inline bool IsEven( int value );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вычисление объема n-мерной сферы радиуса r
/// \details Вычисляется по формуле для четных/нечетных n, нет ограничения на n
/// \param n - количество измерений пространства
/// \param r - радиус сферы (по умолчанию r = 1, единичная сфера)
/// \return Объем n-мерной сферы радиуса r
///
double VolumeOfSphere_formula( unsigned int n, double r = 1.0 );

///
/// \brief Вычисление объема n-мерной сферы радиуса r
/// \details Заданы предопределенные значения при n=0..10
/// \param n - количество измерений пространства
/// \param r - радиус сферы (по умолчанию r = 1, единичная сфера)
/// \return Объем n-мерной сферы радиуса r
///
double VolumeOfSphere_hardcoded( unsigned int n, double r = 1.0 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вычисление относительного объема строба (RVSTR - relative volume of strobe)
/// \param[in] det - определитель матрицы S
/// \param[in] strobe - строб
/// \param[in] volume - объем, относительно которого вычисляется RVSTR (может быть объем зоны обзора, эл-та разрешения и т.д.)
/// \param[in] sizeY - размерность вектора измерений
/// \return Относительный объем строба
///
double CalcRVSTR( double det, double strobe, double volume, int sizeY );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Гауссовская (нормальная) функция распределения (CDF - cumulative distribution function)
/// \details Возвращает значение ( 1.0 + std::erf( ( x - m ) / ( s * std::sqrt( 2.0 ) ) ) ) / 2.0
/// \param[in] arg - аргумент
/// \param[in] m - математическое ожидание
/// \param[in] s - среднеквадратическое отклонение (СКО)
/// \return Значение функции распределения при заданных аргументах
///
inline double GaussCDF( double arg, double m, double s )
{
    double gaussCDF = ( 1.0 + std::erf( ( arg - m ) / ( s * std::sqrt( 2.0 ) ) ) ) / 2.0;
    return gaussCDF;
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Функция, обратная функции erf
/// \details https://github.com/lakshayg/erfinv
/// \param x - аргумент
/// \return Возвращает erf^-1( x )
///
double ErfInv( double x );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Нахождение аргумента по значению Гауссовской (нормальной) функции распределения (CDF - cumulative distribution function)
/// \param gaussCDF - значение Гауссовской (нормальной) функции распределения
/// \param m - математическое ожидание
/// \param s - среднеквадратическое отклонение (СКО)
/// \return Возвращает аргумент, соответствующий значению Гауссовской (нормальной) функции распределения
///
inline double ArgGaussCDF( double gaussCDF, double m, double s )
{
    double arg = ( std::sqrt( 2.0 ) * s * ErfInv( ( 2.0 * gaussCDF ) - 1.0 ) ) + m;
    return arg;
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Аппроксимация Пейзера-Пратта для функции распределения хи-квадрат (CDF - cumulative distribution function)
/// \details По материалам "Аппроксимация центрального распределения хи-квадрат для оперативного расчета
/// вероятности ложной тревоги энергетического обнаражителя" Лесников, Наумович, Частиков, Дубовцев,
/// Вятский государственный университет, Киров, МЭС-2016
/// \param[in] arg - аргумент
/// \param[in] degree - степень свободы
/// \return Возвращает F = ChiSquared(arg^2, degree)
///
double ChiSquaredCDF_PeizerPratt( double arg, unsigned int degree );

///
/// \brief Функция распределения Хи-квадрат (CDF - cumulative distribution function) для четных степеней свободы (точная, не приближенная формула)
/// \details По материалам "Аппроксимация центрального распределения хи-квадрат для оперативного расчета
///          вероятности ложной тревоги энергетического обнаражителя" Лесников, Наумович, Частиков, Дубовцев,
///          Вятский государственный университет, Киров, МЭС-2016
/// \param[in] arg - аргумент
/// \param[in] degree - степень свободы
/// \return Возвращает F = ChiSquared(x^2, degree)
///
double ChiSquaredCDF_EvenDegree( double arg, unsigned int degree );
/*
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Функция распределения F( x^2 ) = ( 1 - CDF( x^2 ) ), где CDF - cumulative distribution function.
/// \param[in] arg    - аргумент
/// \param[in] degree - степень свободы
/// \return Возвращает значение F( x^2 ) = 1 - ChiSquared( x^2, degree )
///
double Freqs( double arg, int degree );
double FreqsNew( double arg, int degree );
double FreqsNew2( double arg, int degree );
*/
///
/// \brief Метод вычисления значения функции распределения Хи-квадрат
///
enum TChiSquaredMethod : int
{
    CSM_Boost = 0, ///< Обращение к библиотеке Boost
    CSM_EvenDegree = 1, ///< Вычисление прямой формулой для четных значений степени свободы (если нечетная, то аппроксимацией CSM_Approx)
    CSM_Approx = 2 ///< Вычисление аппроксимацией Пейзера-Пратта Гауссовской функцией распределения
};

///
/// \brief Функция распределения Хи-квадрат (CDF - cumulative distribution function)
/// \param[in] arg - аргумент
/// \param[in] degree - число степеней свободы
/// \param[in] method - метод вычисления функции хи-квадрат
/// \return Возвращает значение F( x^2 ) = 1 - ChiSquared( x^2, degree )
///
double ChiSquaredCDF( double arg, unsigned int degree, const TChiSquaredMethod method = TChiSquaredMethod::CSM_EvenDegree );

///
/// \brief Функция распределения "1 - Хи-квадрат": F( x^2 ) = ( 1 - ChiSquaredCDF( x^2 ) ), где CDF - cumulative distribution function.
/// \param[in] arg - аргумент
/// \param[in] degree - число степеней свободы
/// \param[in] method - метод вычисления функции хи-квадрат
/// \return Возвращает значение F( x^2 ) = 1 - ChiSquared( x^2, degree )
///
double InvChiSquaredCDF( double arg, unsigned int degree, const TChiSquaredMethod method = TChiSquaredMethod::CSM_EvenDegree );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Экспоненциальное скользящее среднее (Exponential Moving Average)
/// \param[in] ema_prev - предыдущее значение EMA (в момент времени t-1)
/// \param[in] alpha - сглаживающая переменная (0..1)
/// \param[in] x - сглаживаемое значение (в текущий момент времени t)
/// \return Экспоненциально усредненное значение аргумента x
///
inline double EMA( double ema_prev, double alpha, double x )
{
    assert( alpha > 0.0 && alpha < 1.0 );
    return ( ( alpha * x ) + ( ( 1.0 - alpha ) * ema_prev ) );
}

///
/// \brief Экспоненциальное скользящее среднее (Exponential Moving Average)
/// \param[in] ema_prev - предыдущее значение EMA (в момент времени t-1)
/// \param[in] alpha - сглаживающая переменная (0..1)
/// \param[in] x - сглаживаемое значение (в текущий момент времени t)
/// \return Экспоненциально усредненное значение аргумента x
///
inline double EMA( double ema_prev, double alpha, int x )
{
    assert( alpha > 0.0 && alpha < 1.0 );
    return ( ( alpha * x ) + ( ( 1.0 - alpha ) * ema_prev ) );
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Критерий Вальда
/// \details Модифицирован по источнику: Абчук В.А. Справочник по исследованию операций/ Под общ. ред. Ф.А.Матвейчука - М.:Воентиздат, 1979, стр 29
/// \details и первоисточнику: Вальд А. Последовательный анализ под ред. А.Ф. Лапко М.:Гос. изд-во физ.-мат. лит. 1960, стр 125
/// \param[in] pl - вероятность ложной отметки в стробе
/// \param[out] pc - вероятность целевой отметки в стробе
/// \param[in] d - заданная вероятность правильного обнаружения
/// \param[out] dcorr - скорректированная вероятность правильного обнаружения
/// \param[in] f - поток ложных тревог, приведенный к такту измерений
/// \param[in] ktsol_max - максимальное число тактов на принятие решений (вычисляется из времени обнаружения и темпа поступления отметок)
/// \param[out] ktsol_corr - скорректированное число тактов на принятие решений
/// \param[out] ad - коэффициент порога обнаружения
/// \param[out] ar - коэффициент порога сброса
/// \param[out] b - коэффициент b
/// \return 0 в случае успешных вычислений, 1 в случае неудачи
///
bool WaldCriterion( double pl, double &pc, double d, double &dcorr, double f, double ktsol_max, double &ktsol_corr,
    double &ad, double &ar, double &b );

int KfromNcriterion( double f, int n, double pl );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вычислить расстояние Махаланобиса для матрицы, обратной ковариационной матрице S, и вектора невязки V
/// \param V - вектор невязки
/// \param Sinv - матрица, обратная ковариационной матрице вектора невязки V
/// \return расстояние Махаланобиса md = sqrt( abs( V.transposed * Sinv * V ) )
///
inline double CalculateMahalanobisDistance( const arma::vec &V, const arma::mat &Sinv )
{
    arma::mat m11 = ( V.t() * Sinv * V );
    double mahalanobisDistance = std::sqrt( std::abs( m11(0,0) ) );
    return mahalanobisDistance;
}


// TODO: Решить надо ли это?
inline double JaccardIndex( double p1, double p2 )
{
    return ( ( p1 * p2 ) / ( p1 + p2 - ( p1 * p2 ) ) );
}

inline double GetPfromW( double w1, double w2 )
{
    double w1_norm = w1 / std::max( w1, w2 );
    double w2_norm = w2 / std::max( w1, w2 );
    double p = JaccardIndex( w1_norm, w2_norm );
    return p;
}

inline double ProbabilityDensityFunction( const arma::vec &Z, const arma::mat &R, const arma::vec &X )
{
    arma::vec delta = ( Z - X );
    arma::mat m11 = ( arma::trans( delta ) * arma::inv( R ) * delta );
    return ( ( 1.0 / ( ( 2.0 * Consts::PI ) * std::sqrt( arma::det( R ) ) ) ) * std::exp( -0.5 * m11(0,0) ) );
}

inline double ProbabilityDensityFunction( const arma::vec &Z, const arma::mat &R )
{
    arma::mat m11 = ( arma::trans( Z ) * arma::inv( R ) * Z );
    return ( ( 1.0 / ( ( 2.0 * Consts::PI ) * std::sqrt( arma::det( R ) ) ) ) * std::exp( -0.5 * m11(0,0) ) );
}

} // end Statistics
} // end SPML
#endif // SPML_STATISTICS_H
/// \}
