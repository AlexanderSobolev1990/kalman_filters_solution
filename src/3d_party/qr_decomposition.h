//------------------------------------------------------------------------------
///
/// \file       qr_decomposition.h
/// \brief      QR-разложение матриц
/// \date       01.04.21 - создан, 13.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_QR_H
#define SPML_QR_H

// System includes:
#include <armadillo>
#include <deque>
#include <cassert>
#include <float.h> // для DBL_EPSILON

namespace SPML /// Специальная библиотека программных модулей (СБ ПМ)
{

//------------------------------------------------------------------------------
///
/// \brief Функция знака числа sgn
/// \return 1 : число положительное, -1 : число отрицательное, 0 : ноль
///
template <typename T>
int sgn( T val ) {
    return ( T(0) < val ) - ( val < T(0) );
}

///
/// \brief Функция знака числа sgn2
/// \return 1 : число положительное или ноль, -1 : число отрицательное
///
template <typename T>
int sgn2( T val, T tol ) {
    if( std::abs( val ) < DBL_EPSILON ) {
        return 1;
    } else {
        return sgn( val );
    }
}

arma::vec sgn2( arma::vec v );

//------------------------------------------------------------------------------    
namespace QR /// QR-разложение матриц
{
///
/// \brief Модифицированное QR-разложение Грама-Шмидта
/// \details Приводится по G. W. Stewart, “Matrix Algorithms, Volume 1: 
///          Basic Decompositions”, SIAM, 1998, стр.300(279)
/// \param[out] Q    - выходная Q матрица
/// \param[out] R    - выходная R матрица
/// \param[in]  A    - входная матрица
///
void MGS( arma::mat &Q, arma::mat &R, const arma::mat &A );

//------------------------------------------------------------------------------
///
/// \brief Модифицированное QR-разложение Грама-Шмидта,
///        возвращающее только матрицу R
/// \details Приводится по G. W. Stewart, “Matrix Algorithms, Volume 1: 
///          Basic Decompositions”, SIAM, 1998, стр.300(279)
/// \attention Возвращает только матрицу R
/// \param[out] R    - выходная R матрица
/// \param[in]  A    - входная матрица
///
void MGS_1( arma::mat &R, const arma::mat &A );

//------------------------------------------------------------------------------
///
/// \brief Модифицированное QR-разложение Шварца-Рутисхаузера (Schwarz-Rutishauser)
/// \details Приводится по Walter Gander, Algorithms for the QR-Decomposition, RESEARCH REPORT NO. 80-02, APRIL 1980,
/// SEMINAR FUER ANGEWANDTE MATHEMATIK EIDGENOESSISCHE TECHNISCHE HOCHSCHULE CH-8092 ZUERICH
/// \details Модифицировано под C++:
/// https://github.com/SPancratz/flint2/blob/69395ab1a939212533ac9cdbff69f67e0b35da8f/mpf_mat/qr.c
/// \param[out] Q    - выходная Q матрица
/// \param[out] R    - выходная R матрица
/// \param[in]  A    - входная матрица
/// \param[in]  econ - признак экономичного возврата (true по-умолчанию)
/// \param[in]  tol  - точность при сравнениях double
/// \return 0 - успех, 1 - ошибка
///
int SchwarzRutishauser( arma::mat &Q, arma::mat &R, const arma::mat &A, 
    bool econ = true, double tol = 1.0e-9 );

} // end QR

//------------------------------------------------------------------------------
namespace LQ /// LQ-разложение матриц
{
///
/// \brief Модифицированное LQ-разложение Грама-Шмидта
/// \details Приводится по G. W. Stewart, “Matrix Algorithms, Volume 1: 
///          Basic Decompositions”, SIAM, 1998, стр.300(279)
/// \details Экономичный возврат
/// \param[out] Q    - выходная Q матрица
/// \param[out] L    - выходная L матрица
/// \param[in]  A    - входная матрица
///
void MGS( arma::mat &Q, arma::mat &L, const arma::mat &A );

//------------------------------------------------------------------------------
///
/// \brief Модифицированное LQ-разложение Грама-Шмидта, 
///        возвращающее только матрицу L
/// \details Приводится по G. W. Stewart, “Matrix Algorithms, Volume 1: 
///          Basic Decompositions”, SIAM, 1998, стр.300(279)
/// \details Экономичный возврат
/// \param[out] L    - выходная L матрица
/// \param[in]  A    - входная матрица
/// \return 0 - успех, 1 - ошибка
///
void MGS_1( arma::mat &L, const arma::mat &A );

}

//------------------------------------------------------------------------------
namespace JQR /// JQR-разложение матриц
{
///
/// \brief J-ортогональное QR разложение
/// \details Приводится по:
/// \param[out] Q    - выходная Q матрица
/// \param[out] R    - выходная R матрица
/// \param[out] Jp   - выходная матрица знаков J
/// \param[in]  A    - входная матрица
/// \param[in]  J    - входная матрица знаков J
/// \param[in]  econ - признак экономичного возврата (true по-умолчанию)
/// \param[in]  tol  - точность при сравнениях double
/// \return 0 - успех, 1 - ошибка
///
int JQR( arma::mat &Q, arma::mat &R, arma::vec &Jp, const arma::mat &A, 
    const arma::vec &J, bool econ = true, double tol = 1.0e-16 );

//------------------------------------------------------------------------------
///
/// \brief J-ортогональное QR разложение, возвращающее только матрицу R
/// \details Приводится по:
/// \attention Возвращает только матрицу R
/// \param[out] R   - выходная R матрица
/// \param[in]  A   - входная матрица
/// \param[in]  J   - входная матрица знаков J
/// \return 0 - успех, 1 - ошибка
///
int JQR_1( arma::mat &R, const arma::mat &A, const arma::vec &J );    

}

//------------------------------------------------------------------------------
namespace JLQ /// JLQ-разложение матриц
{
///
/// \brief J-ортогональное QR разложение
/// \details Приводится по:
/// \param[out] Q    - выходная Q матрица
/// \param[out] L    - выходная L матрица
/// \param[out] Jp   - выходная матрица знаков J
/// \param[in]  A    - входная матрица
/// \param[in]  J    - входная матрица знаков J
/// \param[in]  econ - признак экономичного возврата (true по-умолчанию)
/// \param[in]  tol  - точность при сравнениях double
/// \return 0 - успех, 1 - ошибка
///
int JLQ( arma::mat &Q, arma::mat &L, arma::vec &Jp, const arma::mat &A, 
    const arma::vec &J, bool econ = true, double tol = 1.0e-16 );

//------------------------------------------------------------------------------
///
/// \brief J-ортогональное QR разложение, возвращающее только матрицу L
/// \details Приводится по:
/// \attention Возвращает только матрицу L
/// \param[out] L   - выходная L матрица
/// \param[in]  A   - входная матрица
/// \param[in]  J   - входная матрица знаков J
/// \return 0 - успех, 1 - ошибка
///
int JLQ_1( arma::mat &L, const arma::mat &A, const arma::vec &J );    

}

} // end SPML
#endif // SPML_QR_H
/// \}
