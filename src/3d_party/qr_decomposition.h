//----------------------------------------------------------------------------------------------------------------------
///
/// \file       qr_decomposition.h
/// \brief      QR-разложение матриц
/// \date       01.04.21 - создан
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
namespace QR /// QR-разложение матриц
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Функция знака числа sgn
/// \return 1 : число положительное, -1 : число отрицательное, 0 : ноль
///
template <typename T>
int sgn( T val ) {
    return ( T(0) < val ) - ( val < T(0) );
}
//----------------------------------------------------------------------------------------------------------------------
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

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief J-ортогональное QR разложение
/// \details Приводится по:
/// \param[out] Q    - выходная Q матрица
/// \param[out] R    - выходная R матрица
/// \param[out] J1   - выходная матрица знаков J
/// \param[in]  A    - входная матрица
/// \param[in]  J0   - входная матрица знаков J
/// \param[in]  econ - признак экономичного возврата (true по-умолчанию)
/// \param[in]  tol  - точность при сравнениях double
/// \return 0 - успех, 1 - ошибка
///
int J_orthogonal( arma::mat &Q, arma::mat &R, arma::mat &J1, const arma::mat &A, const arma::mat &J0,
    bool econ = true, double tol = 1.0e-16 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Модифицированное QR-разложение Грама-Шмидта
/// \details Приводится по G. W. Stewart, “Matrix Algorithms, Volume 1: Basic Decompositions”, SIAM, 1998, стр.300(279)
/// \details Экономичный возврат
/// \param[out] Q    - выходная Q матрица
/// \param[out] R    - выходная R матрица
/// \param[in]  A    - входная матрица
/// \return 0 - успех, 1 - ошибка
///
int ModifiedGramSchmidt( arma::mat &Q, arma::mat &R, const arma::mat &A );

int ModifiedGramSchmidtRowByRow( arma::mat &Q, arma::mat &R, const arma::mat &A );

//----------------------------------------------------------------------------------------------------------------------
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
int SchwarzRutishauser( arma::mat &Q, arma::mat &R, const arma::mat &A, bool econ = true, double tol = 1.0e-9 );

} // end QR
} // end SPML
#endif // SPML_QR_H
/// \}
