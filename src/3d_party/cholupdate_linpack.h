//----------------------------------------------------------------------------------------------------------------------
///
/// \file       cholupdate_linpack.h
/// \brief      Процедура cholupdate (осуществляет вызов linpack процедур: dchud для update и dchdd для downdate)
/// \date       17.02.22 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef CHOLUPDATE_LINPACK_H
#define CHOLUPDATE_LINPACK_H

#include <armadillo>
#include <cassert>

#include <blas0.h>
#include <blas1_d.h>
#include <linpack_d.h>

///
/// \brief Процедура cholupdate
/// \details Осуществляет вызов процедур dchud для update и dchdd для downdate;
/// процедуры dchud и dchdd взяты из пакета linpack ()
/// \attention Вычисляет A~ = A + ( sqrt(v) * X ) * ( sqrt(v) * X )'
/// \attention В случае, если фактора масштабирования нет, надо задать v = 1 или v = -1 (для update и downdate соответственно)
/// \param[in out]  A - матрица входа/выхода
/// \param[in]      X - вектор
/// \param[in]      v - фактор масштабирования/управления: ( v > 0 - это update, v < 0 - это downdate )
/// \return 0 - успех, 1 - неудача (в случае downdate)
///
int cholupdate_linpack( arma::mat &A, const arma::vec &X, double v );

#endif // CHOLUPDATE_LINPACK_H
/// \}
