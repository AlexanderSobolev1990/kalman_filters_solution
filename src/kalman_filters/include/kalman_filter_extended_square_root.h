//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_extended_square_root.h
/// \brief      Шаблонный класс квадратно-корневого  расширенного фильтра Калмана, КК-РФК
///             (Square Root Extended Kalman Filter, SREKF)
/// \date       11.01.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_SQUARE_ROOT_H
#define KALMAN_FILTER_EXTENDED_SQUARE_ROOT_H

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_extended.h>
#include <qr_decomposition.h> // JQR разложение матрицы
#include <cholupdate_linpack.h>

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Класс квадратно-корневого расширенного фильтра Калмана, КК-РФК (Square Root Extended Kalman Filter, EKF)
/// \details Источник: NASA Technical report R-135, Application of statistical filter theory
/// to the optimal estimation of position and velocity on board a circumlunar vehicle,
/// Gerald L. Smith, Stanley F. Schmidt and Leonard A. McGee, 1962
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSREKF : virtual public CKalmanEKF<SizeX, SizeY>
{
public:
    using CKalmanEKF<SizeX, SizeY>::CKalmanEKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSREKF() : CKalmanEKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SREKF" );
#endif
    }
    // default copy/move/assignment semantic:
//    CKalmanSREKF( const CKalmanSREKF& ) = default;
//    CKalmanSREKF& operator=( const CKalmanSREKF& ) = default;
//    CKalmanSREKF( CKalmanSREKF&& ) = default;
//    CKalmanSREKF& operator=( CKalmanSREKF&& ) = default;
    virtual ~CKalmanSREKF() = default;

    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
        PredictionSREKF( dt );
    }

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Correction started" << std::endl;
        Y_msd.print( this->filterName_ + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла:
        this->prediction_isDone = false; // Сразу же снять признак

        // 1. Вычисление ковариационной матрицы S        
        arma::mat S_tmp = this->H_ * this->P_;
        arma::mat QR_input_S_corr = arma::trans( arma::join_horiz( S_tmp, this->R_ ) ); // QR_input = [ S_, R ], Транспонировать, т.к. в QR разложение так надо
        arma::mat Q_qr_S_corr, R_qr_S_corr;
        int res_QR_S_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_S_corr, R_qr_S_corr, QR_input_S_corr ); // Вызов QR разложения
        assert( res_QR_S_corr == 0 );
        this->S_ = arma::trans( R_qr_S_corr ); // Матрица S считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( this->filterName_ + " Correction, S:" );
#endif
        this->checkMatrixDiagPositive( this->S_ );

        // 2. Вычисление коэффициента усиления фильтра K
        arma::mat P_xy = this->P_ * arma::trans( this->P_ ) * arma::trans( this->H_ ); // Матрица кросс-коварации векторов Х и Y, размерность [SizeX * SizeY]
        this->K_ = ( P_xy * arma::inv( arma::trans( this->S_ ) ) ) * arma::inv( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->K_ ).print( this->filterName_ + " Correction, K:" );
#endif
        // 3. Вычисление невязки Delta
        if( !this->deltaY_isSet ) {
            assert( this->Y_msd_isSet ); // Если не установлен deltaY, то Y_msd_ обязан быть установлен
            this->Y_msd_isSet = false; // Сразу же снять признак
            this->DeltaY_ = this->Y_msd_ - this->Y_pred_;
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->DeltaY_ = this->checkDeltaMeasurement_( this->DeltaY_ );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->DeltaY_ ).print( this->filterName_ + " Correction, DeltaY:" );
#endif
        // 4. Вычисление ковариационной матрицы P
        arma::mat A = this->K_ * this->S_;
        for( unsigned i = 0; i < A.n_cols; i++ ) {
            cholupdate_linpack( this->P_, A.col(i), -1.0 );
        }
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Correction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );
        // 5. Вычисление X_est, Y_est
        this->X_est_ = this->X_pred_ + ( this->K_ * this->DeltaY_ );
        if( this->checkBordersStateAfterCorrection_ != nullptr ) {
            this->X_est_ = this->checkBordersStateAfterCorrection_( this->X_est_ );
        }
        this->Y_est_ = this->observationModel_( this->X_est_ );
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_est_ = this->checkBordersMeasurement_( this->Y_est_ );
        }
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( this->filterName_ + " Correction, X_est:" );
        ( this->Y_est_ ).print( this->filterName_ + " Correction, Y_est:" );
#endif
        this->deltaY_isSet = false; // Снять признак (выставляется в true в сеттере deltaY)
    }

protected:
    ///
    /// \brief Прогноз SREKF
    /// \param dt - Время прогноза, [с]
    ///
    void PredictionSREKF( double dt )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Prediction started, dt = " << dt << std::endl;
#endif
        assert( dt > 0.0 ); // Прогноз имеет смысл только при положительном времени

        // 1. Вычисление X_pred, F
        this->X_pred_ = this->stateTransitionModel_( this->X_est_, dt );
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }
        this->F_ = this->stateTransitionJacobianF_( this->X_est_, dt );
#ifdef DEBUG_KALMAN
        ( this->F_ ).print( this->filterName_ + " Prediction, F:" );
        ( this->P_ ).print( this->filterName_ + " Prediction, P before:" );
#endif
        // 2. Вычисление ковариационной матрицы P
        arma::mat Qdt = this->Q_ * std::sqrt( std::abs( dt ) );
        this->P_ = this->F_ * this->P_;
        arma::mat QR_input_P_pred = arma::trans( arma::join_horiz( this->P_, Qdt ) ); // QR_input = [ P_, Qdt ], Транспонировать, т.к. в QR разложение так надо
        arma::mat Q_qr_P_pred, R_qr_P_pred;
        int res_JQR_P_pred = SPML::QR::ModifiedGramSchmidt( Q_qr_P_pred, R_qr_P_pred, QR_input_P_pred ); // Вызов QR разложения
        assert( res_JQR_P_pred == 0 );
        this->P_ = arma::trans( R_qr_P_pred ); // Матрица Р считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Prediction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );
        // 3. Вычисление матрицы перехода H
        this->H_ = this->observationJacobianH_( this->X_pred_ );

        // 4. Вычисление X_est, Y_est
        this->Y_pred_ = this->observationModel_( this->X_pred_ );
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_pred_ = this->checkBordersMeasurement_( this->Y_pred_ );
        }
        this->X_est_ = this->X_pred_;
        this->Y_est_ = this->Y_pred_;
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( this->filterName_ + " Prediction, X_est:" );
        ( this->Y_est_ ).print( this->filterName_ + " Prediction, Y_est:" );
#endif
        this->prediction_isDone = true; // Выставить признак состоявшегося прогноза
    }
};

}

#endif // KALMAN_FILTER_EXTENDED_SQUARE_ROOT_H
/// \}
