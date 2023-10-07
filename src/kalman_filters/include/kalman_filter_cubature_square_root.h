//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_cubature_square_root.h
/// \brief      Шаблонный класс квадратно-корневого кубатурного фильтра Калмана, КК-КФК
///             (Square Root Cubature Kalman Filter, SR-CKF)
/// \date       24.03.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_CUBATURE_SQUARE_ROOT_H
#define KALMAN_FILTER_CUBATURE_SQUARE_ROOT_H

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_cubature.h>
#include <qr_decomposition.h> // JQR разложение матрицы

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс квадратно-корневого кубатурного фильтра Калмана, КК-КФК
/// (Square Root Cubature Kalman Filter, SR-CKF)
/// \details Источники:
/// \n[1] Cubature Kalman Filters, Ienkaran Arasaratnam and Simon Haykin, Life Fellow, IEEE
/// \attention Фильтр построен по классическому его варианту НИЖНЕтреугольного разложения Холецкого!
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRCKF : public CKalmanCKF<SizeX, SizeY>
{
public:
    using CKalmanCKF<SizeX, SizeY>::CKalmanCKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSRCKF() : CKalmanCKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SRCKF" );
#endif
    }
    // default copy/move/assignment semantic:
//    CKalmanSRCKF( const CKalmanSRCKF& ) = default;
//    CKalmanSRCKF& operator=( const CKalmanSRCKF& ) = default;
//    CKalmanSRCKF( CKalmanSRCKF&& ) = default;
//    CKalmanSRCKF& operator=( CKalmanSRCKF&& ) = default;
//    virtual ~CKalmanSRCKF() = default;

    //------------------------------------------------------------------------------------------------------------------
    // Методы-сеттеры:

    ///
    /// \brief Установка функции вычисления матрицы перехода состояния F (makeMatrixF)
    /// \sa stateTransitionJacobianF_
    ///
    void SetStateTransitionJacobianF( std::function<arma::mat( const arma::vec &X, double dt )> stateTransitionJacobianF )
    {
        this->stateTransitionJacobianF_ = stateTransitionJacobianF;
    }

    ///
    /// \brief Установка функции вычисления матрицы перехода измерений H (makeMatrixH)
    /// \sa observationJacobianH_
    ///
    void SetObservationJacobianH( std::function<arma::mat( const arma::vec &X )> observationJacobianH )
    {
        this->observationJacobianH_ = observationJacobianH;
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Prediction started, dt = " << dt << std::endl;
#endif
        assert( dt > 0.0 ); // Прогноз имеет смысл только при положительном времени

        // 1. Создание сигма-точек пространства X
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Prediction, P before:" );
#endif        
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->P_.col(i) );
            this->x_est_sigma_points_.col(i) = this->X_est_ + add;
            this->x_est_sigma_points_.col(i + SizeX) = this->X_est_ - add;
            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
                this->x_est_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_est_sigma_points_.col(i) );
                this->x_est_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_est_sigma_points_.col(i + SizeX) );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->x_est_sigma_points_ ).print( this->filterName_ + " Prediction, x_est_sigma_points_:" );
#endif
        // 2. Прогноз сигма-точек пространства X на текущий такт
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->x_pred_sigma_points_.col(i) = this->stateTransitionModel_( this->x_est_sigma_points_.col(i), dt );
        }
#ifdef DEBUG_KALMAN
        ( this->x_pred_sigma_points_ ).print( this->filterName_ + " Prediction, x_pred_sigma_points_:" );
#endif
        // 3. Вычисление X_pred по сигма-точкам пространства Х
        if( this->weightedSumStateSigmas_ == nullptr ) {
            this->X_pred_ = this->x_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->X_pred_ = this->weightedSumStateSigmas_( this->weights_mean_, this->x_pred_sigma_points_ );
        }
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }
#ifdef DEBUG_KALMAN
        ( this->X_pred_ ).print( this->filterName_ + " Prediction, X_pred_:" );
#endif
        // 4. Вычисление ковариационной матрицы Р
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = ( this->x_pred_sigma_points_.col(i) - this->X_pred_ );
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->dXcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }
        arma::mat Qdt = this->Q_ * std::sqrt( std::abs( dt ) );
        arma::mat QR_input_P_pred = arma::trans( arma::join_horiz( this->dXcal_, Qdt ) ); // [ dXcal, Qdt ]', Транспонировать, т.к. в QR разложение так надо
#ifdef DEBUG_KALMAN
        QR_input_P_pred.print( this->filterName_ + " Prediction, QR_input_P_pred:" );
#endif
        arma::mat Q_qr_P_pred, R_qr_P_pred;
        int res_QR_P_pred = SPML::QR::ModifiedGramSchmidt( Q_qr_P_pred, R_qr_P_pred, QR_input_P_pred ); // Вызов QR разложения
        assert( res_QR_P_pred == 0 );
        this->P_ = arma::trans( R_qr_P_pred ); // Матрица Р считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        Q_qr_P_pred.print( this->filterName_ + " Prediction, Q_qr_P_pred after:" );
        ( this->P_ ).print( this->filterName_ + " Prediction, P after:" );
        ( this->P_ * arma::trans( this->P_ ) ).print( this->filterName_ + " Prediction, Pfull after:" );
#endif
        // 5. Пересоздание сигма-точек после прогноза по новой P
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->P_.col(i) );
            this->x_pred_sigma_points_.col(i) = this->X_pred_ + add;
            this->x_pred_sigma_points_.col(i + SizeX) = this->X_pred_ - add;
            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
                this->x_pred_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i) );
                this->x_pred_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i + SizeX) );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->x_pred_sigma_points_ ).print( this->filterName_ + " Prediction, x_pred_sigma_points_:" );
#endif
        // 6. Вычисление X_pred заново по пересозданным сигма-точкам пространства Х (в [1] не указано, но это подразумевается)
        if( this->weightedSumStateSigmas_ == nullptr ) {
            this->X_pred_ = this->x_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->X_pred_ = this->weightedSumStateSigmas_( this->weights_mean_, this->x_pred_sigma_points_ );
        }
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }

        // 7. Вычисление сигма-точек пространства Y
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->y_pred_sigma_points_ ).print( this->filterName_ + " Prediction, y_pred_sigma_points_:" );
#endif
        // 8. Вычисление предсказанного вектора Y_pred на текущий такт по сигма-точкам пространства Х
        if( this->weightedSumMeasurementSigmas_ == nullptr ) {
            this->Y_pred_ = this->y_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->Y_pred_ = this->weightedSumMeasurementSigmas_( this->weights_mean_, this->y_pred_sigma_points_ );
        }
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_pred_ = this->checkBordersMeasurement_( this->Y_pred_ );
        }

        // 9. Обновление X_est, Y_est
        this->X_est_ = this->X_pred_;
        this->Y_est_ = this->Y_pred_;
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( this->filterName_ + " Prediction, X_est:" );
        ( this->Y_est_ ).print( this->filterName_ + " Prediction, Y_est:" );
#endif
        this->prediction_isDone = true; // Выставить признак состоявшегося прогноза
    }

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
        this->CorrectionSRCKF( Y_msd );
    }

protected:
    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Коррекция SRCKF
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    void CorrectionSRCKF( const arma::vec &Y_msd )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Correction started" << std::endl;
        Y_msd.print( this->filterName_ + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла
        this->prediction_isDone = false; // Сразу же снять признак

        // 1. Вычисление ковариационной матрицы S (P_yy)
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dYcal_.col(i) = ( this->y_pred_sigma_points_.col(i) - this->Y_pred_ );
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->dYcal_.col(i) = this->checkDeltaMeasurement_( this->dYcal_.col(i) );
            }
            this->dYcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }
        arma::mat QR_input_S_corr = arma::trans( arma::join_horiz( this->dYcal_, this->R_ ) ); //[ R, YpredW ]', Транспонировать, т.к. в QR разложение так надо
#ifdef DEBUG_KALMAN
        QR_input_S_corr.print( this->filterName_ + " Correction, QR_input_S_corr:" );
#endif
        arma::mat Q_qr_S_corr, R_qr_S_corr;
        int res_QR_S_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_S_corr, R_qr_S_corr, QR_input_S_corr ); // Вызов QR разложения
        assert( res_QR_S_corr == 0 );
        this->S_ = arma::trans( R_qr_S_corr ); // Матрица S считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        Q_qr_S_corr.print( this->filterName_ + " Correction, Q_qr_S_corr after:" );
        ( this->S_ ).print( this->filterName_ + " Correction, S after:" );
        ( this->S_ * arma::trans( this->S_ ) ).print( this->filterName_ + " Correction, Sfull after:" );
#endif
        // 2. Вычисление кросс-ковариационной матрицы P_xy
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = ( this->x_pred_sigma_points_.col(i) - this->X_pred_ );
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->dXcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }
        this->P_xy_ = this->dXcal_ * arma::trans( this->dYcal_ );
#ifdef DEBUG_KALMAN
        ( this->P_xy_ ).print( this->filterName_ + " Correction, P_xy:" );
#endif
        // 3. Вычисление коэффициента усиления фильтра K

        // 1 способ
        this->K_ = ( this->P_xy_ * arma::inv( arma::trans( this->S_ ) ) ) * arma::inv( this->S_ );

        // 2 способ
       // this->K_ = arma::trans( arma::solve( arma::trans( this->S_ ), arma::solve( this->S_, arma::trans( P_xy ) ) ) );

#ifdef DEBUG_KALMAN
        ( this->K_ ).print( this->filterName_ + " Correction, K:" );
#endif
        // 4. Вычисление невязки Delta
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
        // 5. Вычисление ковариационной матрицы P
        arma::mat P1 = this->K_ * this->R_;
        arma::mat P2 = this->dXcal_ - ( this->K_ * this->dYcal_ );
        if( this->checkDeltaState_ != nullptr ) {
            for( int i = 0; i < this->k_sigma_points_; i++ ) {
                P2.col(i) = this->checkDeltaState_( P2.col(i) );
            }
        }
        arma::mat QR_input_P_corr = arma::trans( arma::join_horiz( P1, P2 ) ); // Транспонировать, т.к. в QR разложение так надо
#ifdef DEBUG_KALMAN
        QR_input_P_corr.print( this->filterName_ + " Correction, QR_input_P_corr:" );
#endif
        arma::mat Q_qr_P_corr, R_qr_P_corr;
        int res_JQR_P_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_P_corr, R_qr_P_corr, QR_input_P_corr ); // Вызов QR разложения
        assert( res_JQR_P_corr == 0 );
        this->P_ = arma::trans( R_qr_P_corr ); // Матрица P считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        Q_qr_P_corr.print( this->filterName_ + " Correction, Q_qr_P_corr after:" );
        ( this->P_ ).print( this->filterName_ + " Correction, P after:" );
        ( this->P_ * arma::trans( this->P_ ) ).print( this->filterName_ + " Correction, Pfull after:" );
#endif
        // 6. Вычисление X_est, Y_est
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
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс квадратно-корневого кубатурного фильтра Калмана (блочная фильтрация),
/// КК-КФКБ (Square Root Cubature Kalman Filter Block, SR-CKFB)
/// \details См. CKalmanSRCKF и
/// The J-Orthogonal Square-Root Euler-Maruyama-Based Unscented Kalman Filter for Nonliear Stohastic Systems,
/// Gennady Yu. Kulikov, Maria V. Kulikova, CEMAT, Instituto Superior Técnico, Universidade de Lisboa, Av.
/// Rovisco Pais 1, 1049-001 LISBOA, Portugal (emails: gennady.kulikov[at]tecnico.ulisboa.pt, maria.kulikova[at]ist.utl.pt)
/// \attention Фильтр построен по классическому его варианту НИЖНЕтреугольного разложения Холецкого!
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRCKFB : public CKalmanSRCKF<SizeX, SizeY>
{
public:
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSRCKFB() : CKalmanSRCKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SRCKFB" );
#endif
    }
    // default copy/move/assignment semantic:
//    CKalmanSRCKFB( const CKalmanSRCKFB& ) = default;
//    CKalmanSRCKFB& operator=( const CKalmanSRCKFB& ) = default;
//    CKalmanSRCKFB( CKalmanSRCKFB&& ) = default;
//    CKalmanSRCKFB& operator=( CKalmanSRCKFB&& ) = default;
//    virtual ~CKalmanSRCKFB() = default;

    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    // Prediction как в НЕблочном фильтре

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

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла
        this->prediction_isDone = false; // Сразу же снять признак

        // 1. Составление блочной матрицы:
        // B = [ R    dYcal * WcovSqrtDiag ]
        //     [ 0    dXcal * WcovSqrtDiag ]
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dYcal_.col(i) = ( this->y_pred_sigma_points_.col(i) - this->Y_pred_ );
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->dYcal_.col(i) = this->checkDeltaMeasurement_( this->dYcal_.col(i) );
            }
            this->dYcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = ( this->x_pred_sigma_points_.col(i) - this->X_pred_ );
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->dXcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }
#ifdef DEBUG_KALMAN
        ( this->dXcal_ ).print( this->filterName_ + " Correction, dXcal_" );
        ( this->dYcal_ ).print( this->filterName_ + " Correction, dYcal_" );
#endif
        size_t rowsB = ( SizeY + SizeX );
        size_t colsB = ( SizeY + ( 2 * SizeX ) );
        arma::mat B = arma::mat( rowsB, colsB, arma::fill::zeros );
        for( size_t i = 0; i < rowsB; i++ ) {
            for( size_t j = 0; j < colsB; j++ ) {
                if( ( i < SizeY ) && ( j < SizeY ) ) { // Блок R
                    B( i, j ) = this->R_( i, j );
                } else if( ( i < SizeY ) && ( j >= SizeY ) ) { // Блок dYcal
                    B( i, j ) = this->dYcal_( i, j - SizeY );
                } else if( ( i >= SizeY ) && ( j < SizeY ) ) { // Нулевой блок
                    B( i, j ) = 0.0;
                } else { // Блок dXcal if( ( i >= SizeY ) && ( j >= SizeY ) )
                    B( i, j ) = this->dXcal_( i - SizeY, j - SizeY );
                }
            }
        }
        B = arma::trans( B ); // Транспонировать, т.к. в JQR разложение так надо
#ifdef DEBUG_KALMAN
        B.print( this->filterName_ + " Correction, B' ( block of R_, dYcal_, dXcal_ ):" );
#endif
        // 2. Выполнить QR разложение блочной матрицы B и считать результат:
        // [ S     0 ]
        // [ Pxy*  P ]        
        arma::mat Q_qr_corr, R_qr_corr;
        int res_JQR_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_corr, R_qr_corr, B );
        assert( res_JQR_corr == 0 );
        R_qr_corr = arma::trans( R_qr_corr ); // Транспонировать, т.к. в QR разложение так надо
        for( size_t i = 0; i < rowsB; i++ ) {
            for( size_t j = 0; j < rowsB; j++ ) { //for( int j = 0; j < colsB; j++ ) {
                if( ( i < SizeY ) && ( j < SizeY ) ) { // Матрица S
                   this->S_(i, j) = R_qr_corr(i, j);
                } else if( ( i >= SizeY ) && ( j < SizeY ) ) { // Матрица P_xy*
                    this->P_xy_(i - SizeY, j) = R_qr_corr(i, j);
                } else if( ( i >= SizeY ) && ( ( j >= SizeY ) && ( j < SizeY + SizeX ) ) ) { // Матрица P
                    this->P_(i - SizeY, j - SizeY) = R_qr_corr(i, j);
                }
            }
        }
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( this->filterName_ + " Correction, S_:" );
        ( this->P_xy_ ).print( this->filterName_ + " Correction, P_xy:" );
        ( this->P_ ).print( this->filterName_ + " Correction, P_:" );
#endif
        // 3. Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_xy_ * arma::inv( this->S_ ); // В блочном случае - это ВЕРНО!
#ifdef DEBUG_KALMAN
        ( this->K_ ).print( this->filterName_ + " Correction, K:" );
#endif
        // 4. Вычисление невязки Delta
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
        // 6. Вычисление X_est, Y_est
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
};

}

#endif // KALMAN_FILTER_CUBATURE_H
/// \}
