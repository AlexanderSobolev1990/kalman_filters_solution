//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_extended_cubature_square_root.h
/// \brief      Шаблонный класс квадратно-корневого расширенного кубатурного фильтра Калмана, КК-РКФК
/// (Square Root Extended Cubature Kalman Filter, SR-ECKF)
/// \date       30.03.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_CUBATURE_SQUARE_ROOT_H
#define KALMAN_FILTER_EXTENDED_CUBATURE_SQUARE_ROOT_H

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_extended.h>
#include <kalman_filter_extended_square_root.h>
#include <kalman_filter_cubature_square_root.h>
#include <qr_decomposition.h> // JQR разложение матрицы

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс квадратно-корневого расширенного кубатурного фильтра Калмана, КК-РКФК
/// (Square Root Extended Cubature Kalman Filter, SR-ECKF)
/// \details Источники:
/// \n[1] Cubature Kalman Filters, Ienkaran Arasaratnam and Simon Haykin, Life Fellow, IEEE
/// \attention Фильтр построен по классическому его варианту НИЖНЕтреугольного разложения Холецкого!
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRECKF : public CKalmanSREKF<SizeX, SizeY>, public CKalmanSRCKF<SizeX, SizeY>
{
public:
    using CKalmanSREKF<SizeX, SizeY>::CKalmanSREKF;
    using CKalmanSRCKF<SizeX, SizeY>::CKalmanSRCKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSRECKF() : CKalmanSREKF<SizeX, SizeY>(), CKalmanSRCKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SRECKF" );
#endif
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз квадратно-корневого расширенного кубатурного фильтра Калмана (КК-РКФК, SR-ECKF)
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
        this->PredictionSREKF( dt );
    }

    ///
    /// \brief Коррекция квадратно-корневого расширенного кубатурного фильтра Калмана (КК-РКФК, SR-ECKF)
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
        // 1a. Создание сигма-точек
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->P_.col(i) );
            this->x_pred_sigma_points_.col(i) = this->X_pred_ + add;
            this->x_pred_sigma_points_.col(i + SizeX) = this->X_pred_ - add;
            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
                this->x_pred_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i) );
                this->x_pred_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i + SizeX) );
            }
        }
        // 1b. Вычисление сигма-точек пространства Y
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }

        this->CorrectionSRCKF( Y_msd );
    }
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс квадратно-корневого расширенного кубатурного фильтра Калмана (блочная фильтрация), КК-РКФКБ
/// (Square Root Extended Cubature Kalman Filter Block, SR-ECKFB)
/// \details См. CKalmanSRECKF и
/// The J-Orthogonal Square-Root Euler-Maruyama-Based Unscented Kalman Filter for Nonliear Stohastic Systems,
/// Gennady Yu. Kulikov, Maria V. Kulikova, CEMAT, Instituto Superior Técnico, Universidade de Lisboa, Av.
/// Rovisco Pais 1, 1049-001 LISBOA, Portugal (emails: gennady.kulikov[at]tecnico.ulisboa.pt, maria.kulikova[at]ist.utl.pt)
/// \attention Фильтр построен по классическому его варианту НИЖНЕтреугольного разложения Холецкого!
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRECKFB : public CKalmanSRECKF<SizeX, SizeY>
{
public:
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSRECKFB() : CKalmanSRECKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SRECKFB" );
#endif
    }

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

        // 1a. Создание сигма-точек
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->P_.col(i) );
            this->x_pred_sigma_points_.col(i) = this->X_pred_ + add;
            this->x_pred_sigma_points_.col(i + SizeX) = this->X_pred_ - add;
            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
                this->x_pred_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i) );
                this->x_pred_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i + SizeX) );
            }
        }
        // 1b. Вычисление сигма-точек пространства Y
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }

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
        // 2. Выполнить JQR разложение блочной матрицы B и считать результат:
        // [ S     0 ]
        // [ Pxy*  P ]
        size_t QRsizeY = ( SizeY + ( 2 * SizeX ) );
        arma::mat Q_qr_corr, R_qr_corr;
        arma::mat JcorrectBlock( QRsizeY, QRsizeY, arma::fill::eye );
        int res_JQR_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_corr, R_qr_corr, B );
        assert( res_JQR_corr == 0 );
        R_qr_corr = arma::trans( R_qr_corr ); // Транспонировать, т.к. в JQR разложение так надо
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

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс квадратно-корневого гибридного расширенно-кубатурного фильтра Калмана с блочной коррекцией, КК-РКФКБ
/// (Square Root Extended Cubature Kalman Filter Block, SR-ECKFB)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRECKFBpure : public CKalmanEKF<SizeX, SizeY>
{
public:
    using CKalmanEKF<SizeX, SizeY>::CKalmanEKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSRECKFBpure() : CKalmanEKF<SizeX, SizeY>()
        , k_sigma_points_( 2 * SizeX )
        , weights_mean_( arma::zeros( this->k_sigma_points_ ) )
        , weights_covariance_( arma::zeros( this->k_sigma_points_ ) )
        , x_est_sigma_points_( arma::zeros( SizeX, this->k_sigma_points_ ) )
        , x_pred_sigma_points_( arma::zeros( SizeX, this->k_sigma_points_ ) )
        , y_pred_sigma_points_( arma::zeros( SizeY, this->k_sigma_points_ ) )
        , dXcal_( arma::zeros( SizeX, this->k_sigma_points_ ) )
        , dYcal_( arma::zeros( SizeY, this->k_sigma_points_ ) )
        , P_xy_( arma::zeros( SizeX, SizeY ) )
        , sqrt_P_chol_( arma::zeros( SizeX, SizeX ) )
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SRECKFB" );
#endif
        this->SetupDesignParametersCubatureBaseSet();
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы-сеттеры:

    ///
    /// \brief Установка функции вычисления взвешенной суммы сигма-точек пространства Х
    /// \param weightedSumStateSigmas - Функция вычисления взвешенной суммы сигма-точек пространства Х
    ///
    void SetWeightedSumStateSigmas( std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumStateSigmas )
    {
        weightedSumStateSigmas_ = weightedSumStateSigmas;
    }

    ///
    /// \brief Установка функции вычисления взвешенной суммы сигма-точек пространства Y
    /// \param weightedSumMeasurementSigmas - Функция вычисления взвешенной суммы сигма-точек пространства Y
    ///
    void SetWeightedSumMeasurementSigmas( std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumMeasurementSigmas )
    {
        weightedSumMeasurementSigmas_ = weightedSumMeasurementSigmas;
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
        arma::mat Qdt = this->Q_ * std::sqrt( dt );
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

        // 1a. Создание сигма-точек
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->P_.col(i) );
            this->x_pred_sigma_points_.col(i) = this->X_pred_ + add;
            this->x_pred_sigma_points_.col(i + SizeX) = this->X_pred_ - add;
            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
                this->x_pred_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i) );
                this->x_pred_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i + SizeX) );
            }
        }
        // 1b. Вычисление сигма-точек пространства Y
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }

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
        // 2. Выполнить JQR разложение блочной матрицы B и считать результат:
        // [ S     0 ]
        // [ Pxy*  P ]
        size_t QRsizeY = ( SizeY + ( 2 * SizeX ) );
        arma::mat Q_qr_corr, R_qr_corr;
        arma::mat JcorrectBlock( QRsizeY, QRsizeY, arma::fill::eye );
        int res_JQR_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_corr, R_qr_corr, B );
        assert( res_JQR_corr == 0 );
        R_qr_corr = arma::trans( R_qr_corr ); // Транспонировать, т.к. в JQR разложение так надо
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

protected:
    //------------------------------------------------------------------------------------------------------------------
    // Параметры, зависящие от SizeX, SizeY:
    int k_sigma_points_;            ///< Число сигма-точек
    arma::vec weights_mean_;        ///< Веса среднего
    arma::vec weights_covariance_;  ///< Веса ковариации
    arma::mat x_est_sigma_points_;  ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве X на текущем такте, размерность [SizeX,k_sigma_points_]
    arma::mat x_pred_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве X, экстраполированный на текущий такт, размерность [SizeX,k_sigma_points_]
    arma::mat y_pred_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве Y, экстраполированный на текущий такт, размерность [SizeX,k_sigma_points_]

    arma::mat dXcal_;       ///< Матрица Х-каллиграфическое (матрица сигма-точек - столбцов)
    arma::mat dYcal_;       ///< Матрица Y-каллиграфическое (матрица сигма-точек - столбцов)
    arma::mat P_xy_;        ///< Матрица кросс-коварации векторов Х и Y, размерность [SizeX * SizeY]
    arma::mat sqrt_P_chol_; ///< Корень из матрицы P
    //------------------------------------------------------------------------------------------------------------------
    // Параметры выбора сигма-точек, определяемые в методах SetDesignParameters*
    double gamma_;  ///< Автоматически вычисляемый (в методах SetDesignParameters*) параметр (множитель при корне из P при создании сигма-точек)

    //------------------------------------------------------------------------------------------------------------------
    // Обертки функций вычисления взвешенной суммы:
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumStateSigmas_; ///< Вычисление взвешенной суммы сигма-точек пространства Х
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumMeasurementSigmas_; ///< Вычисление взвешенной суммы сигма-точек пространства Y

    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Установка кубатурных весов (базовый вариант ансцентного преобразования)
    /// \details Смотри [1] и [2]
    ///
    void SetupDesignParametersCubatureBaseSet()
    {
        double gammaSq = static_cast<double>( SizeX );
        this->gamma_ = std::sqrt( gammaSq );

        // Set the weights for sigma points
        double all_points = 1.0 / ( 2.0 * gammaSq );
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->weights_mean_( i ) = all_points;
            this->weights_covariance_( i ) = all_points;
        }

#ifdef DEBUG_KALMAN
        std::cout << this->filterName_ + " SetDesignParameters_Cubature:" << std::endl;
        std::cout << "weights_mean_:" << std::endl;
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            std::cout << this->weights_mean_( i ) << std::endl;
        }
        std::cout << "weights_covariance_:" << std::endl;
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            std::cout << this->weights_covariance_( i ) << std::endl;
        }
#endif
    }

};

}

#endif // KALMAN_FILTER_EXTENDED_CUBATURE_SQUARE_ROOT_H
/// \}
