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
    // default copy/move/assignment semantic:
//    CKalmanSRECKF( const CKalmanSRECKF& ) = default;
//    CKalmanSRECKF& operator=( const CKalmanSRECKF& ) = default;
//    CKalmanSRECKF( CKalmanSRECKF&& ) = default;
//    CKalmanSRECKF& operator=( CKalmanSRECKF&& ) = default;
//    virtual ~CKalmanSRECKF() = default;

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
    // default copy/move/assignment semantic:
//    CKalmanSRECKFB( const CKalmanSRECKFB& ) = default;
//    CKalmanSRECKFB& operator=( const CKalmanSRECKFB& ) = default;
//    CKalmanSRECKFB( CKalmanSRECKFB&& ) = default;
//    CKalmanSRECKFB& operator=( CKalmanSRECKFB&& ) = default;
//    virtual ~CKalmanSRECKFB() = default;

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

}

#endif // KALMAN_FILTER_EXTENDED_CUBATURE_SQUARE_ROOT_H
/// \}
