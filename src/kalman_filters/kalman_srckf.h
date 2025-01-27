//------------------------------------------------------------------------------
///
/// \file       kalman_srckf.h
/// \brief      Класс квадратно-корневого кубатурного фильтра Калмана,
///             КК-КФК (square root cubature Kalman filter, SR-CKF)
/// \date       24.03.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_CUBATURE_SQUARE_ROOT_H
#define KALMAN_FILTER_CUBATURE_SQUARE_ROOT_H

// Project includes:
#include <kalman_debug.h>
#include <kalman_ckf.h>
#include <qr_decomposition.h>

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief Класс квадратно-корневого кубатурного фильтра Калмана, КК-КФК
///        (square root cubature Kalman filter, SR-CKF)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRCKF : public CKalmanCKF<SizeX, SizeY>
{
public:
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanSRCKF() = default;
    CKalmanSRCKF( const CKalmanSRCKF& ) = default;
    CKalmanSRCKF& operator=( const CKalmanSRCKF& ) = default;
    CKalmanSRCKF( CKalmanSRCKF&& ) = default;
    CKalmanSRCKF& operator=( CKalmanSRCKF&& ) = default;
    virtual ~CKalmanSRCKF() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "SRCKF"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Cubature; }
    virtual TCovType GetCovType() const override { return TCovType::CT_SquareRootMatrices; }    

    //--------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt ) override
    {
#ifdef DEBUG_KALMAN
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << GetFilterName() + " Prediction started, dt = " << dt << std::endl;
#endif
        // Создание сигма-точек пространства X
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( GetFilterName() + " Prediction, P before:" );
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
        ( this->x_est_sigma_points_ ).print( GetFilterName() + " Prediction, x_est_sigma_points_:" );
#endif
        // Прогноз сигма-точек пространства X на текущий такт
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->x_pred_sigma_points_.col(i) = this->stateTransitionModel_( this->x_est_sigma_points_.col(i), dt );
        }
#ifdef DEBUG_KALMAN
        ( this->x_pred_sigma_points_ ).print( GetFilterName() + " Prediction, x_pred_sigma_points_:" );
#endif
        // Вычисление X_pred по сигма-точкам пространства Х
        if( this->weightedSumStateSigmas_ == nullptr ) {
            this->X_pred_ = this->x_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->X_pred_ = this->weightedSumStateSigmas_( this->weights_mean_, this->x_pred_sigma_points_ );
        }
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }
#ifdef DEBUG_KALMAN
        ( this->X_pred_ ).print( GetFilterName() + " Prediction, X_pred_:" );
#endif
        // Вычисление матрицы Р        
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = this->x_pred_sigma_points_.col(i) - this->X_pred_;
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->dXcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }
        arma::mat Qdt = this->Q_ * std::sqrt( std::abs( dt ) );
        arma::mat B = arma::join_horiz( Qdt, this->dXcal_ ); // Блочная матрица
#ifdef DEBUG_KALMAN        
        B.print( GetFilterName() + " Prediction, B:" );
#endif  
        SPML::LQ::MGS_1( this->P_, B );
#ifdef DEBUG_KALMAN        
        ( this->P_ ).print( GetFilterName() + " Prediction, P after:" );
        ( this->P_ * arma::trans( this->P_ ) ).print( GetFilterName() + " Prediction, Pfull after:" );
#endif
        assert( this->indexOfNegativeDiagonalElement( this->P_ ) < 0 );

        // Пересоздание сигма-точек пространства Х после прогноза по новой матрице P
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = this->gamma_ * this->P_.col(i);
            this->x_pred_sigma_points_.col(i) = this->X_pred_ + add;
            this->x_pred_sigma_points_.col(i + SizeX) = this->X_pred_ - add;
            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
                this->x_pred_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i) );
                this->x_pred_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i + SizeX) );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->x_pred_sigma_points_ ).print( GetFilterName() + " Prediction, x_pred_sigma_points_:" );
#endif
        // Вычисление X_pred заново по пересозданным сигма-точкам пространства Х
        if( this->weightedSumStateSigmas_ == nullptr ) {
            this->X_pred_ = this->x_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->X_pred_ = this->weightedSumStateSigmas_( this->weights_mean_, this->x_pred_sigma_points_ );
        }
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }

        // Вычисление сигма-точек пространства Y по сигма-точкам пространства Х
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->y_pred_sigma_points_ ).print( GetFilterName() + " Prediction, y_pred_sigma_points_:" );
#endif
        // Вычисление Y_pred
        if( this->weightedSumMeasurementSigmas_ == nullptr ) {
            this->Y_pred_ = this->y_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->Y_pred_ = this->weightedSumMeasurementSigmas_( this->weights_mean_, this->y_pred_sigma_points_ );
        }
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_pred_ = this->checkBordersMeasurement_( this->Y_pred_ );
        }

        // Обновление X_est, Y_est
        this->X_est_ = this->X_pred_;
        this->Y_est_ = this->Y_pred_;
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( GetFilterName() + " Prediction, X_est:" );
        ( this->Y_est_ ).print( GetFilterName() + " Prediction, Y_est:" );
#endif
        this->prediction_isDone = true; // Выставить признак состоявшегося прогноза
    }

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd ) override
    {
        this->correctionSRCKF( Y_msd );
    }

protected:
    //--------------------------------------------------------------------------
    ///
    /// \brief Коррекция SRCKF
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    void correctionSRCKF( const arma::vec &Y_msd )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << GetFilterName() + " Correction started" << std::endl;
        Y_msd.print( GetFilterName() + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone );
        this->prediction_isDone = false;

        // Вычисление матрицы dXcal
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = this->x_pred_sigma_points_.col(i) - this->X_pred_;
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->dXcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }

        // Вычисление матрицы dYcal
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dYcal_.col(i) = ( this->y_pred_sigma_points_.col(i) - this->Y_pred_ );
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->dYcal_.col(i) = this->checkDeltaMeasurement_( this->dYcal_.col(i) );
            }
            this->dYcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }

        // Вычисление матрицы S
        arma::mat B1 = arma::join_horiz( this->R_, this->dYcal_ ); // Блочная матрица
#ifdef DEBUG_KALMAN
        B1.print( GetFilterName() + " Correction, B1:" );        
#endif
        SPML::LQ::MGS_1( this->S_, B1 );                        
#ifdef DEBUG_KALMAN        
        ( this->S_ ).print( GetFilterName() + " Correction, S after:" );
        ( this->S_ * arma::trans( this->S_ ) ).print( GetFilterName() + " Correction, Sfull after:" );
#endif
        assert( this->indexOfNegativeDiagonalElement( this->S_ ) < 0 );

        // Вычисление матрицы P_xy        
        this->P_xy_ = this->dXcal_ * arma::trans( this->dYcal_ );
#ifdef DEBUG_KALMAN
        ( this->dXcal_ ).print( GetFilterName() + " Prediction, dXcal_:" );
        ( this->dYcal_ ).print( GetFilterName() + " Prediction, dYcal_:" );
        ( this->P_xy_ ).print( GetFilterName() + " Correction, P_xy:" );
#endif
        // Вычисление коэффициента усиления фильтра K
        // Внимание! Библиотека Armadillo сама оптимизирует inv через solve!
        this->K_ = ( this->P_xy_ * arma::inv( arma::trans( this->S_ ) ) ) * arma::inv( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->K_ ).print( GetFilterName() + " Correction, K:" );
#endif
        // Вычисление невязки Delta
        if( !this->deltaY_isSet ) {
            assert( this->Y_msd_isSet ); // Если не установлен deltaY, то Y_msd_ обязан быть установлен
            this->Y_msd_isSet = false; // Сразу же снять признак
            this->DeltaY_ = this->Y_msd_ - this->Y_pred_;
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->DeltaY_ = this->checkDeltaMeasurement_( this->DeltaY_ );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->DeltaY_ ).print( GetFilterName() + " Correction, DeltaY:" );
#endif
        // Вычисление матрицы P
        arma::mat P1 = this->K_ * this->R_;
        arma::mat P2 = this->dXcal_ - ( this->K_ * this->dYcal_ );
        if( this->checkDeltaState_ != nullptr ) {
            for( int i = 0; i < this->k_sigma_points_; i++ ) {
                P2.col(i) = this->checkDeltaState_( P2.col(i) );
            }
        }
        arma::mat B2 = arma::join_horiz( P1, P2 ); // Блочная матрица
#ifdef DEBUG_KALMAN
        B2.print( GetFilterName() + " Correction, B2:" );
#endif
        SPML::LQ::MGS_1( this->P_, B2 );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( GetFilterName() + " Correction, P after:" );
        ( this->P_ * arma::trans( this->P_ ) ).print( GetFilterName() + " Correction, Pfull after:" );
#endif
        assert( this->indexOfNegativeDiagonalElement( this->P_ ) < 0 );

        // Вычисление X_est, Y_est
        this->X_est_ = this->X_pred_ + ( this->K_ * this->DeltaY_ );
        if( this->checkBordersStateAfterCorrection_ != nullptr ) {
            this->X_est_ = this->checkBordersStateAfterCorrection_( this->X_est_ );
        }
        this->Y_est_ = this->observationModel_( this->X_est_ );
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_est_ = this->checkBordersMeasurement_( this->Y_est_ );
        }
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( GetFilterName() + " Correction, X_est:" );
        ( this->Y_est_ ).print( GetFilterName() + " Correction, Y_est:" );
#endif
        this->deltaY_isSet = false; // Снять признак (выставляется в true в сеттере deltaY)
    }

    ///
    /// \brief Коррекция SRCKFB
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    void correctionSRCKFB( const arma::vec &Y_msd )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << GetFilterName() + " Correction started" << std::endl;
        Y_msd.print( GetFilterName() + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла
        this->prediction_isDone = false; // Сразу же снять признак

        // Вычисление матрицы dXcal
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = ( this->x_pred_sigma_points_.col(i) - this->X_pred_ );
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->dXcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }

        // Вычисление матрицы dYcal
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dYcal_.col(i) = ( this->y_pred_sigma_points_.col(i) - this->Y_pred_ );
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->dYcal_.col(i) = this->checkDeltaMeasurement_( this->dYcal_.col(i) );
            }
            this->dYcal_.col(i) *= std::sqrt( this->weights_covariance_(i) );
        }        

        // Составление блочной матрицы B:
        // B = [ R  dYcal ]
        //     [ 0  dXcal ]
        arma::mat B = arma::join_vert(
            arma::join_horiz( this->R_, this->dYcal_ ),
            arma::join_horiz( arma::mat( SizeX, SizeY, arma::fill::zeros ), this->dXcal_ )
        );
#ifdef DEBUG_KALMAN
        B.print( GetFilterName() + " Correction, B:" );
#endif
        // QR-разложение составленной блочной матрицы B       
        arma::mat B_output;
        SPML::LQ::MGS_1( B_output, B );

        // Считывание результата из матрицы B_output: 
        //  X.submat( first_row, first_col, last_row, last_col )
        // [ S     0 ] = B_output
        // [ Pxy*  P ]
        int qrSize = SizeX + SizeY;
        this->S_ = B_output.submat( 0, 0, ( SizeY - 1 ), ( SizeY - 1 ) );
        this->P_xy_ = B_output.submat( SizeY, 0, ( qrSize - 1 ), ( SizeY - 1 ) );
        this->P_ = B_output.submat( SizeY, SizeY, ( qrSize - 1 ), ( qrSize - 1 ) );
#ifdef DEBUG_KALMAN        
        B_output.print( GetFilterName() + " Correction, B_output:" );
        ( this->S_ ).print( GetFilterName() + " Correction, S_:" );
        ( this->P_xy_ ).print( GetFilterName() + " Correction, P_xy:" );
        ( this->P_ ).print( GetFilterName() + " Correction, P_:" );
#endif
        assert( this->indexOfNegativeDiagonalElement( this->P_ ) < 0 );
        assert( this->indexOfNegativeDiagonalElement( this->S_ ) < 0 );

        // Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_xy_ * arma::inv( this->S_ ); // В блочном случае - это ВЕРНО!
#ifdef DEBUG_KALMAN
        ( this->K_ ).print( GetFilterName() + " Correction, K:" );
#endif
        // Вычисление невязки Delta
        if( !this->deltaY_isSet ) {
            assert( this->Y_msd_isSet ); // Если не установлен deltaY, то Y_msd_ обязан быть установлен
            this->Y_msd_isSet = false; // Сразу же снять признак
            this->DeltaY_ = this->Y_msd_ - this->Y_pred_;
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->DeltaY_ = this->checkDeltaMeasurement_( this->DeltaY_ );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->DeltaY_ ).print( GetFilterName() + " Correction, DeltaY:" );
#endif
        // Вычисление X_est, Y_est
        this->X_est_ = this->X_pred_ + ( this->K_ * this->DeltaY_ );
        if( this->checkBordersStateAfterCorrection_ != nullptr ) {
            this->X_est_ = this->checkBordersStateAfterCorrection_( this->X_est_ );
        }
        this->Y_est_ = this->observationModel_( this->X_est_ );
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_est_ = this->checkBordersMeasurement_( this->Y_est_ );
        }
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( GetFilterName() + " Correction, X_est:" );
        ( this->Y_est_ ).print( GetFilterName() + " Correction, Y_est:" );
#endif
        this->deltaY_isSet = false; // Снять признак (выставляется в true в сеттере deltaY)
    }
};

//------------------------------------------------------------------------------
///
/// \brief Класс квадратно-корневого кубатурного фильтра Калмана
///        с блочной коррекцией, КК-КФКБ (square root cubature Kalman filter
///        with block correction, SR-CKFB)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRCKFB : public CKalmanSRCKF<SizeX, SizeY>
{
public:
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanSRCKFB() = default;
    CKalmanSRCKFB( const CKalmanSRCKFB& ) = default;
    CKalmanSRCKFB& operator=( const CKalmanSRCKFB& ) = default;
    CKalmanSRCKFB( CKalmanSRCKFB&& ) = default;
    CKalmanSRCKFB& operator=( CKalmanSRCKFB&& ) = default;
    virtual ~CKalmanSRCKFB() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "SRCKFB"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Cubature; }
    virtual TCovType GetCovType() const override { return TCovType::CT_SquareRootMatrices; }    

    //--------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    // Prediction как в НЕблочном фильтре

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd ) override
    {
        this->correctionSRCKFB( Y_msd );
    }
};

}
#endif // KALMAN_FILTER_CUBATURE_H
/// \}
