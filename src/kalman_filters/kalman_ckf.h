//------------------------------------------------------------------------------
///
/// \file       kalman_ckf.h
/// \brief      Класс кубатурного фильтра Калмана, КФК 
///             (cubature Kalman filter, CKF)
/// \date       24.03.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_CUBATURE_H
#define KALMAN_FILTER_CUBATURE_H

// Project includes:
#include <kalman_debug.h>
#include <kalman_ekf.h>
#include <qr_decomposition.h>

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief Класс кубатурного фильтра Калмана, КФК
///        (cubature Kalman filter, CKF)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanCKF : virtual public CKalmanEKF<SizeX, SizeY>
{
public:
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanCKF() = default;
    CKalmanCKF( const CKalmanCKF& ) = default;
    CKalmanCKF& operator=( const CKalmanCKF& ) = default;
    CKalmanCKF( CKalmanCKF&& ) = default;
    CKalmanCKF& operator=( CKalmanCKF&& ) = default;
    virtual ~CKalmanCKF() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "CKF"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Cubature; }
    virtual TCovType GetCovType() const override { return TCovType::CT_FullMatrices; }    

    //--------------------------------------------------------------------------
    // Методы-сеттеры:

    virtual void SetDesignParametersCubatureBaseSet() override
    {
        const double gammaSq = this->gamma_ * this->gamma_;

        // Set the weights for sigma points
        const double all_points = 1.0 / ( 2.0 * gammaSq );
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->weights_mean_( i ) = all_points;
            this->weights_covariance_( i ) = all_points;
        }

#ifdef DEBUG_KALMAN
        std::cout << GetFilterName() + " SetDesignParametersCubature:" << std::endl;
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
        this->sqrt_P_chol_ = arma::chol( this->P_, "lower" );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( GetFilterName() + " Prediction, P before:" );
        ( this->sqrt_P_chol_ ).print( GetFilterName() + " Prediction, sqrt_P_chol:" );
#endif
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->sqrt_P_chol_.col(i) );
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
        this->P_ = this->Q_ * std::abs( dt );
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = this->x_pred_sigma_points_.col(i) - this->X_pred_;
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->P_ += this->weights_covariance_(i) * this->dXcal_.col(i) * arma::trans( this->dXcal_.col(i) );
        }
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->dXcal_ ).print( GetFilterName() + " Prediction, dXcal_:" );
        ( this->P_ ).print( GetFilterName() + " Prediction, P after:" );
#endif        
        assert( this->indexOfNegativeDiagonalElement( this->P_ ) < 0 );

        // Пересоздание сигма-точек пространства Х после прогноза по новой матрице P
        this->sqrt_P_chol_ = arma::chol( this->P_, "lower" );
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->sqrt_P_chol_.col(i) );
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
        this->correctionCKF( Y_msd );
    }

protected:
    //--------------------------------------------------------------------------
    ///
    /// \brief Коррекция CKF
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    void correctionCKF( const arma::vec &Y_msd )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << GetFilterName() + " Correction started" << std::endl;
        Y_msd.print( GetFilterName() + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone );
        this->prediction_isDone = false;

        // Вычисление матрицы S
        this->S_ = this->R_;
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dYcal_.col(i) = this->y_pred_sigma_points_.col(i) - this->Y_pred_;
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->dYcal_.col(i) = this->checkDeltaMeasurement_( this->dYcal_.col(i) );
            }
            this->S_ += this->weights_covariance_(i) * this->dYcal_.col(i) * arma::trans( this->dYcal_.col(i) );
        }
        this->fixMatrixMainDiagonalSymmetry( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( GetFilterName() + " Correction, S:" );
#endif        
        assert( this->indexOfNegativeDiagonalElement( this->S_ ) < 0 );

        // Вычисление матрицы P_xy
        this->P_xy_.zeros();
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = this->x_pred_sigma_points_.col(i) - this->X_pred_;
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }            
            this->P_xy_ += this->weights_covariance_(i) * this->dXcal_.col(i) * arma::trans( this->dYcal_.col(i) );
        }
#ifdef DEBUG_KALMAN
        ( this->dXcal_ ).print( GetFilterName() + " Prediction, dXcal_:" );
        ( this->dYcal_ ).print( GetFilterName() + " Prediction, dYcal_:" );
        ( this->P_xy_ ).print( GetFilterName() + " Correction, P_xy:" );
#endif
        // Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_xy_ * arma::inv( this->S_ );
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
        this->P_ = this->P_ - ( this->K_ * this->S_ * arma::trans( this->K_ ) );
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( GetFilterName() + " Correction, P after:" );
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

    //--------------------------------------------------------------------------
    // Параметры, зависящие от SizeX, SizeY:
    static const int k_sigma_points_ = ( 2 * SizeX ); ///< Число сигма-точек
    static constexpr double gamma_ = std::sqrt( SizeX );  ///< Автоматически вычисляемый (в методах SetDesignParameters*) параметр (множитель при корне из P при создании сигма-точек)
    
    arma::vec::fixed<k_sigma_points_> weights_mean_; ///< Веса среднего
    arma::vec::fixed<k_sigma_points_> weights_covariance_; ///< Веса ковариации
    arma::mat::fixed<SizeX, k_sigma_points_> x_est_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве X на текущем такте, размерность [SizeX,k_sigma_points_]
    arma::mat::fixed<SizeX, k_sigma_points_> x_pred_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве X, экстраполированный на текущий такт, размерность [SizeX,k_sigma_points_]
    arma::mat::fixed<SizeY, k_sigma_points_> y_pred_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве Y, экстраполированный на текущий такт, размерность [SizeX,k_sigma_points_]

    arma::mat::fixed<SizeX, k_sigma_points_> dXcal_; ///< Матрица Х-каллиграфическое (матрица сигма-точек - столбцов)
    arma::mat::fixed<SizeY, k_sigma_points_> dYcal_; ///< Матрица Y-каллиграфическое (матрица сигма-точек - столбцов)
    arma::mat::fixed<SizeX, SizeY> P_xy_; ///< Матрица кросс-коварации векторов Х и Y, размерность [SizeX * SizeY]
    arma::mat::fixed<SizeX, SizeX> sqrt_P_chol_; ///< Корень из матрицы P
    
    //--------------------------------------------------------------------------
    // Обертки функций вычисления взвешенной суммы:
    // std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumStateSigmas_; ///< Вычисление взвешенной суммы сигма-точек пространства Х
    // std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumMeasurementSigmas_; ///< Вычисление взвешенной суммы сигма-точек пространства Y

    virtual void SetDesignParametersMeanSet( double w0 ){}
    virtual void SetDesignParametersScaledSet( double alpha, double beta, double kappa ){}
};

}
#endif // KALMAN_FILTER_CUBATURE_H
/// \}
