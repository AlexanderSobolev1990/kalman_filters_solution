//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_unscented.h
/// \brief      Шаблонный класс сигма-точечного (ансцентного) фильтра Калмана, СТФК (АФК) (Unscented Kalman Filter, UKF)
/// \date       17.01.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_UNSCENTED_H
#define KALMAN_FILTER_UNSCENTED_H

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_extended.h>

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс сигма-точечного (ансцентного) фильтра Калмана, СТФК (АФК) (Unscented Kalman Filter, UKF)
/// \details Источники:
/// \n[1] A New Extension of the Kalman Filter to Nonlinear Systems, Simon J. Julier, Jeffrey K. Uhlmann,
/// The Robotics Research Group, Department of Engineering Science, The University of Oxford, 1997
/// \n[2] Julier S.J., Uhlmann J.K. Unscented filtering and nonlinear estiation // Proc. of the IEEE, 2004, №3.P.401-422
/// \n[3] The Unscented Kalman Filter for Nonlinear Estimation, Eric A. Wan and Rudolph van der Merwe
/// Oregon Graduate Institute of Science & Technology 20000 NW Walker Rd, Beaverton, Oregon 97006, 2000
/// ericwan[at]ece.ogi.edu, rvdmerwe[at]ece.ogi.edu
/// \n[4] Sigma-Point Kalman Filters for Probabilistic Inference in Dynamic State-Space Models, Rudolph van der Merwe &
/// Eric Wan, OGI School of Science & Engineering Oregon Health & Science University Beaverton, Oregon, 97006, USA, 2003
/// {rvdmerwe,ericwan}[at]ece.ogi.edu
/// \n[5] THE SQUARE-ROOT UNSCENTED KALMAN FILTER FOR STATE AND PARAMETER-ESTIMATION, Rudolph van der Merwe and
/// Eric A. Wan, Oregon Graduate Institute of Science and Technology 20000 NW Walker Road, Beaverton, Oregon 97006, USA
/// rvdmerwe,ericwan[at]ece.ogi.edu
/// \n[6] Sebastian Bitzer, Technische Universität Dresden, https://github.com/sbitzer/UKF-exposed/blob/master/UKF.pdf
/// \attention Внимание! В реализации UKF вес нулевой сигма-точки Wcov не может быть отрицательным! (Wmean - может)
/// \attention Фильтр построен по классическому его варианту НИЖНЕтреугольного разложения Холецкого!
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanUKF : virtual public CKalmanEKF<SizeX, SizeY>
{
public:
    using CKalmanEKF<SizeX, SizeY>::CKalmanEKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanUKF() : CKalmanEKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "UKF" );
#endif
    }
    // default copy/move/assignment semantic:
    CKalmanUKF( const CKalmanUKF& ) = default;
    CKalmanUKF& operator=( const CKalmanUKF& ) = default;
    CKalmanUKF( CKalmanUKF&& ) = default;
    CKalmanUKF& operator=( CKalmanUKF&& ) = default;
    virtual ~CKalmanUKF() = default;

    //------------------------------------------------------------------------------------------------------------------
    // Методы-сеттеры:

    ///
    /// \brief Установка параметра w0 ансцентного фильтра (MeanSet)
    /// \details При выборе w0 = [0...1) обеспечиваются положительные веса. Смотри [1], [2]
    /// \attention Нельзя выбирать w0 так, чтобы нулевой вес был > 0, а остальные меньше нуля. Наоборот - МОЖНО, т.е.
    /// нулевой вес может быть отрицательным.
    /// \param w0 - параметр разброса сигма точек (w0 = [0...1) типичная рекомендация для положительных весов)
    ///
    virtual void SetupDesignParametersMeanSet( double w0 )
    {
        this->w0_ = w0;
        this->kappa_ = ( w0 * static_cast<double>( SizeX ) ) / ( 1.0 - w0 );
        double gammaSq = static_cast<double>( SizeX ) + this->kappa_;
        this->gamma_ = std::sqrt( gammaSq );

        // Set the weights for zero sigma point:
        double zero_point = this->kappa_ / gammaSq;        
        assert( zero_point >= 0.0 ); // Но это не точно
        this->weights_mean_( 0 ) = zero_point;
        this->weights_covariance_( 0 ) = zero_point;

        // Set the weights for other sigma points
        double other_points = 1.0 / ( 2.0 * gammaSq );
        assert( other_points > 0.0 );
        for( int i = 1; i < this->k_sigma_points_; i++ ) {
            this->weights_mean_( i ) = other_points;
            this->weights_covariance_( i ) = other_points;
        }

        // Параметры alpha_, beta_, lambda_ при данном способе создания сигма точек не используются:
        this->alpha_ = 0.0;
        this->beta_ = 0.0;
        this->lambda_ = 0.0;

#ifdef DEBUG_KALMAN
        std::cout << this->filterName_ + " SetupDesignParametersMeanSet:" << std::endl;
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

    ///
    /// \brief Установка параметров масштабируемого ансцентного преобразования (Scaled UT)
    /// \details Смотри [2], [5], [6]
    /// \param alpha - параметр разброса сигма-точек (alpha = 10^-3 - типичная рекомендация по van der Merwe)
    /// \param beta - параметр, отвечающий за характер распредеелния (beta = 2 - нормальное)
    /// \param kappa - параметр, отвечающий за разброс сигма-точек (kappa = 0 или (3 - SizeX) - типичная рекомендация по van der Merwe)
    ///
    virtual void SetupDesignParametersScaledSet( double alpha, double beta, double kappa )
    {
        assert( alpha > 0.0 && alpha <= 1.0 );
        assert( beta >= 0.0 );

        this->alpha_ = alpha;
        this->beta_ = beta;
        this->kappa_ = kappa;        
        double gammaSq = ( alpha * alpha ) * ( static_cast<double>( SizeX ) + kappa );
        assert( gammaSq > 0.0 );
        this->gamma_ = std::sqrt( gammaSq );
        this->lambda_ = gammaSq - static_cast<double>( SizeX );

        // Set the weights for zero sigma point
        this->weights_mean_( 0 ) = this->lambda_ / gammaSq;
        this->weights_covariance_( 0 ) = this->weights_mean_( 0 ) + ( 1.0 - ( alpha * alpha ) + beta );
//        assert( this->weights_mean_( 0 ) >= 0.0 );
        assert( this->weights_covariance_( 0 ) >= 0.0 );

        // Set the weights for other sigma points
        double other_points = 1.0 / ( 2.0 * gammaSq );
        assert( other_points > 0.0 );
        for( int i = 1; i < this->k_sigma_points_; i++ ) {
            this->weights_mean_( i ) = other_points;
            this->weights_covariance_( i ) = other_points;
        }

        // Параметр w0_ при данном способе создания сигма точек не используется:
        this->w0_ = 0.0;

#ifdef DEBUG_KALMAN
        std::cout << this->filterName_ + " SetupDesignParametersScaledSet:" << std::endl;
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

    ///
    /// \brief Установка параметра фильтра по рекомендации Central Difference Kalman Filter (CDKF)
    /// \details Смотри [4]
    /// \param h2 - параметр разброса сигма-точек (h^2 = 3 типичная рекомендация для гауссовых шумов, [3])
    ///
    virtual void SetupDesignParametersCDKF( double h2 )
    {
        assert( h2 > 0.0 );
        double gammaSq = h2;
        this->gamma_ = std::sqrt( gammaSq );

        // Set the weights for zero sigma point
        double zero_point = ( gammaSq - SizeX ) / gammaSq;
        assert( zero_point >= 0.0 );
        this->weights_mean_( 0 ) = zero_point;
        this->weights_covariance_( 0 ) = zero_point;

        // Set the weights for other sigma points
        double other_points = 1.0 / ( 2.0 * gammaSq );
        assert( other_points > 0.0 );
        for( int i = 1; i < this->k_sigma_points_; i++ ) {
            this->weights_mean_( i ) = other_points;
            this->weights_covariance_( i ) = other_points;
        }

        // Параметры alpha_, beta_, kappa_, w0_, lambda_ при данном способе создания сигма точек не используются:
        this->alpha_ = 0.0;
        this->beta_ = 0.0;
        this->kappa_ = 0.0;
        this->w0_ = 0.0;
        this->lambda_ = 0.0;

#ifdef DEBUG_KALMAN
        std::cout << this->filterName_ + " SetupDesignParametersCDKF:" << std::endl;
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
        this->PredictionUKF( dt );
    }

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
        this->CorrectionUKF( Y_msd );
    }

protected:
    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Прогноз UKF
    /// \param dt - Время прогноза, [с]
    ///
    void PredictionUKF( double dt )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Prediction started, dt = " << dt << std::endl;
#endif
        assert( dt > 0.0 ); // Прогноз имеет смысл только при положительном времени

        // 1. Создание сигма-точек пространства X
        this->sqrt_P_chol_ = arma::chol( this->P_, "lower" ); // Должно быть взято НИЖНЕЕ РАЗЛОЖЕНИЕ!
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Prediction, P before:" );
        ( this->sqrt_P_chol_ ).print( this->filterName_ + " Prediction, sqrt_P_chol:" );
#endif
        this->x_est_sigma_points_.col(0) = this->X_est_; // Нулевая сигма-точка - это вектор состояния
        for( size_t i = 1; i < ( SizeX + 1 ); i++ ) {
            arma::vec add = ( this->gamma_ * this->sqrt_P_chol_.col(i - 1) );
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
        this->P_ = ( this->Q_ * std::abs( std::abs( dt ) ) );
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = this->x_pred_sigma_points_.col(i) - this->X_pred_;
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->P_ += this->weights_covariance_(i) * this->dXcal_.col(i) * arma::trans( this->dXcal_.col(i) );
        }
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->dXcal_ ).print( this->filterName_ + " Prediction, dXcal_:" );
        ( this->P_ ).print( this->filterName_ + " Prediction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );
        // 5. Пересоздание сигма-точек после прогноза по новой P
        this->sqrt_P_chol_ = arma::chol( this->P_, "lower" ); // Должно быть взято НИЖНЕЕ РАЗЛОЖЕНИЕ!
        this->x_pred_sigma_points_.col(0) = this->X_pred_; // Нулевая сигма-точка - это вектор состояния
        for( size_t i = 1; i < ( SizeX + 1 ); i++ ) {
            arma::vec add = ( this->gamma_ * this->sqrt_P_chol_.col(i - 1) );
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
        // 6. Вычисление X_pred заново по пересозданным сигма-точкам пространства Х (в [4] не указано, но это подразумевается)
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
    /// \brief Коррекция UKF
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    void CorrectionUKF( const arma::vec &Y_msd )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Correction started" << std::endl;
        Y_msd.print( this->filterName_ + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла:
        this->prediction_isDone = false; // Сразу же снять признак

        // 1. Вычисление ковариационной матрицы S (P_yy)
        this->S_ = this->R_;
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dYcal_.col(i) = this->y_pred_sigma_points_.col( i ) - this->Y_pred_;
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->dYcal_.col(i) = this->checkDeltaMeasurement_( this->dYcal_.col(i) );
            }
            this->S_ += this->weights_covariance_(i) * this->dYcal_.col(i) * arma::trans( this->dYcal_.col(i) );
        }
        this->fixMatrixMainDiagonalSymmetry( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( this->filterName_ + " Correction, S:" );
#endif
        this->checkMatrixDiagPositive( this->S_ );
        // 2. Вычисление кросс-ковариационной матрицы P_xy
        this->P_xy_.zeros();
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = this->x_pred_sigma_points_.col(i) - this->X_pred_;
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            // dYcal_ - вычислено выше в 1, поэтому заново - не надо
            this->P_xy_ += this->weights_covariance_(i) * this->dXcal_.col(i) * arma::trans( this->dYcal_.col(i) );
        }
#ifdef DEBUG_KALMAN
        ( this->dXcal_ ).print( this->filterName_ + " Prediction, dXcal_:" );
        ( this->dYcal_ ).print( this->filterName_ + " Prediction, dYcal_:" );
        ( this->P_xy_ ).print( this->filterName_ + " Correction, P_xy:" );
#endif
        // 3. Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_xy_ * arma::inv( this->S_ );
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
        this->P_ = this->P_ - ( this->K_ * this->S_ * arma::trans( this->K_ ) );
//        this->P_ = this->P_ - P_xy * arma::trans( this->K_ );
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Correction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );
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

    //------------------------------------------------------------------------------------------------------------------
    // Параметры, зависящие от SizeX, SizeY:
    static const int k_sigma_points_ = ( 2 * SizeX ) + 1; ///< Число сигма-точек

    arma::vec::fixed<k_sigma_points_> weights_mean_; ///< Веса среднего
    arma::vec::fixed<k_sigma_points_> weights_covariance_; ///< Веса ковариации
    arma::mat::fixed<SizeX, k_sigma_points_> x_est_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве X на текущем такте, размерность [SizeX,k_sigma_points_]
    arma::mat::fixed<SizeX, k_sigma_points_> x_pred_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве X, экстраполированный на текущий такт, размерность [SizeX,k_sigma_points_]
    arma::mat::fixed<SizeY, k_sigma_points_> y_pred_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве Y, экстраполированный на текущий такт, размерность [SizeX,k_sigma_points_]

    arma::mat::fixed<SizeX, k_sigma_points_> dXcal_; ///< Матрица Х-каллиграфическое (матрица сигма-точек - столбцов)
    arma::mat::fixed<SizeY, k_sigma_points_> dYcal_; ///< Матрица Y-каллиграфическое (матрица сигма-точек - столбцов)
    arma::mat::fixed<SizeX, SizeY> P_xy_; ///< Матрица кросс-коварации векторов Х и Y, размерность [SizeX * SizeY]
    arma::mat::fixed<SizeX, SizeX> sqrt_P_chol_; ///< Корень из матрицы P

    //------------------------------------------------------------------------------------------------------------------
    // Параметры выбора сигма-точек, определяемые в методах SetDesignParameters* :
    //
    // Параметры сигма-точек согласно van der Merwe:
    double alpha_;  ///< Параметр разброса сигма-точек (alpha = 10^-3 типичная рекомендация)
    double kappa_;  ///< Параметр разброса сигма-точек (kappa = 3 - SizeX типичная рекомендация)
    double beta_;   ///< Параметр разброса сигма-точек (beta = 2 - нормальное, 0 - нет сведений о распределении)
    double lambda_; ///< Автоматически вычисляемый параметр, равный ( alpha * alpha ) * ( SizeX + kappa ) - SizeX;

    // Параметры сигма-точек согласно Julier, Uhlmann:
    double w0_;     ///< Параметр разброса сигма-точек (0..1)

    double gamma_;  ///< Автоматически вычисляемый (в методах SetDesignParameters*) параметр (множитель при корне из P при создании сигма-точек)
    //------------------------------------------------------------------------------------------------------------------
    // Обертки функций вычисления взвешенной суммы:
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumStateSigmas_; ///< Вычисление взвешенной суммы сигма-точек пространства Х
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumMeasurementSigmas_; ///< Вычисление взвешенной суммы сигма-точек пространства Y

    //------------------------------------------------------------------------------------------------------------------
    // Запрет доступа к методам родительского класса:
    void SetStateTransitionJacobianF(){} ///< Запрет доступа
    void SetObservationJacobianH(){} ///< Запрет доступа
};
}
#endif // KALMAN_FILTER_UNSCENTED_H
/// \}
