//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_unscented_square_root.h
/// \brief      Шаблонный класс квадратно-корневого сигма-точечного (ансцентного) фильтра Калмана, КК-СТФК (КК-АФК)
/// (Square Root Unscented Kalman Filter, SR-UKF)
/// \date       11.01.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_UNSCENTED_SQUARE_ROOT_H
#define KALMAN_FILTER_UNSCENTED_SQUARE_ROOT_H

//#ifdef DEBUG_KALMAN
//#undef DEBUG_KALMAN
//#endif
//#define DEBUG_KALMAN // Включение отладочной печати в консоль

// System includes:
#include <float.h>

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_unscented.h>
#include <qr_decomposition.h> // JQR разложение матрицы

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс квадратно-корневого сигма-точечного (ансцентного) фильтра Калмана, КК-СТФК (КК-АФК)
/// (Square Root Unscented Kalman Filter, SR-UKF)
/// \details Источники:
/// \n[1] A New Extension of the Kalman Filter to Nonlinear Systems, Simon J. Julier, Jeffrey K. Uhlmann,
/// The Robotics Research Group, Department of Engineering Science, The University of Oxford, 1997
/// \n[2] The Unscented Kalman Filter for Nonlinear Estimation, Eric A. Wan and Rudolph van der Merwe
/// Oregon Graduate Institute of Science & Technology 20000 NW Walker Rd, Beaverton, Oregon 97006, 2000
/// ericwan@ece.ogi.edu, rvdmerwe@ece.ogi.edu
/// \n[3] Sigma-Point Kalman Filters for Probabilistic Inference in Dynamic State-Space Models, Rudolph van der Merwe &
/// Eric Wan, OGI School of Science & Engineering Oregon Health & Science University Beaverton, Oregon, 97006, USA, 2003
/// {rvdmerwe,ericwan}[at]ece.ogi.edu
/// \n[4] THE SQUARE-ROOT UNSCENTED KALMAN FILTER FOR STATE AND PARAMETER-ESTIMATION, Rudolph van der Merwe and
/// Eric A. Wan, Oregon Graduate Institute of Science and Technology 20000 NW Walker Road, Beaverton, Oregon 97006, USA
/// rvdmerwe,ericwan [at]ece.ogi.edu
/// \n[5] The J-Orthogonal Square-Root Euler-Maruyama-Based Unscented Kalman Filter for Nonliear Stohastic Systems,
/// Gennady Yu. Kulikov, Maria V. Kulikova, CEMAT, Instituto Superior Técnico, Universidade de Lisboa, Av.
/// Rovisco Pais 1, 1049-001 LISBOA, Portugal (emails: gennady.kulikov[at]tecnico.ulisboa.pt, maria.kulikova[at]ist.utl.pt)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRUKF : public CKalmanUKF<SizeX, SizeY>
{
public:
    using CKalmanUKF<SizeX, SizeY>::CKalmanUKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSRUKF() : CKalmanUKF<SizeX, SizeY>(),
        negativeZeroCovWeight_( false )
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SRUKF" );
#endif
        createSignMatrices();
    }
    // default copy/move/assignment semantic:
    CKalmanSRUKF( const CKalmanSRUKF& ) = default;
    CKalmanSRUKF& operator=( const CKalmanSRUKF& ) = default;
    CKalmanSRUKF( CKalmanSRUKF&& ) = default;
    CKalmanSRUKF& operator=( CKalmanSRUKF&& ) = default;
    virtual ~CKalmanSRUKF() = default;

    //------------------------------------------------------------------------------------------------------------------
    // Методы-сеттеры:

    ///
    /// \brief Установка параметра w0 ансцентного фильтра (MeanSet)
    /// \details При выборе w0 = [0...1) обеспечиваются положительные веса. Смотри [1], [2]
    /// \attention Нельзя выбирать w0 так, чтобы нулевой вес был > 0, а остальные меньше нуля. Наоборот - МОЖНО, т.е.
    /// нулевой вес может быть отрицаительным.
    /// \param w0 - параметр разброса сигма точек (w0 = [0...1) типичная рекомендация для положительных весов)
    ///
    virtual void SetupDesignParametersMeanSet( double w0 )
    {
        this->w0_ = w0;
        this->kappa_ = ( w0 * static_cast<double>( SizeX ) ) / ( 1.0 - w0 );
        double gammaSq = static_cast<double>( SizeX ) + this->kappa_;
        assert( gammaSq > 0.0 );
        this->gamma_ = std::sqrt( gammaSq );

        // В SR_UKF нулевую точку положим "назад":

        // Set the weights for zero sigma point:
        double zero_point = this->kappa_ / gammaSq;
        if( zero_point < 0.0 ) {
            this->negativeZeroCovWeight_ = true;
        } else {
            this->negativeZeroCovWeight_ = false;
        }
        this->weights_mean_( this->k_sigma_points_ - 1 ) = zero_point;
        this->weights_covariance_( this->k_sigma_points_ - 1 ) = zero_point;

        // Set the weights for other sigma points
        double other_points = 1.0 / ( 2.0 * gammaSq );
        for( int i = 0; i < ( this->k_sigma_points_ - 1 ); i++ ) {
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

        // В SR_UKF нулевую точку положим "назад":

        // Set the weights for zero sigma point
        this->weights_mean_( this->k_sigma_points_ - 1 ) = this->lambda_ / gammaSq;
        this->weights_covariance_( this->k_sigma_points_ - 1 ) =
            this->weights_mean_( this->k_sigma_points_ - 1 ) + ( 1.0 - ( alpha * alpha ) + beta );
        if( this->weights_covariance_( this->k_sigma_points_ - 1 ) < 0.0 ) {
            this->negativeZeroCovWeight_ = true;
        } else {
            this->negativeZeroCovWeight_ = false;
        }
        // Set the weights for other sigma points
        double other_points = 1.0 / ( 2.0 * gammaSq );
        for( int i = 0; i < ( this->k_sigma_points_ - 1 ); i++ ) {
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
    /// \details Смотри [3]
    /// \param h2 - параметр разброса сигма-точек (h^2 = 3 типичная рекомендация для гауссовых шумов, [3])
    ///
    virtual void SetupDesignParametersCDKF( double h2 )
    {
        assert( h2 > 0.0 );
        double gammaSq = h2;
        this->gamma_ = std::sqrt( gammaSq );

        // В SR_UKF нулевую точку положим "назад":

        // Set the weights for zero sigma point
        double zero_point = ( gammaSq - SizeX ) / gammaSq;
        if( zero_point < 0.0 ) {
            this->negativeZeroCovWeight_ = true;
        } else {
            this->negativeZeroCovWeight_ = false;
        }
        this->weights_mean_( this->k_sigma_points_ - 1 ) = zero_point;
        this->weights_covariance_( this->k_sigma_points_ - 1 ) = zero_point;

        // Set the weights for other sigma points
        double other_points = 1.0 / ( 2.0 * gammaSq );
        assert( other_points > 0.0 );
        for( int i = 0; i < ( this->k_sigma_points_ - 1 ); i++ ) {
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
        this->x_est_sigma_points_.col( this->k_sigma_points_ - 1 ) = this->X_est_; // Нулевая точка - сзади!
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
            this->dXcal_.col(i) *= std::sqrt( std::abs( this->weights_covariance_(i) ) );
        }
        arma::mat Qdt = this->Q_ * std::sqrt( std::abs( dt ) );
        arma::mat JQR_input_P_pred = arma::trans( arma::join_horiz( Qdt, this->dXcal_ ) ); // JQR_input = [ Qdt, dXcal ], Транспонировать, т.к. в JQR разложение так надо
#ifdef DEBUG_KALMAN
        JQR_input_P_pred.print( this->filterName_ + " Prediction, JQR_input_P_pred:" );
#endif
        arma::mat Q_qr_P_pred, R_qr_P_pred, J_qr_P_pred;
        int res_JQR_P_pred;
        if( this->negativeZeroCovWeight_ ) {
            res_JQR_P_pred = SPML::QR::J_orthogonal( Q_qr_P_pred, R_qr_P_pred, J_qr_P_pred, JQR_input_P_pred, this->Jpredict_, true ); // Вызов JQR разложения
        } else {
            res_JQR_P_pred = SPML::QR::ModifiedGramSchmidt( Q_qr_P_pred, R_qr_P_pred, JQR_input_P_pred ); // Вызов QR разложения
        }
        assert( res_JQR_P_pred == 0 );
        this->P_ = arma::trans( R_qr_P_pred ); // Матрица Р считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        J_qr_P_pred.print( this->filterName_ + " Prediction, J_qr_P_pred after:" );
        Q_qr_P_pred.print( this->filterName_ + " Prediction, Q_qr_P_pred after:" );
        ( this->P_ ).print( this->filterName_ + " Prediction, P after:" );
        ( this->P_ * arma::trans( this->P_ ) ).print( this->filterName_ + " Prediction, Pfull after:" );
#endif
        // 5. Пересоздание сигма-точек после прогноза по новой P
        this->x_pred_sigma_points_.col(this->k_sigma_points_ - 1) = this->X_pred_; // Нулевая точка - сзади!
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
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
        this->CorrectionSRUKF( Y_msd );
    }

protected:
    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Коррекция SRUKF
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    void CorrectionSRUKF( const arma::vec &Y_msd )
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
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dYcal_.col(i) = ( this->y_pred_sigma_points_.col(i) - this->Y_pred_ );
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->dYcal_.col(i) = this->checkDeltaMeasurement_( this->dYcal_.col(i) );
            }
            this->dYcal_.col(i) *= std::sqrt( std::abs( this->weights_covariance_(i) ) );
        }
        arma::mat JQR_input_S_corr = arma::trans( arma::join_horiz( this->R_, this->dYcal_ ) ); // JQR_input_S_corr = [ R, dYcal ], Транспонировать, т.к. в JQR разложение так надо
#ifdef DEBUG_KALMAN
        JQR_input_S_corr.print( this->filterName_ + " Correction, JQR_input_S_corr:" );
#endif
        arma::mat Q_qr_S_corr, R_qr_S_corr, J_qr_S_corr;
        int res_JQR_S_corr;
        if( this->negativeZeroCovWeight_ ) {
            res_JQR_S_corr = SPML::QR::J_orthogonal( Q_qr_S_corr, R_qr_S_corr, J_qr_S_corr, JQR_input_S_corr, this->Jcorrect_, true ); // Вызов JQR разложения
        } else {
            res_JQR_S_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_S_corr, R_qr_S_corr, JQR_input_S_corr ); // Вызов QR разложения
        }
        assert( res_JQR_S_corr == 0 );
        this->S_ = arma::trans( R_qr_S_corr ); // Матрица S считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        J_qr_S_corr.print( this->filterName_ + " Correction, J_qr_S_corr after:" );
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
            this->dXcal_.col(i) *= std::sqrt( std::abs( this->weights_covariance_(i) ) );
        }
        this->P_xy_ = this->dXcal_ * this->J_ * arma::trans( this->dYcal_ );
#ifdef DEBUG_KALMAN
        ( this->P_xy_ ).print( this->filterName_ + " Correction, P_xy:" );
#endif
        // 3. Вычисление коэффициента усиления фильтра K

        // 1 способ
        this->K_ = ( this->P_xy_ * arma::inv( arma::trans( this->S_ ) ) ) * arma::inv( this->S_ );

        // 2 способ
//        this->K_ = arma::trans( arma::solve( arma::trans( this->S_ ), arma::solve( this->S_, arma::trans( P_xy ) ) ) );

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
        arma::mat JQR_input_P_corr = arma::trans( arma::join_horiz( P1, P2 ) ); // Транспонировать, т.к. в JQR разложение так надо
#ifdef DEBUG_KALMAN
        JQR_input_P_corr.print( this->filterName_ + " Correction, JQR_input_P_corr:" );
#endif
        arma::mat Q_qr_P_corr, R_qr_P_corr, J_qr_P_corr;
        int res_JQR_P_corr;
        if( this->negativeZeroCovWeight_ ) {
            res_JQR_P_corr = SPML::QR::J_orthogonal( Q_qr_P_corr, R_qr_P_corr, J_qr_P_corr, JQR_input_P_corr, this->Jcorrect_, true ); // Вызов JQR разложения
        } else {
            res_JQR_P_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_P_corr, R_qr_P_corr, JQR_input_P_corr ); // Вызов QR разложения
        }
        assert( res_JQR_P_corr == 0 );
        this->P_ = arma::trans( R_qr_P_corr ); // Матрица P считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        J_qr_P_corr.print( this->filterName_ + " Correction, J_qr_P_corr after:" );
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

    //------------------------------------------------------------------------------------------------------------------
    // Матрицы знаков

    static const int JQR_predict_size = ( 2 * SizeX + 1 ) + SizeX;
    static const int JQR_correct_size = ( 2 * SizeX + 1 ) + SizeY;

    arma::mat::fixed<CKalmanSRUKF::k_sigma_points_, CKalmanSRUKF::k_sigma_points_> J_ =
        arma::mat::fixed<CKalmanSRUKF::k_sigma_points_, CKalmanSRUKF::k_sigma_points_>( arma::fill::eye ); ///< Матрица знаков при Pxy
    arma::mat::fixed<JQR_predict_size, JQR_predict_size> Jpredict_ =
        arma::mat::fixed<JQR_predict_size, JQR_predict_size>( arma::fill::eye ); ///< Матрица знаков при прогнозе
    arma::mat::fixed<JQR_correct_size, JQR_correct_size> Jcorrect_ =
        arma::mat::fixed<JQR_correct_size, JQR_correct_size>( arma::fill::eye ); ///< Матрица знаков при коррекции

    bool negativeZeroCovWeight_; ///< Признак отрицательного веса "нулевой" сигма-точки Wcov
    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Создание матриц знаков
    ///
    void createSignMatrices()
    {
        if( this->weights_covariance_( this->k_sigma_points_ - 1 ) < 0.0 ) {
            this->Jpredict_( JQR_predict_size - 1, JQR_predict_size - 1 ) = -1.0;
            this->Jcorrect_( JQR_correct_size - 1, JQR_correct_size - 1 ) = -1.0;
        }

        if( this->weights_covariance_( this->k_sigma_points_ - 1 ) < 0.0 ) {
            this->J_( this->k_sigma_points_ - 1, this->k_sigma_points_ - 1 ) = -1.0;
        }
    }
};

//----------------------------------------------------------------------------------------------------------------------
// Добавлен для сравнения с CKalmanSRUKF, подразумевающим пересоздание сигма-точек.
template<size_t SizeX, size_t SizeY>
class CKalmanSRUKF2 : public CKalmanSRUKF<SizeX, SizeY>
{
public:
//    using CKalmanUKF<SizeX, SizeY>::CKalmanUKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSRUKF2() : CKalmanSRUKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SRUKF2" );
#endif
    }
    // default copy/move/assignment semantic:
    CKalmanSRUKF2( const CKalmanSRUKF2& ) = default;
    CKalmanSRUKF2& operator=( const CKalmanSRUKF2& ) = default;
    CKalmanSRUKF2( CKalmanSRUKF2&& ) = default;
    CKalmanSRUKF2& operator=( CKalmanSRUKF2&& ) = default;
    virtual ~CKalmanSRUKF2() = default;

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
        this->x_est_sigma_points_.col( this->k_sigma_points_ - 1 ) = this->X_est_; // Нулевая точка - сзади!
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
            this->dXcal_.col(i) *= std::sqrt( std::abs( this->weights_covariance_(i) ) );
        }
        arma::mat Qdt = this->Q_ * std::sqrt( std::abs( dt ) );
        arma::mat JQR_input_P_pred = arma::trans( arma::join_horiz( Qdt, this->dXcal_ ) ); // JQR_input = [ Qdt, dXcal ], Транспонировать, т.к. в JQR разложение так надо
#ifdef DEBUG_KALMAN
        JQR_input_P_pred.print( this->filterName_ + " Prediction, JQR_input_P_pred:" );
#endif
        arma::mat Q_qr_P_pred, R_qr_P_pred, J_qr_P_pred;
        int res_JQR_P_pred;
        if( this->negativeZeroCovWeight_ ) {
            res_JQR_P_pred = SPML::QR::J_orthogonal( Q_qr_P_pred, R_qr_P_pred, J_qr_P_pred, JQR_input_P_pred, this->Jpredict_, true ); // Вызов JQR разложения
        } else {
            res_JQR_P_pred = SPML::QR::ModifiedGramSchmidt( Q_qr_P_pred, R_qr_P_pred, JQR_input_P_pred ); // Вызов QR разложения
        }
        assert( res_JQR_P_pred == 0 );
        this->P_ = arma::trans( R_qr_P_pred ); // Матрица Р считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        J_qr_P_pred.print( this->filterName_ + " Prediction, J_qr_P_pred after:" );
        Q_qr_P_pred.print( this->filterName_ + " Prediction, Q_qr_P_pred after:" );
        ( this->P_ ).print( this->filterName_ + " Prediction, P after:" );
        ( this->P_ * arma::trans( this->P_ ) ).print( this->filterName_ + " Prediction, Pfull after:" );
#endif
//        // 5. Пересоздание сигма-точек после прогноза по новой P
//        this->x_pred_sigma_points_.col(this->k_sigma_points_ - 1) = this->X_pred_; // Нулевая точка - сзади!
//        for( size_t i = 0; i < SizeX; i++ ) {
//            arma::vec add = this->gamma_ * this->P_.col(i);
//            this->x_pred_sigma_points_.col(i) = this->X_pred_ + add;
//            this->x_pred_sigma_points_.col(i + SizeX) = this->X_pred_ - add;
//            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
//                this->x_pred_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i) );
//                this->x_pred_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i + SizeX) );
//            }
//        }
#ifdef DEBUG_KALMAN
        ( this->x_pred_sigma_points_ ).print( this->filterName_ + " Prediction, x_pred_sigma_points_:" );
#endif
//        // 6. Вычисление X_pred заново по пересозданным сигма-точкам пространства Х (в [4] не указано, но это подразумевается)
//        if( this->weightedSumStateSigmas_ == nullptr ) {
//            this->X_pred_ = this->x_pred_sigma_points_ * this->weights_mean_; // В матричной форме
//        } else {
//            this->X_pred_ = this->weightedSumStateSigmas_( this->weights_mean_, this->x_pred_sigma_points_ );
//        }
//        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
//            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
//        }

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
};

//----------------------------------------------------------------------------------------------------------------------
// Добавлен для сравнения с CKalmanSRUKF, подразумевающим пересоздание сигма-точек.
template<size_t SizeX, size_t SizeY>
class CKalmanSRUKF3 : public CKalmanSRUKF<SizeX, SizeY>
{
public:
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSRUKF3() : CKalmanSRUKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SRUKF3" );
#endif
//        createSignMatrices2();
    }
    // default copy/move/assignment semantic:
    CKalmanSRUKF3( const CKalmanSRUKF3& ) = default;
    CKalmanSRUKF3& operator=( const CKalmanSRUKF3& ) = default;
    CKalmanSRUKF3( CKalmanSRUKF3&& ) = default;
    CKalmanSRUKF3& operator=( CKalmanSRUKF3&& ) = default;
    virtual ~CKalmanSRUKF3() = default;

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
        // Сразу увеличим шумы
        arma::mat Qdt = this->Q_ * std::sqrt( std::abs( dt ) );
        this->P_ = this->P_ + Qdt;

        this->x_est_sigma_points_.col( this->k_sigma_points_ - 1 ) = this->X_est_; // Нулевая точка - сзади!
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
            this->dXcal_.col(i) *= std::sqrt( std::abs( this->weights_covariance_(i) ) );
        }
//        arma::mat Qdt = this->Q_ * std::sqrt( dt );
        arma::mat Qdt2 = this->Q_ * 1.0e-16;//std::sqrt( dt );
        arma::mat JQR_input_P_pred = arma::trans( arma::join_horiz( Qdt2, this->dXcal_ ) ); // JQR_input = [ Qdt, dXcal ], Транспонировать, т.к. в JQR разложение так надо
//        arma::mat JQR_input_P_pred = arma::trans( this->dXcal_ );
#ifdef DEBUG_KALMAN
        JQR_input_P_pred.print( this->filterName_ + " Prediction, JQR_input_P_pred:" );
#endif
        arma::mat Q_qr_P_pred, R_qr_P_pred, J_qr_P_pred;
        int res_JQR_P_pred;
        if( this->negativeZeroCovWeight_ ) {
            res_JQR_P_pred = SPML::QR::J_orthogonal( Q_qr_P_pred, R_qr_P_pred, J_qr_P_pred, JQR_input_P_pred, this->Jpredict_, true ); // Вызов JQR разложения
        } else {
            res_JQR_P_pred = SPML::QR::ModifiedGramSchmidt( Q_qr_P_pred, R_qr_P_pred, JQR_input_P_pred ); // Вызов QR разложения
        }
        assert( res_JQR_P_pred == 0 );
        this->P_ = arma::trans( R_qr_P_pred ); // Матрица Р считывается из транспонированной выходной матрицы R
#ifdef DEBUG_KALMAN
        J_qr_P_pred.print( this->filterName_ + " Prediction, J_qr_P_pred after:" );
        Q_qr_P_pred.print( this->filterName_ + " Prediction, Q_qr_P_pred after:" );
        ( this->P_ ).print( this->filterName_ + " Prediction, P after:" );
        ( this->P_ * arma::trans( this->P_ ) ).print( this->filterName_ + " Prediction, Pfull after:" );
#endif
//        // 5. Пересоздание сигма-точек после прогноза по новой P
//        this->x_pred_sigma_points_.col(this->k_sigma_points_ - 1) = this->X_pred_; // Нулевая точка - сзади!
//        for( size_t i = 0; i < SizeX; i++ ) {
//            arma::vec add = this->gamma_ * this->P_.col(i);
//            this->x_pred_sigma_points_.col(i) = this->X_pred_ + add;
//            this->x_pred_sigma_points_.col(i + SizeX) = this->X_pred_ - add;
//            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
//                this->x_pred_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i) );
//                this->x_pred_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i + SizeX) );
//            }
//        }
#ifdef DEBUG_KALMAN
        ( this->x_pred_sigma_points_ ).print( this->filterName_ + " Prediction, x_pred_sigma_points_:" );
#endif
//        // 6. Вычисление X_pred заново по пересозданным сигма-точкам пространства Х (в [4] не указано, но это подразумевается)
//        if( this->weightedSumStateSigmas_ == nullptr ) {
//            this->X_pred_ = this->x_pred_sigma_points_ * this->weights_mean_; // В матричной форме
//        } else {
//            this->X_pred_ = this->weightedSumStateSigmas_( this->weights_mean_, this->x_pred_sigma_points_ );
//        }
//        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
//            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
//        }

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

//    //------------------------------------------------------------------------------------------------------------------
//    ///
//    /// \brief Создание матриц знаков
//    ///
//    void createSignMatrices2()
//    {
//        int JQR_predict_size = ( 2 * SizeX + 1 );// + SizeX;
//        int JQR_correct_size = ( 2 * SizeX + 1 ) + SizeY;
//        this->Jpredict_ = arma::mat( JQR_predict_size, JQR_predict_size, arma::fill::eye );
//        this->Jcorrect_ = arma::mat( JQR_correct_size, JQR_correct_size, arma::fill::eye );
//        if( this->weights_covariance_( this->k_sigma_points_- 1 ) < 0.0 ) {
//            this->Jpredict_( JQR_predict_size-1, JQR_predict_size-1 ) = -1.0;
//            this->Jcorrect_( JQR_correct_size-1, JQR_correct_size-1 ) = -1.0;
//        }

//        this->J_ = arma::mat( this->k_sigma_points_, this->k_sigma_points_, arma::fill::eye );
//        if( this->weights_covariance_( this->k_sigma_points_- 1 ) < 0.0 ) {
//            this->J_( this->k_sigma_points_- 1, this->k_sigma_points_- 1 ) = -1.0;
//        }
//    }

};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс квадратно-корневого сигма-точечного (ансцентного) фильтра Калмана (блочная фильтрация),
/// КК-СТФКБ (КК-АФКБ) (Square Root Unscented Kalman Filter Block, SR-UKFB)
/// \details см. CKalmanSRUKF
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRUKFB : public CKalmanSRUKF<SizeX, SizeY>
{
public:
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSRUKFB() : CKalmanSRUKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SRUKFB" );
#endif
        createSignMatricesBlock();
    }
    // default copy/move/assignment semantic:
    CKalmanSRUKFB( const CKalmanSRUKFB& ) = default;
    CKalmanSRUKFB& operator=( const CKalmanSRUKFB& ) = default;
    CKalmanSRUKFB( CKalmanSRUKFB&& ) = default;
    CKalmanSRUKFB& operator=( CKalmanSRUKFB&& ) = default;
    virtual ~CKalmanSRUKFB() = default;

    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    // Prediction как в НЕблочном фильтре

    ///
    /// \brief Коррекция квадратно-корневого сигма-точечного фильтра Калмана (КК-СТФК, SR-UKF)
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
            this->dYcal_.col(i) *= std::sqrt( std::abs( this->weights_covariance_(i) ) );
        }
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = ( this->x_pred_sigma_points_.col(i) - this->X_pred_ );
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->dXcal_.col(i) *= std::sqrt( std::abs( this->weights_covariance_(i) ) );
        }

        size_t rowsB = ( SizeY + SizeX );
        size_t colsB = ( SizeY + ( 2 * SizeX + 1 ) );
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
        arma::mat Q_qr_corr, R_qr_corr, J_qr_corr;
        int res_JQR_corr;
        if( this->negativeZeroCovWeight_ ) {
            res_JQR_corr = SPML::QR::J_orthogonal( Q_qr_corr, R_qr_corr, J_qr_corr, B, this->JcorrectBlock_, true );
        } else {
            res_JQR_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_corr, R_qr_corr, B );
        }
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
        Q_qr_corr.print( this->filterName_ + " Correction, Q_qr_corr:" );
        R_qr_corr.print( this->filterName_ + " Correction, R_qr_corr:" );

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
    static const int QRsizeY = ( SizeY + ( 2 * SizeX + 1 ) );
    arma::mat::fixed<QRsizeY, QRsizeY> JcorrectBlock_ = arma::mat::fixed<QRsizeY, QRsizeY>( arma::fill::eye ); ///< Матрица знаков для фильтра в блочном виде

    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Создание матриц знаков
    ///
    void createSignMatricesBlock()
    {        
        if( this->weights_covariance_( this->k_sigma_points_- 1 ) < 0.0 ) {
            JcorrectBlock_( QRsizeY - 1, QRsizeY - 1 ) = -1.0;
        }
    }
};

}

#endif // KALMAN_FILTER_UNSCENTED_SQUARE_ROOT_H
/// \}
