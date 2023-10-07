//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_extended_unscented.h
/// \brief      Шаблонный класс расширенного сигма-точечного (ансцентного) фильтра Калмана, РСТФК
///             (Extended Unscented Kalman Filter, EUKF)
/// \date       11.04.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_UNSCENTED_H
#define KALMAN_FILTER_EXTENDED_UNSCENTED_H

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_unscented.h>

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief  Шаблонный класс расширенного сигма-точечного (ансцентного) фильтра Калмана, РСТФК
///         (Extended Unscented Kalman Filter, EUKF)
/// \details Источники:
/// \n[1] A New Extension of the Kalman Filter to Nonlinear Systems, Simon J. Julier, Jeffrey K. Uhlmann,
/// The Robotics Research Group, Department of Engineering Science, The University of Oxford, 1997
/// \n[2] The Unscented Kalman Filter for Nonlinear Estimation, Eric A. Wan and Rudolph van der Merwe
/// Oregon Graduate Institute of Science & Technology 20000 NW Walker Rd, Beaverton, Oregon 97006, 2000
/// ericwan[at]ece.ogi.edu, rvdmerwe[at]ece.ogi.edu
/// \n[3] Sigma-Point Kalman Filters for Probabilistic Inference in Dynamic State-Space Models, Rudolph van der Merwe &
/// Eric Wan, OGI School of Science & Engineering Oregon Health & Science University Beaverton, Oregon, 97006, USA, 2003
/// {rvdmerwe,ericwan}[at]ece.ogi.edu
/// \n[4] THE SQUARE-ROOT UNSCENTED KALMAN FILTER FOR STATE AND PARAMETER-ESTIMATION, Rudolph van der Merwe and
/// Eric A. Wan, Oregon Graduate Institute of Science and Technology 20000 NW Walker Road, Beaverton, Oregon 97006, USA
/// rvdmerwe,ericwan[at]ece.ogi.edu
/// \attention Внимание! В реализации UKF вес нулевой сигма-точки Wcov не может быть отрицательным! (Wmean - может)
/// \attention Фильтр построен по классическому его варианту НИЖНЕтреугольного разложения Холецкого!
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanEUKF : public CKalmanUKF<SizeX, SizeY>
{
public:
    using CKalmanUKF<SizeX, SizeY>::CKalmanUKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanEUKF() : CKalmanUKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "EUKF" );
#endif
    }
    // default copy/move/assignment semantic:
    CKalmanEUKF( const CKalmanEUKF& ) = default;
    CKalmanEUKF& operator=( const CKalmanEUKF& ) = default;
    CKalmanEUKF( CKalmanEUKF&& ) = default;
    CKalmanEUKF& operator=( CKalmanEUKF&& ) = default;
    virtual ~CKalmanEUKF() = default;

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
    /// \brief Прогноз расширенного сигма-точечного фильтра Калмана (РСТФК, EUKF)
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
        this->PredictionEKF( dt );
    }

    ///
    /// \brief Коррекция расширенного сигма-точечного фильтра Калмана (РСТФК, EUKF)
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
        // 1a. Создание сигма-точек пространства X
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
        // 1b. Вычисление сигма-точек пространства Y
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }

        this->CorrectionUKF( Y_msd );
    }
};

}

#endif // KALMAN_FILTER_EXTENDED_UNSCENTED_H
/// \}
