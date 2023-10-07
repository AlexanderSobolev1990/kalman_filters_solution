//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_extended_cubature.h
/// \brief      Шаблонный класс расширенного кубатурного фильтра Калмана, РКФК (Extended Cubature Kalman Filter, ECKF)
/// \date       30.03.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_CUBATURE_H
#define KALMAN_FILTER_EXTENDED_CUBATURE_H

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_cubature.h>

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс расширенного кубатурного фильтра Калмана, РКФК (Extended Cubature Kalman Filter, ECKF)
/// \details Источники:
/// \n[1] Cubature Kalman Filters, Ienkaran Arasaratnam and Simon Haykin, Life Fellow, IEEE
/// \attention Фильтр построен по классическому его варианту НИЖНЕтреугольного разложения Холецкого!
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanECKF : public CKalmanCKF<SizeX, SizeY>
{
public:
    using CKalmanCKF<SizeX, SizeY>::CKalmanCKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanECKF() : CKalmanCKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "ECKF" );
#endif
    }
    // default copy/move/assignment semantic:
//    CKalmanECKF( const CKalmanECKF& ) = default;
//    CKalmanECKF& operator=( const CKalmanECKF& ) = default;
//    CKalmanECKF( CKalmanECKF&& ) = default;
//    CKalmanECKF& operator=( CKalmanECKF&& ) = default;
//    virtual ~CKalmanECKF() = default;

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
    /// \brief Прогноз расширенного кубатурного фильтра Калмана (РКФК, ECKF)
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
        this->PredictionEKF( dt );
    }

    ///
    /// \brief Коррекция расширенного кубатурного фильтра Калмана (РКФК, ECKF)
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
        // 1a. Создание сигма-точек пространства X
        this->sqrt_P_chol_ = arma::chol( this->P_, "lower" ); // Должно быть взято НИЖНЕЕ РАЗЛОЖЕНИЕ!
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->sqrt_P_chol_.col(i) );
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

        this->CorrectionCKF( Y_msd );
    }
};

}

#endif // KALMAN_FILTER_EXTENDED_CUBATURE_H
/// \}
