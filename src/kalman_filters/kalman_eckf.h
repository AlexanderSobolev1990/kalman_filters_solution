//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_eckf.h
/// \brief      Шаблонный класс расширенно-кубатурного фильтра Калмана, РКФК 
///             (extended cubature Kalman filter, ECKF)
/// \date       30.03.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_CUBATURE_H
#define KALMAN_FILTER_EXTENDED_CUBATURE_H

// Project includes:
#include <kalman_debug.h>
#include <kalman_ckf.h>

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief Шаблонный класс расширенн-кубатурного фильтра Калмана, РКФК 
///        (extended cubature Kalman filter, ECKF)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanECKF : public CKalmanCKF<SizeX, SizeY>
{
public:    
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanECKF() = default;
    CKalmanECKF( const CKalmanECKF& ) = default;
    CKalmanECKF& operator=( const CKalmanECKF& ) = default;
    CKalmanECKF( CKalmanECKF&& ) = default;
    CKalmanECKF& operator=( CKalmanECKF&& ) = default;
    virtual ~CKalmanECKF() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "ECKF"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Cubature; }
    virtual TCovType GetCovType() const override { return TCovType::CT_FullMatrices; }    

    //--------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt ) override
    {
        this->predictionEKF( dt );
    }

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd ) override
    {
        this->drawSigmaPoints();
        this->correctionCKF( Y_msd );
    }

protected:
    //--------------------------------------------------------------------------
    ///
    /// \brief Создание сигма-точек пространств X и Y
    ///
    virtual void drawSigmaPoints()
    {
        // Создание сигма-точек пространства X
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
        // Вычисление сигма-точек пространства Y
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }

        // Вычисление X_pred по сигма-точкам пространства Х
        if( this->weightedSumStateSigmas_ == nullptr ) {
            this->X_pred_ = this->x_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->X_pred_ = this->weightedSumStateSigmas_( this->weights_mean_, this->x_pred_sigma_points_ );
        }
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }

        // Вычисление Y_pred по сигма-точкам пространства Y
        if( this->weightedSumMeasurementSigmas_ == nullptr ) {
            this->Y_pred_ = this->y_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->Y_pred_ = this->weightedSumMeasurementSigmas_( this->weights_mean_, this->y_pred_sigma_points_ );
        }
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_pred_ = this->checkBordersMeasurement_( this->Y_pred_ );
        }
    }
};

}

#endif // KALMAN_FILTER_EXTENDED_CUBATURE_H
/// \}
