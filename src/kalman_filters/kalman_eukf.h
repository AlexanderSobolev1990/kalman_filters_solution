//------------------------------------------------------------------------------
///
/// \file       kalman_eukf.h
/// \brief      Шаблонный класс расширенного сигма-точечного (ансцентного) 
///             фильтра Калмана, РСТФК (extended unscented Kalman filter, EUKF)
/// \date       11.04.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_UNSCENTED_H
#define KALMAN_FILTER_EXTENDED_UNSCENTED_H

// Project includes:
#include <kalman_debug.h>
#include <kalman_ukf.h>

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief  Шаблонный класс гибридного расширенно-сигма-точечного (ансцентного)
///         фильтра Калмана, РСТФК (extended unscented Kalman filter, EUKF)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanEUKF : public CKalmanUKF<SizeX, SizeY>
{
public:    
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanEUKF() = default;
    CKalmanEUKF( const CKalmanEUKF& ) = default;
    CKalmanEUKF& operator=( const CKalmanEUKF& ) = default;
    CKalmanEUKF( CKalmanEUKF&& ) = default;
    CKalmanEUKF& operator=( CKalmanEUKF&& ) = default;
    virtual ~CKalmanEUKF() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "EUKF"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Unscented; }
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
        this->correctionUKF( Y_msd );
    }

protected:
    //--------------------------------------------------------------------------
    ///
    /// \brief Создание сигма-точек пространств X и Y
    ///
    void drawSigmaPoints()
    {
        // Создание сигма-точек пространства X
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

#endif // KALMAN_FILTER_EXTENDED_UNSCENTED_H
/// \}
