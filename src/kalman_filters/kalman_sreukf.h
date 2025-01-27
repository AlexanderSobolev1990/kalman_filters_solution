//------------------------------------------------------------------------------
///
/// \file       kalman_eukf.h
/// \brief      Класс квадратно-корневого гибридного 
///             расширенно-сигма-точечного фильтра Калмана, КК-РСТФК 
///             (square root extended-unscented Kalman filter, SR-EUKF)
/// \date       11.04.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_UNSCENTED_SQUARE_ROOT_H
#define KALMAN_FILTER_EXTENDED_UNSCENTED_SQUARE_ROOT_H

// Project includes:
#include <kalman_debug.h>
#include <kalman_srekf.h>
#include <kalman_srukf.h>
#include <qr_decomposition.h>

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief Класс квадратно-корневого гибридного
///        расширенно-сигма-точечного фильтра Калмана, КК-РСТФК
///        (square root extended-unscented Kalman filter, SR-EUKF)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSREUKF : 
    public CKalmanSREKF<SizeX, SizeY>,
    public CKalmanSRUKF<SizeX, SizeY>
{
public:    
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanSREUKF() = default;
    CKalmanSREUKF( const CKalmanSREUKF& ) = default;
    CKalmanSREUKF& operator=( const CKalmanSREUKF& ) = default;
    CKalmanSREUKF( CKalmanSREUKF&& ) = default;
    CKalmanSREUKF& operator=( CKalmanSREUKF&& ) = default;
    virtual ~CKalmanSREUKF() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "SREUKF"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Unscented; }
    virtual TCovType GetCovType() const override { return TCovType::CT_SquareRootMatrices; }    

    //--------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///    
    virtual void Prediction( double dt ) override
    {
        this->predictionSREKF( dt );
    }

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd ) override
    {
        this->drawSigmaPoints();
        this->correctionSRUKF( Y_msd );
    }

protected:
    //--------------------------------------------------------------------------
    ///
    /// \brief Создание сигма-точек пространств X и Y и переоценка X_pred, Y_pred
    ///
    void drawSigmaPoints()
    {
        // Создание сигма-точек пространства X
        this->x_pred_sigma_points_.col(this->k_sigma_points_ - 1) = this->X_pred_; // Нулевая точка - сзади!
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->P_.col(i) );
            this->x_pred_sigma_points_.col(i) = this->X_pred_ + add;
            this->x_pred_sigma_points_.col(i + SizeX) = this->X_pred_ - add;
            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
                this->x_pred_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i) );
                this->x_pred_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i + SizeX) );
            }
        }

        // Вычисление сигма-точек пространства Y по сигма-точкам пространства X
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

//------------------------------------------------------------------------------
///
/// \brief Класс квадратно-корневого гибридного расширенно-сигма-точечного 
///        фильтра Калмана c блочной коррекцией, КК-РСТФК 
///        (square root extended-unscented Kalman filter 
///        with block correction, SR-EUKFB)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSREUKFB : public CKalmanSREUKF<SizeX, SizeY>
{
public:
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanSREUKFB() = default;
    CKalmanSREUKFB( const CKalmanSREUKFB& ) = default;
    CKalmanSREUKFB& operator=( const CKalmanSREUKFB& ) = default;
    CKalmanSREUKFB( CKalmanSREUKFB&& ) = default;
    CKalmanSREUKFB& operator=( CKalmanSREUKFB&& ) = default;
    virtual ~CKalmanSREUKFB() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "SREUKFB"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Unscented; }
    virtual TCovType GetCovType() const override { return TCovType::CT_SquareRootMatrices; }    

    //--------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    // Prediction как в НЕблочном фильтре

    virtual void Correction( const arma::vec &Y_msd ) override
    {
        this->drawSigmaPoints();
        this->correctionSRUKFB( Y_msd );
    }
};

}

#endif // KALMAN_FILTER_EXTENDED_UNSCENTED_SQUARE_ROOT_H
/// \}
