//------------------------------------------------------------------------------
///
/// \file       kalman_eckf.h
/// \brief      Класс квадратно-корневого гибридного 
///             расширенно-кубатурного фильтра Калмана, КК-РКФК 
///             (square root extended-cubature Kalman filter, SR-ECKF)
/// \date       30.03.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_CUBATURE_SQUARE_ROOT_H
#define KALMAN_FILTER_EXTENDED_CUBATURE_SQUARE_ROOT_H

// Project includes:
#include <kalman_debug.h>
#include <kalman_srekf.h>
#include <kalman_srckf.h>
#include <qr_decomposition.h>

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief Класс квадратно-корневого гибридного
///        расширенно-кубатурного фильтра Калмана, КК-РКФК
///        (square root extended-cubature Kalman filter, SR-ECKF)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRECKF : 
    public CKalmanSREKF<SizeX, SizeY>,
    public CKalmanSRCKF<SizeX, SizeY>
{
public:    
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanSRECKF() = default;
    CKalmanSRECKF( const CKalmanSRECKF& ) = default;
    CKalmanSRECKF& operator=( const CKalmanSRECKF& ) = default;
    CKalmanSRECKF( CKalmanSRECKF&& ) = default;
    CKalmanSRECKF& operator=( CKalmanSRECKF&& ) = default;
    virtual ~CKalmanSRECKF() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "SRECKF"; }
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
        this->predictionSREKF( dt );
    }

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd ) override
    {
        this->drawSigmaPoints();
        this->correctionSRCKF( Y_msd );
    }

protected:
    //--------------------------------------------------------------------------
    ///
    /// \brief Создание сигма-точек пространств X и Y и переоценка X_pred, Y_pred
    ///
    void drawSigmaPoints()
    {
        // Создание сигма-точек пространства X
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
/// \brief Класс квадратно-корневого гибридного расширенно-кубатурного 
///        фильтра Калмана c блочной коррекцией, КК-РКФКБ 
///        (square root extended-cubature Kalman filter 
///        with block correction, SR-ECKFB)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSRECKFB : public CKalmanSRECKF<SizeX, SizeY>
{
public:
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanSRECKFB() = default;
    CKalmanSRECKFB( const CKalmanSRECKFB& ) = default;
    CKalmanSRECKFB& operator=( const CKalmanSRECKFB& ) = default;
    CKalmanSRECKFB( CKalmanSRECKFB&& ) = default;
    CKalmanSRECKFB& operator=( CKalmanSRECKFB&& ) = default;
    virtual ~CKalmanSRECKFB() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "SRECKFB"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Cubature; }
    virtual TCovType GetCovType() const override { return TCovType::CT_SquareRootMatrices; }    

    //--------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    // Prediction как в НЕблочном фильтре

    virtual void Correction( const arma::vec &Y_msd ) override
    {
        this->drawSigmaPoints();
        this->correctionSRCKFB( Y_msd );
    }
};

}

#endif // KALMAN_FILTER_EXTENDED_CUBATURE_SQUARE_ROOT_H
/// \}
