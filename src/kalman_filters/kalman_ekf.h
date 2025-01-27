//------------------------------------------------------------------------------
///
/// \file       kalman_ekf.h
/// \brief      Класс расширенного фильтра Калмана, РФК
///             (extended Kalman filter, EKF)
/// \date       11.01.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_H
#define KALMAN_FILTER_EXTENDED_H

// Project includes:
#include <kalman_debug.h>
#include <kalman_base.h>

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief Класс расширенного фильтра Калмана, РФК 
///        (extended Kalman filter, EKF)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanEKF : public CKalmanBase<SizeX, SizeY>
{
public:
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanEKF() = default;
    CKalmanEKF( const CKalmanEKF& ) = default;
    CKalmanEKF& operator=( const CKalmanEKF& ) = default;
    CKalmanEKF( CKalmanEKF&& ) = default;
    CKalmanEKF& operator=( CKalmanEKF&& ) = default;
    virtual ~CKalmanEKF() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "EKF"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Extended; }
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
#ifdef DEBUG_KALMAN
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << GetFilterName() + " Correction started" << std::endl;
        Y_msd.print( GetFilterName() + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла:
        this->prediction_isDone = false; // Сразу же снять признак

        // Вычисление матрицы S
        arma::mat Ht = arma::trans( this->H_ );
        this->S_ = ( this->H_ * this->P_ * Ht ) + this->R_; // Сразу прибавить R (в случае, если матрица R - плотная)
        this->fixMatrixMainDiagonalSymmetry( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( GetFilterName() + " Correction, S:" );
#endif        
        assert( this->indexOfNegativeDiagonalElement( this->S_ ) < 0 );

        // Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_ * Ht * arma::inv( this->S_ );
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
        this->P_ = ( this->I_ - ( this->K_ * this->H_ ) ) * this->P_;
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

    // ///
    // /// \brief Отдельное вычисление ковариационной матрицы S невязки измерений
    // /// \details Отдельное вычисление S применяется при стробировании
    // /// \attention Необходимо учесть 2 обстоятельства:
    // /// 1) метод должен выполняться после вызова метода Prediction, где должна быть вычислена матрица H;
    // /// 2) диагональ матрицы R при вычислении S будет использована та, что передана через Rdiag
    // /// \param PdiagAdd - добавка в диагональ матрицы P, размерность [SizeX]
    // /// \param Rdiag - диагональная матрица R априорных шумов измерений, размерность [SizeY]
    // ///
    // virtual void CalculateInnovationCovarianceS( 
    //     const arma::vec &PdiagAdd, const arma::vec Rdiag ) override
    // {
    //     assert( this->prediction_isDone ); // Перед вы обязательно должен быть выполнен прогноз, иначе не имеет смысла:
    //     arma::mat Ptmp = this->P_;
    //     if( PdiagAdd.size() != 0 ) { // Если размер не нулевой - прибавим,
    //         assert( PdiagAdd.size() == SizeX ); // проверяя при этом, что размеры совпадают!
    //         Ptmp.diag() += PdiagAdd;
    //     }        
    //     this->S_ = ( this->H_ * Ptmp * arma::trans( this->H_ ) );
    //     this->S_.diag() += Rdiag;
    //     this->fixMatrixMainDiagonalSymmetry( this->S_ );
    //     assert( this->indexOfNegativeDiagonalElement( this->S_ ) < 0 );
    // }

protected:
    //--------------------------------------------------------------------------
    ///
    /// \brief Прогноз EKF
    /// \param dt - Время прогноза, [с]
    ///
    void predictionEKF( double dt )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << GetFilterName() + " Prediction started, dt = " << dt << std::endl;
#endif
        // Вычисление X_pred, F
        this->X_pred_ = this->stateTransitionModel_( this->X_est_, dt );
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }
        this->F_ = this->stateTransitionJacobianF_( this->X_est_, dt );
#ifdef DEBUG_KALMAN
        ( this->X_pred_ ).print( GetFilterName() + " Prediction, X_pred_:" );
        ( this->F_ ).print( GetFilterName() + " Prediction, F:" );
        ( this->P_ ).print( GetFilterName() + " Prediction, P before:" );
#endif
        // Вычисление матрицы P
        this->P_ = this->F_ * this->P_ * arma::trans( this->F_ ) + ( this->Q_ * std::abs( dt ) ); // Сразу прибавить Q (в случае, если матрица Q - плотная)
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( GetFilterName() + " Prediction, P after:" );
#endif
        assert( this->indexOfNegativeDiagonalElement( this->P_ ) < 0 );

        // Вычисление матрицы H
        this->H_ = this->observationJacobianH_( this->X_pred_ );

        // Вычисление X_est, Y_est
        this->Y_pred_ = this->observationModel_( this->X_pred_ );
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_pred_ = this->checkBordersMeasurement_( this->Y_pred_ );
        }
        this->X_est_ = this->X_pred_;
        this->Y_est_ = this->Y_pred_;
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( GetFilterName() + " Prediction, X_est:" );
        ( this->Y_est_ ).print( GetFilterName() + " Prediction, Y_est:" );
#endif
        this->prediction_isDone = true; // Выставить признак состоявшегося прогноза
    }

    virtual void SetDesignParametersMeanSet( double w0 ){}
    virtual void SetDesignParametersScaledSet( double alpha, double beta, double kappa ){}
    virtual void SetDesignParametersCubatureBaseSet(){}
};

}

#endif // KALMAN_FILTER_EXTENDED_H
/// \}
