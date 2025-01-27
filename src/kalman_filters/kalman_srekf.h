//------------------------------------------------------------------------------
///
/// \file       kalman_srekf.h
/// \brief      Класс квадратно-корневого расширенного фильтра Калмана, КК-РФК
///             (square root extended Kalman filter, SR-EKF)
/// \date       11.01.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_SQUARE_ROOT_H
#define KALMAN_FILTER_EXTENDED_SQUARE_ROOT_H

// Project includes:
#include <kalman_debug.h>
#include <kalman_ekf.h>
#include <qr_decomposition.h> 

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief Класс квадратно-корневого расширенного фильтра Калмана, КК-РФК 
///        (square root extended Kalman filter, EKF)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSREKF : virtual public CKalmanEKF<SizeX, SizeY>
{
public:
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanSREKF() = default;
    CKalmanSREKF( const CKalmanSREKF& ) = default;
    CKalmanSREKF& operator=( const CKalmanSREKF& ) = default;
    CKalmanSREKF( CKalmanSREKF&& ) = default;
    CKalmanSREKF& operator=( CKalmanSREKF&& ) = default;
    virtual ~CKalmanSREKF() = default;

    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "SREKF"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Extended; }
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
#ifdef DEBUG_KALMAN
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << GetFilterName() + " Correction started" << std::endl;
        Y_msd.print( GetFilterName() + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone );
        this->prediction_isDone = false;

        // Вычисление матрицы S        
        arma::mat B1 = arma::join_horiz( ( this->H_ * this->P_ ), this->R_ );
        SPML::LQ::MGS_1( this->S_, B1 );
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( GetFilterName() + " Correction, S:" );
#endif        
        assert( this->indexOfNegativeDiagonalElement( this->S_ ) < 0 );

        // Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_ * arma::trans( this->P_ ) * arma::trans( this->H_ ) * arma::inv( this->S_ * arma::trans( this->S_ ) );
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
        arma::mat IKH = this->I_ - ( this->K_ * this->H_ );
        arma::mat B2 = arma::join_horiz( ( IKH * this->P_ ), ( this->K_ * this->R_ ) );       
        SPML::LQ::MGS_1( this->P_, B2 );
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

protected:
    ///
    /// \brief Прогноз SREKF
    /// \param dt - Время прогноза, [с]
    ///
    void predictionSREKF( double dt )
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
        ( this->P_ * arma::trans( this->P_ ) ).print( GetFilterName() + " Prediction, Pfull before:" );
#endif
        // 2. Вычисление матрицы P                
        arma::mat B = arma::join_horiz(
            ( this->F_ * this->P_ ), ( this->Q_ * std::sqrt( std::abs( dt ) ) ) 
        );
        SPML::LQ::MGS_1( this->P_, B );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( GetFilterName() + " Prediction, P after:" );
        ( this->P_ * arma::trans( this->P_ ) ).print( GetFilterName() + " Prediction, Pfull after:" );
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
};

//------------------------------------------------------------------------------
///
/// \brief Класс квадратно-корневого расширенного фильтра Калмана 
///        c блочной коррекцией, КК-РФКБ 
///        (square root extended Kalman filter with block correction, SREKFB)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSREKFB : public CKalmanSREKF<SizeX, SizeY>
{
public:    
    //--------------------------------------------------------------------------
    // Конструкторы:
    CKalmanSREKFB() = default;
    CKalmanSREKFB( const CKalmanSREKFB& ) = default;
    CKalmanSREKFB& operator=( const CKalmanSREKFB& ) = default;
    CKalmanSREKFB( CKalmanSREKFB&& ) = default;
    CKalmanSREKFB& operator=( CKalmanSREKFB&& ) = default;
    virtual ~CKalmanSREKFB() = default;
    
    //--------------------------------------------------------------------------
    virtual std::string GetFilterName() const override { return "SREKFB"; }
    virtual TKalmanType GetKalmanType() const override { return TKalmanType::KT_Extended; }
    virtual TCovType GetCovType() const override { return TCovType::CT_SquareRootMatrices; } 

    //--------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    // Прогноз как в НЕблочном фильтре

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

        // Составление блочной матрицы B:
        // [ R, H*P ]
        // [ 0, P   ]        
        arma::mat B = arma::join_vert(        
            arma::join_horiz( this->R_ , ( this->H_ * this->P_ ) ),
            arma::join_horiz( arma::mat( SizeX, SizeY, arma::fill::zeros ), this->P_ )
        );
        
        // LQ-разложение составленной блочной матрицы B
        arma::mat B_output;
        SPML::LQ::MGS_1( B_output, B );
		        
        // Считывание результата из матрицы B_output: 
        //  X.submat( first_row, first_col, last_row, last_col )
        // [ S,  0 ]
        // [ K*, P ]
        const size_t qrSize = SizeY + SizeX;
        this->S_ = B_output.submat( 0, 0, ( SizeY - 1 ), ( SizeY - 1 ) );
        this->K_ = ( B_output.submat( SizeY, 0, ( qrSize - 1 ), ( SizeY - 1 ) ) ) *         
            arma::inv( this->S_ );
        this->P_ = B_output.submat( SizeY, SizeY, ( qrSize - 1 ), ( qrSize - 1 ) );        

#ifdef DEBUG_KALMAN
        ( this->S_ ).print( GetFilterName() + " Correction, S:" );
        ( this->K_ ).print( GetFilterName() + " Correction, K:" );
        ( this->P_ ).print( GetFilterName() + " Correction, P:" );
#endif        
        // assert( this->indexOfNegativeDiagonalElement( this->S_ ) < 0 );
        // assert( this->indexOfNegativeDiagonalElement( this->P_ ) < 0 );

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
};

}

#endif // KALMAN_FILTER_EXTENDED_SQUARE_ROOT_H
/// \}
