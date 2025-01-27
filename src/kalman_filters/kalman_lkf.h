//------------------------------------------------------------------------------
///
/// \file       kalman_lkf.h
/// \brief      Шаблонный класс линейного фильтра Калмана, ЛФК (Linear Kalman Filter, LKF)
/// \date       11.01.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_LINEAR_H
#define KALMAN_FILTER_LINEAR_H

// System includes:
#include <armadillo> // Матрицы
#include <cassert> // assert
#include <functional> // std::function

// Project includes
#include <kalman_debug.h>
#include <kalman_base.h>

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief Шаблонный класс линейного фильтра Калмана, ЛФК (linear Kalman filter, LKF)
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanLKF : public CKalmanBase<SizeX, SizeY>
{
public:    
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanLKF() = default;    
    CKalmanLKF( const CKalmanLKF& ) = default;
    CKalmanLKF& operator=( const CKalmanLKF& ) = default;
    CKalmanLKF( CKalmanLKF&& ) = default;
    CKalmanLKF& operator=( CKalmanLKF&& ) = default;
    virtual ~CKalmanLKF() = default;

    //--------------------------------------------------------------------------
    // Методы-сеттеры:

    ///
    /// \brief Установка матрицы перехода состояния F
    /// \param F - матрица перехода состояния, размерность [SizeX * SizeX]
    ///
    void SetStateTransitionMatrixF( const arma::mat &F )
    {
        if( arma::size( F ) != arma::size( this->F_ ) ) {
            throw std::length_error( "Incorrect dimensions of state transition matrix F" );
        } else {
            this->F_ = F;
        }
    }

    ///
    /// \brief Установка матрицы перехода измерений H
    /// \param H - матрица перехода измерений, размерность [SizeY * SizeX]
    ///
    void SetObservationMatrixH( const arma::mat &H )
    {
        if( arma::size( H ) != arma::size( this->H_ ) ) {
            throw std::length_error( "Incorrect dimensions of оbservation matrix H" );
        } else {
            this->H_ = H;
        }
    }

    ///
    /// \brief Установка функции вычисления матрицы перехода состояния F в случае LKF (makeMatrixF)
    /// \param stateTransitionJacobianLinearF - Функция вычисления матрицы перехода состояния F в случае LKF
    ///
    void SetStateTransitionJacobianLinearF( std::function<arma::mat( double dt )> stateTransitionJacobianLinearF )
    {
        this->stateTransitionJacobianLinearF_ = stateTransitionJacobianLinearF;
    }

    //--------------------------------------------------------------------------
    // Виртуальные методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << GetFilterName() + " Prediction started, dt = " << dt << std::endl;
#endif
        // 1. Вычисление X_pred, F
        this->F_ = this->stateTransitionJacobianLinearF_( dt ); // В общем случае матрица F зависит от времени, поэтому так
        this->X_pred_ = this->F_ * this->X_est_;
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }
#ifdef DEBUG_KALMAN
        this->F_.print( GetFilterName() + " Prediction, F:" );
        this->P_.print( GetFilterName() + " Prediction, P before:" );
#endif
        // 2. Вычисление ковариационной матрицы Р
        this->P_ = ( this->F_ * this->P_ * arma::trans( this->F_ ) ) + ( this->Q_ * std::abs( dt ) ); // Сразу прибавить Q (в случае, если матрица Q - плотная)
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        this->P_.print( GetFilterName() + " Prediction, P after:" );
#endif        
        assert( this->indexOfNegativeDiagonalElement( this->P_ ) < 0 );

        // 3. Вычисление X_est, Y_est
        this->Y_pred_ = this->H_ * this->X_pred_;
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_pred_ = this->checkBordersMeasurement_( this->Y_pred_ );
        }
        this->X_est_ = this->X_pred_;
        this->Y_est_ = this->Y_pred_;
#ifdef DEBUG_KALMAN
        this->X_est_.print( GetFilterName() + " Prediction, X_est:" );
        this->Y_est_.print( GetFilterName() + " Prediction, Y_est:" );
#endif
        this->prediction_isDone = true; // Выставить признак состоявшегося прогноза
    }

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << GetFilterName() + " Correction started" << std::endl;
        Y_msd.print( GetFilterName() + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла:
        this->prediction_isDone = false; // Сразу же снять признак

        // 1. Вычисление ковариационной матрицы S
        arma::mat Ht = arma::trans( this->H_ );
        this->S_ = ( this->H_ * this->P_ * Ht ) + this->R_; // Сразу прибавить R (в случае, если матрица R - плотная)
        this->fixMatrixMainDiagonalSymmetry( this->S_ );
#ifdef DEBUG_KALMAN
        this->S_.print( GetFilterName() + " Correction, S:" );
#endif        
        assert( this->indexOfNegativeDiagonalElement( this->S_ ) < 0 );

        // 2. Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_ * Ht * arma::inv( this->S_ );
#ifdef DEBUG_KALMAN
        this->K_.print( GetFilterName() + " Correction, K:" );
#endif
        // 3. Вычисление невязки Delta
        if( !this->deltaY_isSet ) {
            assert( this->Y_msd_isSet ); // Если не установлен deltaY, то Y_msd_ обязан быть установлен
            this->Y_msd_isSet = false; // Сразу же снять признак
            this->DeltaY_ = this->Y_msd_ - this->Y_pred_;
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->DeltaY_ = this->checkDeltaMeasurement_( this->DeltaY_ );
            }
        }
#ifdef DEBUG_KALMAN
        this->DeltaY_.print( GetFilterName() + " Correction, DeltaY:" );
#endif
        // 4. Вычисление ковариационной матрицы P
        this->P_ = ( this->I_ - ( this->K_ * this->H_ ) ) * this->P_;
//        this->P_ = ( ( ( this->K_ * this->H_ ) * ( -1.0 ) ).diag() + arma::vec( SizeX, arma::fill::ones ) ) * this->P_;

        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        this->P_.print( GetFilterName() + " Correction, P after:" );
#endif        
        assert( this->indexOfNegativeDiagonalElement( this->P_ ) < 0 );

        // 5. Вычисление X_est, Y_est
        this->X_est_ = this->X_pred_ + ( this->K_ * this->DeltaY_ );
        if( this->checkBordersStateAfterCorrection_ != nullptr ) {
            this->X_est_ = this->checkBordersStateAfterCorrection_( this->X_est_ ); // Проверка вектора состояния X на выход за допустимые пределы
        }
        this->Y_est_ = this->H_ * this->X_est_;
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_est_ = this->checkBordersMeasurement_( this->Y_est_ );
        }
#ifdef DEBUG_KALMAN
        this->X_est_.print( GetFilterName() + " Correction, X_est:" );
        this->Y_est_.print( GetFilterName() + " Correction, Y_est:" );
#endif
        this->deltaY_isSet = false; // Снять признак (выставляется в true в сеттере deltaY)
    }

    ///
    /// \brief Отдельное вычисление ковариационной матрицы S невязки измерений
    /// \details Отдельное вычисление S применяется при стробировании
    /// \attention Необходимо учесть 2 обстоятельства:
    /// 1) метод должен выполняться после вызова метода Prediction, где должна быть вычислена матрица H;
    /// 2) диагональ матрицы R при вычислении S будет использована та, что передана через Rdiag
    /// \param PdiagAdd - добавка в диагональ матрицы P, размерность [SizeX]
    /// \param Rdiag - диагональная матрица R априорных шумов измерений, размерность [SizeY]
    ///
    virtual void CalculateInnovationCovarianceS( const arma::vec &PdiagAdd, const arma::vec Rdiag )
    {        
        assert( this->prediction_isDone ); // Перед вы обязательно должен быть выполнен прогноз, иначе не имеет смысла:
        arma::mat Ptmp = this->P_;
        if( PdiagAdd.size() != 0 ) { // Если размер не нулевой - прибавим,
            assert( PdiagAdd.size() == SizeX ); // проверяя при этом, что размеры совпадают!
            Ptmp.diag() += PdiagAdd;
        }        
        this->S_ = ( this->H_ * Ptmp * arma::trans( this->H_ ) );
        this->S_.diag() += Rdiag;
        this->fixMatrixMainDiagonalSymmetry( this->S_ );
        assert( this->indexOfNegativeDiagonalElement( this->S_ ) < 0 );
    }

private:
    //--------------------------------------------------------------------------
    ///
    /// \brief Функция вычисления матрицы перехода состояния F (makeMatrixF)
    /// \attention Используется только в линейном фильтре Калмана
    /// \param X - вектор состояния с прошлого момента времени
    ///
    std::function<arma::mat( double dt )> stateTransitionJacobianLinearF_;
};

}

#endif // KALMAN_FILTER_LINEAR_H
/// \}
