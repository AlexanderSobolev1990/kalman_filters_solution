//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_extended.h
/// \brief      Шаблонный класс расширенного фильтра Калмана, РФК (Extended Kalman Filter, EKF)
/// \date       11.01.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_H
#define KALMAN_FILTER_EXTENDED_H

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_linear.h>

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Класс расширенного фильтра Калмана, РФК (Extended Kalman Filter, EKF)
/// \details Источник: NASA Technical report R-135, Application of statistical filter theory
/// to the optimal estimation of position and velocity on board a circumlunar vehicle,
/// Gerald L. Smith, Stanley F. Schmidt and Leonard A. McGee, 1962
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanEKF : public CKalmanLKF<SizeX, SizeY>
{
public:
    using CKalmanLKF<SizeX, SizeY>::CKalmanLKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanEKF() : CKalmanLKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "EKF" );
#endif
    }

    ///
    /// \brief Конструктор копирования
    /// \param other - экземпляр, с которого делается копия
    ///
    CKalmanEKF( const CKalmanEKF &other ) : CKalmanLKF<SizeX, SizeY>( other )
    {
        this->stateTransitionModel_ = other.stateTransitionModel_;
        this->observationModel_ = other.observationModel_;
        this->stateTransitionJacobianF_ = other.stateTransitionJacobianF_;
        this->observationJacobianH_ = other.observationJacobianH_;
    }

    ///
    /// \brief Перегрузка оператора присвоения
    /// \param other - экземпляр, с которого делается копия
    /// \return *this
    ///
    CKalmanEKF& operator=( const CKalmanEKF &other )
    {
        CKalmanEKF copy( other );
        swap( *this, copy );
        return *this;
    }

    ///
    /// \brief Конструктор перемещения
    /// \param other - экземпляр, с которого делается копия
    ///
    CKalmanEKF( CKalmanEKF &&other ) noexcept
    {
        swap( *this, other );
    }

    ///
    /// \brief Перегрузка оператора перемещения
    /// \param other - экземпляр, с которого делается копия
    /// \return *this
    ///
    CKalmanEKF& operator=( CKalmanEKF &&other ) noexcept
    {
        swap( *this, other );
        return *this;
    }

    virtual ~CKalmanEKF() = default; ///< Дестркутор

    ///
    /// \brief Метод свапа
    /// \param lhs - left hand side instance
    /// \param rhs - right hand side instance
    ///
    friend void swap( CKalmanEKF<SizeX, SizeY> &lhs, CKalmanEKF<SizeX, SizeY> &rhs ) noexcept
    {
        swap( dynamic_cast< CKalmanLKF<SizeX, SizeY> &>( lhs ), dynamic_cast< CKalmanLKF<SizeX, SizeY> &>( rhs ) ); // Свап предка 1 порядка

        std::swap( lhs.stateTransitionModel_, rhs.stateTransitionModel_ );
        std::swap( lhs.observationModel_, rhs.observationModel_ );
        std::swap( lhs.stateTransitionJacobianF_, rhs.stateTransitionJacobianF_ );
        std::swap( lhs.observationJacobianH_, rhs.observationJacobianH_ );
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы-сеттеры:

    ///
    /// \brief Установка функции прогноза состояния (predictState)
    /// \sa stateTransitionModel_
    ///
    void SetStateTransitionModel( std::function<arma::vec( const arma::vec &X, double dt )> stateTransitionModel )
    {
        stateTransitionModel_ = stateTransitionModel;
    }

    ///
    /// \brief Установка функции перевода состояния в измерение (XtoY)
    /// \sa observationModel_
    ///
    void SetObservationModel( std::function<arma::vec( const arma::vec &X )> observationModel )
    {
        observationModel_ = observationModel;
    }

    ///
    /// \brief Установка функции вычисления матрицы перехода состояния F (makeMatrixF)
    /// \sa stateTransitionJacobianF_
    ///
    void SetStateTransitionJacobianF( std::function<arma::mat( const arma::vec &X, double dt )> stateTransitionJacobianF )
    {
        stateTransitionJacobianF_ = stateTransitionJacobianF;
    }

    ///
    /// \brief Установка функции вычисления матрицы перехода измерений H (makeMatrixH)
    /// \sa observationJacobianH_
    ///
    void SetObservationJacobianH( std::function<arma::mat( const arma::vec &X )> observationJacobianH )
    {
        observationJacobianH_ = observationJacobianH;
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
        PredictionEKF( dt );
    }

    ///
    /// \brief Коррекция
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

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла:
        this->prediction_isDone = false; // Сразу же снять признак

        // 1. Вычисление ковариационной матрицы S
        arma::mat Ht = arma::trans( this->H_ );
        this->S_ = ( this->H_ * this->P_ * Ht ) + this->R_; // Сразу прибавить R (в случае, если матрица R - плотная)
        this->fixMatrixMainDiagonalSymmetry( this->S_ );        
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( this->filterName_ + " Correction, S:" );
#endif
        this->checkMatrixDiagPositive( this->S_ );

        // 2. Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_ * Ht * arma::inv( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->K_ ).print( this->filterName_ + " Correction, K:" );
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
        ( this->DeltaY_ ).print( this->filterName_ + " Correction, DeltaY:" );
#endif
        // 4. Вычисление ковариационной матрицы P
        this->P_ = ( this->I_ - ( this->K_ * this->H_ ) ) * this->P_;
        this->fixMatrixMainDiagonalSymmetry( this->P_ );        
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Correction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );

        // 5. Вычисление X_est, Y_est
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
    ///
    /// \brief Прогноз EKF
    /// \param dt - Время прогноза, [с]
    ///
    void PredictionEKF( double dt )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Prediction started, dt = " << dt << std::endl;
#endif
        assert( dt > 0.0 ); // Прогноз имеет смысл только при положительном времени

        // 1. Вычисление X_pred, F
        this->X_pred_ = this->stateTransitionModel_( this->X_est_, dt );
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }
        this->F_ = this->stateTransitionJacobianF_( this->X_est_, dt );
#ifdef DEBUG_KALMAN
        ( this->F_ ).print( this->filterName_ + " Prediction, F:" );
        ( this->P_ ).print( this->filterName_ + " Prediction, P before:" );
#endif
        // 2. Вычисление ковариационной матрицы P
        this->P_ = this->F_ * this->P_ * arma::trans( this->F_ ) + ( this->Q_ * dt ); // Сразу прибавить Q (в случае, если матрица Q - плотная)
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Prediction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );

        // 3. Вычисление матрицы перехода H
        this->H_ = this->observationJacobianH_( this->X_pred_ );

        // 4. Вычисление X_est, Y_est
        this->Y_pred_ = this->observationModel_( this->X_pred_ );
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_pred_ = this->checkBordersMeasurement_( this->Y_pred_ );
        }
        this->X_est_ = this->X_pred_;
        this->Y_est_ = this->Y_pred_;
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( this->filterName_ + " Prediction, X_est:" );
        ( this->Y_est_ ).print( this->filterName_ + " Prediction, Y_est:" );
#endif
        this->prediction_isDone = true; // Выставить признак состоявшегося прогноза
    }

    //------------------------------------------------------------------------------------------------------------------
    // Обертки функций EKF
    ///
    /// \brief Функция прогноза состояния (predictState)
    /// \param X - вектор состояния с прошлого момента времени
    /// \sa SetStateTransitionModel
    ///
    std::function<arma::vec( const arma::vec &X, double dt )> stateTransitionModel_;

    ///
    /// \brief Функция перевода состояния в измерения (XtoY)
    /// \param X - вектор состояния текущего момента времени
    /// \sa SetObservationModel
    ///
    std::function<arma::vec( const arma::vec &X )> observationModel_;

    ///
    /// \brief Функция вычисления матрицы перехода состояния F (makeMatrixF)
    /// \param X - вектор состояния с прошлого момента времени
    /// \sa SetStateTransitionJacobianF
    ///
    std::function<arma::mat( const arma::vec &X, double dt )> stateTransitionJacobianF_;

    ///
    /// \brief Функция вычисления матрицы перехода измерений H (makeMatrixH)
    /// \param X - вектор состояния текущего момента времени
    /// \sa SetObservationJacobianH
    ///
    std::function<arma::mat( const arma::vec &X )> observationJacobianH_;

private:
    //------------------------------------------------------------------------------------------------------------------
    // Запрет доступа к методам родительского класса (поскольку они допустимы только для линейного фильтра LKF):
    void SetStateTransitionJacobianLinearF(){}; ///< Запрет доступа
    void SetStateTransitionMatrixF(){}; ///< Запрет доступа
    void SetObservationMatrixH(){}; ///< Запрет доступа
};

/*
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Класс расширенного фильтра Калмана Бирмана-Тронтона, РФК-БТ (Extended Kalman Filter Bierman-Thornton, EKF-UD)
/// \details Источник: Comparison of Nonlinear Filtering Techniques for Lunar Surface Roving Navigation
/// Kimber Lemon and Bryan W. Welch Glenn Research Center, Cleveland, Ohio, NASA/TM—2008-215152
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanEKFUD : public CKalmanEKF<SizeX, SizeY>
{
public:
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanEKFUD() : CKalmanEKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "EKFUD" );
#endif
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
        PredictionEKFUD( dt );
    }

    ///
    /// \brief Коррекция
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

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла:
        this->prediction_isDone = false; // Сразу же снять признак

        // 1. Вычисление ковариационной матрицы S
        arma::mat Ht = arma::trans( this->H_ );
        this->S_ = ( this->H_ * this->P_ * Ht ) + this->R_; // Сразу прибавить R (в случае, если матрица R - плотная)
        this->fixMatrixMainDiagonalSymmetry( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( this->filterName_ + " Correction, S:" );
#endif
        this->checkMatrixDiagPositive( this->S_ );

        // 2. Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_ * Ht * arma::inv( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->K_ ).print( this->filterName_ + " Correction, K:" );
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
        ( this->DeltaY_ ).print( this->filterName_ + " Correction, DeltaY:" );
#endif
        // 4. Вычисление ковариационной матрицы P
        this->P_ = ( this->I_ - ( this->K_ * this->H_ ) ) * this->P_;
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Correction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );

        // 5. Вычисление X_est, Y_est
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
    ///
    /// \brief Прогноз EKF UD
    /// \param dt - Время прогноза, [с]
    ///
    void PredictionEKFUD( double dt )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Prediction started, dt = " << dt << std::endl;
#endif
        assert( dt > 0.0 ); // Прогноз имеет смысл только при положительном времени

        // 1. Вычисление X_pred, F
        this->X_pred_ = this->stateTransitionModel_( this->X_est_, dt );
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }
        this->F_ = this->stateTransitionJacobianF_( this->X_est_, dt );
#ifdef DEBUG_KALMAN
        ( this->F_ ).print( this->filterName_ + " Prediction, F:" );
        ( this->P_ ).print( this->filterName_ + " Prediction, P before:" );
#endif
        // 2. Вычисление ковариационной матрицы P
        arma::mat W = arma::join_horiz( this->F_ * this->U_, this->I_ );
//        arma::mat WT = arma::join_horiz( arma::trans( this->U_ ) * arma::trans( this->F_), this->I_ );
        arma::mat D = arma::mat( 2 * SizeX, 2 * SizeX );
        arma::mat Qdt = this->Q_ * std::sqrt( dt );
        for( int i = 0; i < ( 2 * SizeX ); i++ ) {
            for( int j = 0; j < ( 2 * SizeX ); j++ ) {
                if( ( i < SizeX ) && ( j < SizeX ) ) {
                    D(i, j) = this->D_;
                } else if( ( i >= SizeX ) && ( j >= SizeX ) ) {
                    D(i, j) = Qdt(i, j);
                }
            }
        }
        this->P_ = W * D * arma::trans( W );// WT;
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Prediction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );

        // 3. Вычисление матрицы перехода H
        this->H_ = this->observationJacobianH_( this->X_pred_ );

        // 4. Вычисление X_est, Y_est
        this->Y_pred_ = this->observationModel_( this->X_pred_ );
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_pred_ = this->checkBordersMeasurement_( this->Y_pred_ );
        }
        this->X_est_ = this->X_pred_;
        this->Y_est_ = this->Y_pred_;
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( this->filterName_ + " Prediction, X_est:" );
        ( this->Y_est_ ).print( this->filterName_ + " Prediction, Y_est:" );
#endif
        this->prediction_isDone = true; // Выставить признак состоявшегося прогноза
    }

private:
    //------------------------------------------------------------------------------------------------------------------
    arma::mat U_; // Матрица UDU' разложения матрицы P
    arma::mat D_; // Матрица UDU' разложения матрицы P
};
*/

}

#endif // KALMAN_FILTER_EXTENDED_H
/// \}
