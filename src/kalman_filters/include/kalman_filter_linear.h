//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_linear.h
/// \brief      Шаблонный класс линейного фильтра Калмана, ЛФК (Linear Kalman Filter, LKF)
/// \date       11.01.21 - создан
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
#include <kalman_filter_debug.h>

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс линейного фильтра Калмана, ЛФК (Linear Kalman Filter, LKF)
/// \details Источник: A New Approach to Linear Filtering and Prediction Problems, R.E.KALMAN,
/// Research Institute for Advanced Study, Baltimore, Md., 1960
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanLKF
{
public:
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanLKF() :
#ifdef DEBUG_KALMAN
        filterName_( "LKF" ),
#endif
        SizeX_( SizeX ),
        SizeY_( SizeY ),

        F_( arma::eye( SizeX_, SizeX_ ) ),
        H_( arma::zeros( SizeY_, SizeX_ ) ),

        K_( arma::zeros( SizeX_, SizeY_ ) ),
        I_( arma::eye( SizeX_, SizeX_ ) ),
        DeltaY_( arma::zeros( SizeY_ ) ),

        P_( arma::zeros( SizeX_, SizeX_ ) ),
        S_( arma::zeros( SizeY_, SizeY_ ) ),
        Q_( arma::zeros( SizeX_, SizeX_ ) ),
        R_( arma::zeros( SizeY_, SizeY_ ) ),

        X_pred_( arma::zeros( SizeX_ ) ),
        Y_pred_( arma::zeros( SizeY_ ) ),
        X_est_( arma::zeros( SizeX_ ) ),
        Y_est_( arma::zeros( SizeY_ ) ),
        Y_msd_( arma::zeros( SizeY_ ) ),

        deltaY_isSet( false ),
        Y_msd_isSet( false ),
        prediction_isDone( false )
    {}

    ///
    /// \brief Конструктор копирования
    /// \param other - экземпляр, с которого делается копия
    ///
    CKalmanLKF( const CKalmanLKF &other )
    {
        this->SizeX_ = other.SizeX_;
        this->SizeY_ = other.SizeY_;

        this->F_ = other.F_;
        this->H_ = other.H_;

        this->K_ = other.K_;
        this->I_ = other.I_;
        this->DeltaY_ = other.DeltaY_;

        this->P_ = other.P_;
        this->S_ = other.S_;
        this->Q_ = other.Q_;
        this->R_ = other.R_;

        this->X_pred_ = other.X_pred_;
        this->Y_pred_ = other.Y_pred_;
        this->X_est_ = other.X_est_;
        this->Y_est_ = other.Y_est_;
        this->Y_msd_ = other.Y_msd_;

        this->deltaY_isSet = other.deltaY_isSet;
        this->Y_msd_isSet = other.Y_msd_isSet;
        this->prediction_isDone = other.prediction_isDone;
    }

    ///
    /// \brief Перегрузка оператора присвоения
    /// \param other - экземпляр, с которого делается копия
    /// \return *this
    ///
    CKalmanLKF& operator=( const CKalmanLKF &other )
    {
        CKalmanLKF copy( other );
        swap( *this, copy );
        return *this;
    }

    ///
    /// \brief Конструктор перемещения
    /// \param other - экземпляр, с которого делается копия
    ///
    CKalmanLKF( CKalmanLKF &&other ) noexcept :
        SizeX_( SizeX ),
        SizeY_( SizeY ),

        F_( arma::eye( SizeX_, SizeX_ ) ),
        H_( arma::zeros( SizeY_, SizeX_ ) ),

        K_( arma::zeros( SizeX_, SizeY_ ) ),
        I_( arma::eye( SizeX_, SizeX_ ) ),
        DeltaY_( arma::zeros( SizeY_ ) ),

        P_( arma::zeros( SizeX_, SizeX_ ) ),
        S_( arma::zeros( SizeY_, SizeY_ ) ),
        Q_( arma::zeros( SizeX_, SizeX_ ) ),
        R_( arma::zeros( SizeY_, SizeY_ ) ),

        X_pred_( arma::zeros( SizeX_ ) ),
        Y_pred_( arma::zeros( SizeY_ ) ),
        X_est_( arma::zeros( SizeX_ ) ),
        Y_est_( arma::zeros( SizeY_ ) ),
        Y_msd_( arma::zeros( SizeY_ ) ),

        deltaY_isSet( false ),
        Y_msd_isSet( false ),
        prediction_isDone( false )
    {
        swap( *this, other );
    }

    ///
    /// \brief Перегрузка оператора перемещения
    /// \param other - экземпляр, с которого делается копия
    /// \return *this
    ///
    CKalmanLKF& operator=( CKalmanLKF &&other ) noexcept
    {
        swap( *this, other );
        return *this;
    }

    virtual ~CKalmanLKF() = default; ///< Дестркутор

    ///
    /// \brief Метод свапа
    /// \param lhs - left hand side instance
    /// \param rhs - right hand side instance
    ///
    friend void swap( CKalmanLKF<SizeX, SizeY> &lhs, CKalmanLKF<SizeX, SizeY> &rhs ) noexcept
    {
        std::swap( lhs.SizeX_, rhs.SizeX_ );
        std::swap( lhs.SizeY_, rhs.SizeY_ );

        std::swap( lhs.F_, rhs.F_ );
        std::swap( lhs.H_, rhs.H_ );

        std::swap( lhs.K_, rhs.K_ );
        std::swap( lhs.I_, rhs.I_ );
        std::swap( lhs.DeltaY_, rhs.DeltaY_ );

        std::swap( lhs.P_, rhs.P_ );
        std::swap( lhs.S_, rhs.S_ );
        std::swap( lhs.Q_, rhs.Q_ );
        std::swap( lhs.R_, rhs.R_ );

        std::swap( lhs.X_pred_, rhs.X_pred_ );
        std::swap( lhs.Y_pred_, rhs.Y_pred_ );
        std::swap( lhs.X_est_, rhs.X_est_ );
        std::swap( lhs.Y_est_, rhs.Y_est_ );
        std::swap( lhs.Y_msd_, rhs.Y_msd_ );

        std::swap( lhs.deltaY_isSet, rhs.deltaY_isSet );
        std::swap( lhs.Y_msd_isSet, rhs.Y_msd_isSet );
        std::swap( lhs.prediction_isDone, rhs.prediction_isDone );
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы-сеттеры:
#ifdef DEBUG_KALMAN
    ///
    /// \brief Установка имени фильтра
    /// \param name
    ///
    void SetFilterName( std::string name )
    {
        this->filterName_ = name;
    }
#endif
    ///
    /// \brief Установка матрицы перехода состояния F
    /// \param F - матрица перехода состояния, размерность [SizeX * SizeX]
    ///
    void SetStateTransitionMatrixF( const arma::mat &F )
    {
        if( arma::size( F ) != arma::size( F_ ) ) {
            throw std::length_error( "Incorrect dimensions of state transition matrix F" );
        } else {
            F_ = F;
        }
    }

    ///
    /// \brief Установка матрицы перехода измерений H
    /// \param H - матрица перехода измерений, размерность [SizeY * SizeX]
    ///
    void SetObservationMatrixH( const arma::mat &H )
    {
        if( arma::size( H ) != arma::size( H_ ) ) {
            throw std::length_error( "Incorrect dimensions of оbservation matrix H" );
        } else {
            H_ = H;
        }
    }

    ///
    /// \brief Установка ковариационной матрицы P состояния X
    /// \param P - матрица шумов состояния X, размерность [SizeX * SizeX]
    ///
    void SetEstimateCovarianceMatrixP( const arma::mat &P )
    {
        if( arma::size( P ) != arma::size( P_ ) ) {
            throw std::length_error( "Incorrect dimensions of estimate covariance matrix P" );
        } else {
            P_ = P;
        }
    }

    ///
    /// \brief Установка диагонали ковариационной матрицы P состояния X
    /// \param Pdiag - вектор главной диагонали шумов состояния X (остальные элементы полагаются равными нулю), размерность [SizeX]
    ///
    void SetEstimateCovarianceMatrixPdiag( const arma::vec &Pdiag )
    {
        if( arma::size( Pdiag ) != arma::size( P_.diag() ) ) {
            throw std::length_error( "Incorrect dimensions of estimate covariance matrix P" );
        } else {
            P_.diag() = Pdiag;
        }
    }

    ///
    /// \brief Установка ковариационной матрицы Q шумов состояния X
    /// \param Q - матрица Q, размерность [SizeX, SizeX]
    ///
    void SetProcessCovarianceMatrixQ( const arma::mat &Q )
    {
        if( arma::size( Q ) != arma::size( Q_ ) ) {
            throw std::length_error( "Incorrect dimensions of process covariance matrix Q" );
        } else {
            Q_ = Q;
        }
    }

    ///
    /// \brief Установка диагонали ковариационной матрицы Q шумов состояния X
    /// \param Qdiag - вектор элементов главное диагонали матрицы Q (остальные элементы полагаются равными нулю), размерность [SizeX]
    ///
    void SetProcessCovarianceMatrixQdiag( const arma::vec &Qdiag )
    {
        if( arma::size( Qdiag ) != arma::size( Q_.diag() ) ) {
            throw std::length_error( "Incorrect dimensions of process covariance matrix Q" );
        } else {
            Q_.diag() = Qdiag;
        }
    }

    ///
    /// \brief Установка ковариационной матрицы R шумов измерений Y
    /// \param R - матрица R, размерность [SizeY, SizeY]
    ///
    void SetObservationCovarianceMatrixR( const arma::mat &R )
    {
        if( arma::size( R ) != arma::size( R_ ) ) {
            throw std::length_error( "Incorrect dimensions of оbservation covariance matrix R" );
        } else {
            R_ = R;
        }
    }

    ///
    /// \brief Установка ковариационной матрицы R шумов измерений Y
    /// \param Rdiag - вектор элементов главной диагонали матрицы R (остальные элементы полагаются равными нулю), размерность [SizeY]
    ///
    void SetObservationCovarianceMatrixRdiag( const arma::vec &Rdiag )
    {
        if( arma::size( Rdiag ) != arma::size( R_.diag() ) ) {
            throw std::length_error( "Incorrect dimensions of оbservation covariance matrix R" );
        } else {
            R_.diag() = Rdiag;
        }
    }

    ///
    /// \brief Установка оценки вектора состояния X
    /// \details Требуется, например, при начальной установке фильтра
    /// \param X_est - вектор состояния, размерность [SizeX]
    ///
    void SetEstimatedVectorX( const arma::vec &X_est )
    {
        if( arma::size( X_est ) != arma::size( X_est_ ) ) {
            throw std::length_error( "Incorrect dimensions of state vector X_est" );
        } else {
            X_est_ = X_est;
        }
    }

    ///
    /// \brief Установка оценки вектора состояния Y
    /// \details Требуется, например, при начальной установке фильтра
    /// \param Y_est - вектор состояния, размерность [SizeY]
    ///
    void SetEstimatedVectorY( const arma::vec &Y_est )
    {
        if( arma::size( Y_est ) != arma::size( Y_est_ ) ) {
            throw std::length_error( "Incorrect dimensions of state vector Y_est" );
        } else {
            Y_est_ = Y_est;
        }
    }

    ///
    /// \brief Установка измереннго вектора измерений Y
    /// \details Требуется после получения новых измерений
    /// \param Y_msd - вектор измерений, размерность [SizeY]
    ///
    void SetMeasuredVectorY( const arma::vec &Y_msd )
    {
        if( arma::size( Y_msd ) != arma::size( Y_msd_ ) ) {
            throw std::length_error( "Incorrect dimensions of measurement vector Y_msd" );
        } else {
            Y_msd_ = Y_msd;
            Y_msd_isSet = true;
        }
    }

    ///
    /// \brief Установка вектора невязки измерений DeltaY
    /// \param DeltaY - вектор невязки измерений, размерность [SizeY]
    ///
    void SetDeltaY( const arma::vec &DeltaY )
    {
        if( arma::size( DeltaY ) != arma::size( DeltaY_ ) ) {
            throw std::length_error( "Incorrect dimensions of vector DeltaY" );
        } else {
            DeltaY_ = DeltaY;
            deltaY_isSet = true;
        }
    }

    ///
    /// \brief Установка функции проверки вектора состояния X после прогноза
    /// \param checkBordersStateAfterPrediction - Функция проверки вектора состояния после прогноза
    void SetCheckBordersStateAfterPrediction( std::function<arma::vec( const arma::vec &X )> checkBordersStateAfterPrediction )
    {
        checkBordersStateAfterPrediction_ = checkBordersStateAfterPrediction;
    }

    ///
    /// \brief Установка функции проверки вектора состояния X после коррекции
    /// \param checkBordersStateAfterCorrection - Функция проверки вектора состояния после коррекции
    ///
    void SetCheckBordersStateAfterCorrection( std::function<arma::vec( const arma::vec &X )> checkBordersStateAfterCorrection )
    {
        checkBordersStateAfterCorrection_ = checkBordersStateAfterCorrection;
    }

    ///
    /// \brief Установка функции проверки вектора измерений Y
    /// \param checkBordersMeasurement - Функция проверки вектора измерений
    ///
    void SetCheckBordersMeasurement( std::function<arma::vec( const arma::vec &Y )> checkBordersMeasurement )
    {
        checkBordersMeasurement_ = checkBordersMeasurement;
    }

    ///
    /// \brief Установка функции проверки разности векторов состояний X
    /// \param checkDeltaState - Функция проверки разности векторов состояний
    ///
    void SetCheckDeltaState( std::function<arma::vec( const arma::vec &DeltaX )> checkDeltaState )
    {
        checkDeltaState_ = checkDeltaState;
    }

    ///
    /// \brief Установка функции проверки разности векторов измерений Y
    /// \param checkDeltaMeasurement - Функция проверки разности векторов измерений
    ///
    void SetCheckDeltaMeasurement( std::function<arma::vec( const arma::vec &DeltaY )> checkDeltaMeasurement )
    {
        checkDeltaMeasurement_ = checkDeltaMeasurement;
    }

    ///
    /// \brief Установка функции вычисления матрицы перехода состояния F в случае LKF (makeMatrixF)
    /// \param stateTransitionJacobianLinearF - Функция вычисления матрицы перехода состояния F в случае LKF
    ///
    void SetStateTransitionJacobianLinearF( std::function<arma::mat( double dt )> stateTransitionJacobianLinearF )
    {
        stateTransitionJacobianLinearF_ = stateTransitionJacobianLinearF;
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы-геттеры:
#ifdef DEBUG_KALMAN
    ///
    /// \brief Получить название фильтра
    /// \return Строка с именем фильтра
    ///
    const std::string GetFilterName() const { return filterName_; }
#endif
    ///
    /// \brief Получить размерность вектора состояния X
    /// \return Размерность вектора состояния X
    ///
    const int GetSizeX() const { return SizeX_; }

    ///
    /// \brief Получить размерность вектора состояния Y
    /// \return Размерность вектора состояния Y
    ///
    const int GetSizeY() const { return SizeY_; }

    ///
    /// \brief Получить ковариационную матрицу P состояния X
    /// \return Текущая ковариационная матрица P
    ///
    const arma::mat &GetEstimatedCovarianceMatrixP() const { return P_; }

    ///
    /// \brief Получить ковариационную матрицу S вектора невязки DeltaY
    /// \return Текущая ковариационная матрица S
    ///
    const arma::mat &GetInnovationCovarianceMatrixS() const { return S_; }

    ///
    /// \brief Получить вектор невязки измерений DeltaY
    /// \return Вектор невязки измерений DeltaY
    ///
    const arma::mat &GetDeltaY() const { return DeltaY_; }

    ///
    /// \brief Получить матрицу коэффициентов усиления фильтра K
    /// \return Текущая матрица коэффициентов усиления фильтра K
    ///
    const arma::mat &GetKalmanGainMatrixK() const { return K_; }

    ///
    /// \brief Получить уточненный вектор состояния X
    /// \return Текущий уточненный вектор состояния X
    ///
    const arma::mat &GetEstimatedVectorX() const { return X_est_; }

    ///
    /// \brief Получить уточненный вектор состояния Y
    /// \return Текущий уточненный вектор состояния Y
    ///
    const arma::mat &GetEstimatedVectorY() const { return Y_est_; }

    //------------------------------------------------------------------------------------------------------------------
    // Виртуальные методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Prediction started, dt = " << dt << std::endl;
#endif
        assert( dt > 0.0 ); // Прогноз имеет смысл только при положительном времени

        // 1. Вычисление X_pred, F
        F_ = stateTransitionJacobianLinearF_( dt ); // В общем случае матрица F зависит от времени, поэтому так
        X_pred_ = F_ * X_est_;
        if( checkBordersStateAfterPrediction_ != nullptr ) {
            X_pred_ = checkBordersStateAfterPrediction_( X_pred_ );
        }
#ifdef DEBUG_KALMAN
        F_.print( this->filterName_ + " Prediction, F:" );
        P_.print( this->filterName_ + " Prediction, P before:" );
#endif
        // 2. Вычисление ковариационной матрицы Р
        P_ = ( F_ * P_ * arma::trans( F_ ) ) + ( Q_ * dt ); // Сразу прибавить Q (в случае, если матрица Q - плотная)
        fixMatrixMainDiagonalSymmetry( P_ );        
#ifdef DEBUG_KALMAN
        P_.print( this->filterName_ + " Prediction, P after:" );
#endif
        checkMatrixDiagPositive( P_ );

        // 3. Вычисление X_est, Y_est
        Y_pred_ = H_ * X_pred_;
        if( checkBordersMeasurement_ != nullptr ) {
            Y_pred_ = checkBordersMeasurement_( Y_pred_ );
        }
        X_est_ = X_pred_;
        Y_est_ = Y_pred_;
#ifdef DEBUG_KALMAN
        X_est_.print( this->filterName_ + " Prediction, X_est:" );
        Y_est_.print( this->filterName_ + " Prediction, Y_est:" );
#endif
        prediction_isDone = true; // Выставить признак состоявшегося прогноза
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
        SetMeasuredVectorY( Y_msd );

        assert( prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла:
        prediction_isDone = false; // Сразу же снять признак

        // 1. Вычисление ковариационной матрицы S
        arma::mat Ht = arma::trans( H_ );
        S_ = ( H_ * P_ * Ht ) + R_; // Сразу прибавить R (в случае, если матрица R - плотная)
        fixMatrixMainDiagonalSymmetry( S_ );        
#ifdef DEBUG_KALMAN
        S_.print( this->filterName_ + " Correction, S:" );
#endif
        checkMatrixDiagPositive( S_ );

        // 2. Вычисление коэффициента усиления фильтра K
        K_ = P_ * Ht * arma::inv( S_ );
#ifdef DEBUG_KALMAN
        K_.print( this->filterName_ + " Correction, K:" );
#endif
        // 3. Вычисление невязки Delta
        if( !deltaY_isSet ) {
            assert( Y_msd_isSet ); // Если не установлен deltaY, то Y_msd_ обязан быть установлен
            Y_msd_isSet = false; // Сразу же снять признак
            DeltaY_ = Y_msd_ - Y_pred_;
            if( checkDeltaMeasurement_ != nullptr ) {
                DeltaY_ = checkDeltaMeasurement_( DeltaY_ );
            }
        }
#ifdef DEBUG_KALMAN
        DeltaY_.print( this->filterName_ + " Correction, DeltaY:" );
#endif
        // 4. Вычисление ковариационной матрицы P
        P_ = ( I_ - ( K_ * H_ ) ) * P_;
        fixMatrixMainDiagonalSymmetry( P_ );        
#ifdef DEBUG_KALMAN
        P_.print( this->filterName_ + " Correction, P after:" );
#endif
        checkMatrixDiagPositive( P_ );

        // 5. Вычисление X_est, Y_est
        X_est_ = X_pred_ + ( K_ * DeltaY_ );
        if( checkBordersStateAfterCorrection_ != nullptr ) {
            X_est_ = checkBordersStateAfterCorrection_( X_est_ ); // Проверка вектора состояния X на выход за допустимые пределы
        }
        Y_est_ = H_ * X_est_;
        if( checkBordersMeasurement_ != nullptr ) {
            Y_est_ = checkBordersMeasurement_( Y_est_ );
        }
#ifdef DEBUG_KALMAN
        X_est_.print( this->filterName_ + " Correction, X_est:" );
        Y_est_.print( this->filterName_ + " Correction, Y_est:" );
#endif
        deltaY_isSet = false; // Снять признак (выставляется в true в сеттере deltaY)
    }

    ///
    /// \brief Отдельное вычисление ковариационной матрицы S невязки измерений
    /// \details Отдельное вычисление S применяется при стробировании
    /// \attention Необходимо учесть 2 обстоятельства:
    /// 1) метод должен выполняться после вызова метода Prediction, где должна быть вычислена матрица H;
    /// 2) матрица R при вычислении S будет использована та, которая уже имеется в фильтре и если необходима иная
    /// матрица R, то следует её установить, вызвав метод SetObservationCovarianceMatrixR
    /// \param PdiagAdd - добавка в диагональ матрицы P, размерность [SizeX]
    /// \param Rdiag - диагональная матрица R априорных шумов измерений, размерность [SizeY]
    ///
    virtual void CalculateInnovationCovarianceS( const arma::vec &PdiagAdd, const arma::vec Rdiag )
    {        
        assert( prediction_isDone ); // Перед вы обязательно должен быть выполнен прогноз, иначе не имеет смысла:
        arma::mat Ptmp = P_;
        if( PdiagAdd.size() != 0 ) { // Если размер не нулевой - прибавим,
            assert( PdiagAdd.size() == SizeX_ ); // проверяя при этом, что размеры совпадают!
            Ptmp.diag() += PdiagAdd;
        }        
//        S_ = ( H_ * Ptmp * arma::trans( H_ ) ) + R_; // Сразу прибавить R (в случае, если матрица R - плотная)
        S_ = ( H_ * Ptmp * arma::trans( H_ ) );
        S_.diag() += Rdiag;
        fixMatrixMainDiagonalSymmetry( S_ );
        checkMatrixDiagPositive( S_ );
    }

protected:
    //------------------------------------------------------------------------------------------------------------------
#ifdef DEBUG_KALMAN
    std::string filterName_; ///< Название фильтра
#endif

    // Размерности:
    size_t SizeX_;      ///< Размерность вектора состояния Х (state)
    size_t SizeY_;      ///< Размернсоть вектора измУстановка матрицы SetEstimateCovarianceMatrixPерений Y (measurement)

    // Матрицы системы:
    arma::mat F_;       ///< Матрица эволюции системы (перехода состояния) (state-transition model), размерность [SizeX * SizeX]
    arma::mat H_;       ///< Матрица измерений (перехода измерений) (observation model), размерность [SizeY * SizeX]

    // Дополнительные матрицы:
    arma::mat K_;       ///< Коэффициент усиления фильтра Калмана (Kalman gain), размерность [SizeX * SizeY]
    arma::mat I_;       ///< Единичная матрица, размерность [SizeX * SizeX]
    arma::vec DeltaY_;  ///< Вектор невязки измерений, размерность [SizeY]

    // Ковариационные матрицы:
    arma::mat P_;       ///< Ковариационная матрица вектора состояния X (estimate covariance matrix), размерность [SizeX * SizeX]
    arma::mat S_;       ///< Ковариационая матрица вектора невязки (innovation covariance), размерность [SizeY * SizeY]
    arma::mat Q_;       ///< Ковариационая матрица (обычно диагональная) шумов вектора состояния Х НА 1 СЕКУНДЕ (covariance of the process noise), размерность [SizeX * SizeX]
    arma::mat R_;       ///< Ковариационая матрица (обычно диагональная) шумов вектора измерений Y (covariance of the observation noise), размерность [SizeY * SizeY]

    // Вектора состояния и измерения:
    arma::vec X_pred_;  ///< Экстраполированный (predicted) вектор состояния X, размерность [SizeX]
    arma::vec Y_pred_;  ///< Экстраполированный (predicted) вектор измерений Y, размерность [SizeY]
    arma::vec X_est_;   ///< Скорректированный (estimated) вектор состояния X, размерность [SizeX]
    arma::vec Y_est_;   ///< Скорректированный (estimated) вектор измерений Y, размерность [SizeY]
    arma::vec Y_msd_;   ///< Измеренный (measured) вектор Y (отметка), размерность [SizeY]

    // Защитные признаки:
    bool deltaY_isSet;  ///< Признак установки вектора невязки DeltaY (нужно в методе Correction)
    bool Y_msd_isSet;   ///< Признак установки вектора измерений Y_msd;
    bool prediction_isDone; ///< Признак выполненного прогноза (без этого нельзя делать фильтрацию)

    // Обертки функций проверок (могут быть переданы в класс при необходимости)
    std::function<arma::vec( const arma::vec &X )> checkBordersStateAfterPrediction_; ///< Проверка границ вектора состояния X после прогноза
    std::function<arma::vec( const arma::vec &X )> checkBordersStateAfterCorrection_; ///< Проверка границ вектора состояния X после фильтрации
    std::function<arma::vec( const arma::vec &Y )> checkBordersMeasurement_; ///< Проверка границ вектора измерения Y
    std::function<arma::vec( const arma::vec &DeltaX )> checkDeltaState_; ///< Проверка разности векторов состояния X
    std::function<arma::vec( const arma::vec &DeltaY )> checkDeltaMeasurement_; ///< Проверка разности векторов измерения Y

    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Исправление симметричности матрицы относительно главной диагонали
    /// \details Элементам вне главной диагонали присваивается полусумма между соответствующими элементами
    /// \param A - проверяемая матрица
    ///
    void fixMatrixMainDiagonalSymmetry( arma::mat &A )
    {
        if( !A.is_square() ) {
            throw std::length_error( "Matrix is not square" );
        } else {
            double new_value = 0.0;
            // Пройдемся только по треугольнику над главной диагональю
            for( unsigned long long j = 1; j < A.n_cols; j++ ) {
                for( unsigned long long i = 0; i < j; i++ ) {
                    new_value = ( A(i,j) + A(j,i) ) * 0.5;
                    A(i,j) = new_value;
                    A(j,i) = new_value;
                }
            }
        }
    }

    ///
    /// \brief Проверка что в диагонали матрицы лежат только положительные элементы
    /// \param A - проверяемая матрица
    /// \attention Метод вызывает assert!
    ///
    void checkMatrixDiagPositive( const arma::mat &A ) const
    {
        if( !A.is_square() ) {
            throw std::length_error( "Matrix is not square" );
        } else {
            for( unsigned long long i = 0; i < A.n_cols; i++ ) {
                assert( A( i, i ) > 0.0 );
            }
        }
    }

private:
    //------------------------------------------------------------------------------------------------------------------
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
