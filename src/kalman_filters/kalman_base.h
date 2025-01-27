//------------------------------------------------------------------------------
///
/// \file       kalman_filter_base.h
/// \brief      Базовый шаблонный класс фильтра Калмана
/// \date       08.11.23 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_BASE_H
#define KALMAN_FILTER_BASE_H

// System includes:
#include <armadillo> // Матрицы
#include <cassert> // assert
#include <functional> // std::function

// Project includes
#include <kalman_debug.h>

namespace KalmanFilters /// Фильтры Калмана
{
//------------------------------------------------------------------------------
///
/// \brief Тип фильтра
/// \details Расширенный, сигма-точечный, кубатурный
///
enum TKalmanType : int
{    
    KT_Extended = 1,
    KT_Unscented = 2,
    KT_Cubature = 3    
};

//------------------------------------------------------------------------------
///
/// \brief Тип ковариационных матриц фильтра (полные/квадратно-корневые)
///
enum TCovType : int
{
    CT_FullMatrices = 0, ///< Фильтр в полных ковариационных матрицах
    CT_SquareRootMatrices = 1 ///< Фильтр в квадратно-корневых ковариационных матрицах
};

//------------------------------------------------------------------------------
///
/// \brief Базовый шаблонный класс фильтра Калмана
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanBase
{
public:        
    //--------------------------------------------------------------------------
    // Конструкторы:

    CKalmanBase() = default;
    CKalmanBase( const CKalmanBase& ) = default;
    CKalmanBase& operator=( const CKalmanBase& ) = default;
    CKalmanBase( CKalmanBase&& ) = default;
    CKalmanBase& operator=( CKalmanBase&& ) = default;
    virtual ~CKalmanBase() = default;

    ///
    /// \brief Название фильтра
    ///
    virtual std::string GetFilterName() const = 0;

    ///
    /// \brief Тип фильтра (расширенный/сигма-точечный/кубатурный)
    ///
    virtual TKalmanType GetKalmanType() const = 0;

    ///
    /// \brief Тип ковариационных матриц фильтра (полные/квадратно-корневые)
    ///
    virtual TCovType GetCovType() const = 0;

    //--------------------------------------------------------------------------
    // Методы-сеттеры:

    virtual void SetDesignParametersMeanSet( double w0 ) = 0;
    virtual void SetDesignParametersScaledSet( double alpha, double beta, double kappa ) = 0;
    virtual void SetDesignParametersCubatureBaseSet() = 0;

    ///
    /// \brief Установка функции прогноза состояния (predictState)
    /// \sa stateTransitionModel_
    ///
    void SetStateTransitionModel( std::function<arma::vec( const arma::vec &X, double dt )> stateTransitionModel )
    {
        this->stateTransitionModel_ = stateTransitionModel;
    }

    ///
    /// \brief Установка функции перевода состояния в измерение (XtoY)
    /// \sa observationModel_
    ///
    void SetObservationModel( std::function<arma::vec( const arma::vec &X )> observationModel )
    {
        this->observationModel_ = observationModel;
    }

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

    //--------------------------------------------------------------------------
    ///
    /// \brief Установка функции вычисления взвешенной суммы сигма-точек пространства Х
    /// \param weightedSumStateSigmas - Функция вычисления взвешенной суммы сигма-точек пространства Х
    ///
    void SetWeightedSumStateSigmas( std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumStateSigmas )
    {
        weightedSumStateSigmas_ = weightedSumStateSigmas;
    }

    ///
    /// \brief Установка функции вычисления взвешенной суммы сигма-точек пространства Y
    /// \param weightedSumMeasurementSigmas - Функция вычисления взвешенной суммы сигма-точек пространства Y
    ///
    void SetWeightedSumMeasurementSigmas( std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumMeasurementSigmas )
    {
        weightedSumMeasurementSigmas_ = weightedSumMeasurementSigmas;
    }

    //--------------------------------------------------------------------------
    // Методы-геттеры:

    // ///
    // /// \brief Получить название фильтра
    // /// \return Строка с именем фильтра
    // ///
    // const std::string GetFilterName() const { return filterName_; }

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

    //--------------------------------------------------------------------------
    // Отдельно геттеры/сеттеры для X и Y
    arma::vec &X() { return X_est_; }
    const arma::vec &X() const { return X_est_; }
    double &X( int index ) { return X_est_( index ); }
    const double &X( int index ) const { return X_est_( index ); }
    //---
    arma::vec &Y() { return Y_est_; }
    const arma::vec &Y() const { return Y_est_; }
    double &Y( int index ) { return Y_est_( index ); }
    const double &Y( int index ) const { return Y_est_( index ); }
    //---
    arma::mat &P() { return P_; }
    const arma::mat &P() const { return P_; }
    double &P( int i, int j ) { return P_( i, j ); }
    const double &P( int i, int j ) const { return P_( i, j ); }
    //---
//    arma::mat &S() { return S_; }
    const arma::mat &S() const { return S_; }
    double &S( int i, int j ) { return S_( i, j ); }
    const double &S( int i, int j ) const { return S_( i, j ); }
    //---
//    arma::mat &H() { return H_; }
    const arma::mat &H() const { return H_; }
    double &H( int i, int j ) { return H_( i, j ); }
    const double &H( int i, int j ) const { return H_( i, j ); }
    //---
//    arma::mat &Q() { return Q_; }
    const arma::mat &Q() const { return Q_; }
    double &Q( int i, int j ) { return Q_( i, j ); }
    const double &Q( int i, int j ) const { return Q_( i, j ); }

    //--------------------------------------------------------------------------
    // Виртуальные методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt ) = 0;

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd ) = 0;
    
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
    // virtual void CalculateInnovationCovarianceS( const arma::vec &PdiagAdd, const arma::vec Rdiag ) = 0;
    virtual void CalculateInnovationCovarianceS( 
        const arma::vec &PdiagAdd, const arma::vec Rdiag )
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
    
protected:
    //--------------------------------------------------------------------------  
    // Размерности:
    const static size_t SizeX_ = SizeX; ///< Размерность вектора состояния Х (state)
    const static size_t SizeY_ = SizeY; ///< Размернсоть вектора измУстановка матрицы SetEstimateCovarianceMatrixPерений Y (measurement)

    // Матрицы системы:
    arma::mat::fixed<SizeX, SizeX> F_ = arma::mat::fixed<SizeX, SizeX>( arma::fill::eye ); ///< Матрица эволюции системы (перехода состояния) (state-transition model), размерность [SizeX * SizeX]
    arma::mat::fixed<SizeY, SizeX> H_; ///< Матрица измерений (перехода измерений) (observation model), размерность [SizeY * SizeX]

    // Дополнительные матрицы:
    arma::mat::fixed<SizeX, SizeY> K_; ///< Коэффициент усиления фильтра Калмана (Kalman gain), размерность [SizeX * SizeY]
    const arma::mat::fixed<SizeX, SizeX> I_ = arma::mat::fixed<SizeX, SizeX>( arma::fill::eye ); ///< Единичная матрица, размерность [SizeX * SizeX]
    arma::vec::fixed<SizeY> DeltaY_; ///< Вектор невязки измерений, размерность [SizeY]

    // Ковариационные матрицы:
    arma::mat::fixed<SizeX, SizeX> P_; ///< Ковариационная матрица вектора состояния X (estimate covariance matrix), размерность [SizeX * SizeX]
    arma::mat::fixed<SizeY, SizeY> S_; ///< Ковариационая матрица вектора невязки (innovation covariance), размерность [SizeY * SizeY]
    arma::mat::fixed<SizeX, SizeX> Q_; ///< Ковариационая матрица (обычно диагональная) шумов вектора состояния Х НА 1 СЕКУНДЕ (covariance of the process noise), размерность [SizeX * SizeX]
    arma::mat::fixed<SizeY, SizeY> R_; ///< Ковариационая матрица (обычно диагональная) шумов вектора измерений Y (covariance of the observation noise), размерность [SizeY * SizeY]

    // Вектора состояния и измерения:
    arma::vec::fixed<SizeX> X_pred_; ///< Экстраполированный (predicted) вектор состояния X, размерность [SizeX]
    arma::vec::fixed<SizeY> Y_pred_; ///< Экстраполированный (predicted) вектор измерений Y, размерность [SizeY]
    arma::vec::fixed<SizeX> X_est_; ///< Скорректированный (estimated) вектор состояния X, размерность [SizeX]
    arma::vec::fixed<SizeY> Y_est_; ///< Скорректированный (estimated) вектор измерений Y, размерность [SizeY]
    arma::vec::fixed<SizeY> Y_msd_; ///< Измеренный (measured) вектор Y (отметка), размерность [SizeY]

    // Защитные признаки:
    bool deltaY_isSet = false; ///< Признак установки вектора невязки DeltaY (нужно в методе Correction)
    bool Y_msd_isSet = false; ///< Признак установки вектора измерений Y_msd;
    bool prediction_isDone = false; ///< Признак выполненного прогноза (без этого нельзя делать фильтрацию)

    // Обертки функций проверок (могут быть переданы в класс при необходимости)
    std::function<arma::vec( const arma::vec &X )> checkBordersStateAfterPrediction_ = nullptr; ///< Проверка границ вектора состояния X после прогноза
    std::function<arma::vec( const arma::vec &X )> checkBordersStateAfterCorrection_ = nullptr; ///< Проверка границ вектора состояния X после фильтрации
    std::function<arma::vec( const arma::vec &Y )> checkBordersMeasurement_ = nullptr; ///< Проверка границ вектора измерения Y
    std::function<arma::vec( const arma::vec &DeltaX )> checkDeltaState_ = nullptr; ///< Проверка разности векторов состояния X
    std::function<arma::vec( const arma::vec &DeltaY )> checkDeltaMeasurement_ = nullptr; ///< Проверка разности векторов измерения Y

    //------------------------------------------------------------------------------------------------------------------
    // Обертки функций EKF
    ///
    /// \brief Функция прогноза состояния (predictState)
    /// \param X - вектор состояния с прошлого момента времени
    /// \sa SetStateTransitionModel
    ///
    std::function<arma::vec( const arma::vec &X, double dt )> stateTransitionModel_ = nullptr;

    ///
    /// \brief Функция перевода состояния в измерения (XtoY)
    /// \param X - вектор состояния текущего момента времени
    /// \sa SetObservationModel
    ///
    std::function<arma::vec( const arma::vec &X )> observationModel_ = nullptr;

    ///
    /// \brief Функция вычисления матрицы перехода состояния F (makeMatrixF)
    /// \param X - вектор состояния с прошлого момента времени
    /// \sa SetStateTransitionJacobianF
    ///
    std::function<arma::mat( const arma::vec &X, double dt )> stateTransitionJacobianF_ = nullptr;

    ///
    /// \brief Функция вычисления матрицы перехода измерений H (makeMatrixH)
    /// \param X - вектор состояния текущего момента времени
    /// \sa SetObservationJacobianH
    ///
    std::function<arma::mat( const arma::vec &X )> observationJacobianH_ = nullptr;

    //--------------------------------------------------------------------------
    // Обертки функций вычисления взвешенной суммы:
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumStateSigmas_ = nullptr; ///< Вычисление взвешенной суммы сигма-точек пространства Х
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumMeasurementSigmas_ = nullptr; ///< Вычисление взвешенной суммы сигма-точек пространства Y

    //--------------------------------------------------------------------------
    ///
    /// \brief Исправление симметричности полной (не нижне/верхне треугольной, естественно!) матрицы относительно главной диагонали
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
            for( arma::uword j = 1; j < A.n_cols; j++ ) {
                for( arma::uword i = 0; i < j; i++ ) {
                    new_value = ( A(i,j) + A(j,i) ) * 0.5;
                    A(i,j) = new_value;
                    A(j,i) = new_value;
                }
            }
        }
    }

    ///
    /// \brief Возвращает индекс первого найденного отрицательного диагонального элемента при наличии
    /// \param A - проверяемая матрица
    /// \return -1 если отрицательных элементов в диагонали нетб иначе - индекс первого найденного (начиная с нуля) отрицательного элемента
    ///
    int indexOfNegativeDiagonalElement( const arma::mat &A ) const
    {
        for( arma::uword i = 0; i < A.n_cols; i++ ) {
            if( A( i, i ) < 0.0 ) {
                return i;
            }
        }
        return -1;
    }
};

}

#endif // KALMAN_FILTER_BASE_H
/// \}
