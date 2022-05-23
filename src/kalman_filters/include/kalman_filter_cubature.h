//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_cubature.h
/// \brief      Шаблонный класс кубатурного фильтра Калмана, КФК (Cubature Kalman Filter, CKF)
/// \date       24.03.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_CUBATURE_H
#define KALMAN_FILTER_CUBATURE_H

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_extended.h>
#include <qr_decomposition.h> // QR разложение матрицы

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Шаблонный класс кубатурного фильтра Калмана, КФК (Cubature Kalman Filter, CKF)
/// \details Частный случай сигма-точечного фильтра Калмана при параметрах разброса сигма-точек выбранных
/// по рекомендации Merwe alpha=1.0, beta=0.0, kappa=0.0.
/// \details Источники:
/// \n[1] Cubature Kalman Filters, Ienkaran Arasaratnam and Simon Haykin, Life Fellow, IEEE
/// \n[2] Sebastian Bitzer, Technische Universität Dresden, https://github.com/sbitzer/UKF-exposed/blob/master/UKF.pdf
/// \attention Фильтр построен по классическому его варианту НИЖНЕтреугольного разложения Холецкого!
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanCKF : virtual public CKalmanEKF<SizeX, SizeY>
{
public:
    using CKalmanEKF<SizeX, SizeY>::CKalmanEKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanCKF() : CKalmanEKF<SizeX, SizeY>(),
        k_sigma_points_( 2 * SizeX ),
        weights_mean_( arma::zeros( this->k_sigma_points_ ) ),
        weights_covariance_( arma::zeros( this->k_sigma_points_ ) ),
        x_est_sigma_points_( arma::zeros( SizeX, this->k_sigma_points_ ) ),
        x_pred_sigma_points_( arma::zeros( SizeX, this->k_sigma_points_ ) ),
        y_pred_sigma_points_( arma::zeros( SizeY, this->k_sigma_points_ ) ),

        dXcal_( arma::zeros( SizeX, this->k_sigma_points_ ) ),
        dYcal_( arma::zeros( SizeY, this->k_sigma_points_ ) ),
        P_xy_( arma::zeros( SizeX, SizeY ) ),
        sqrt_P_chol_( arma::zeros( SizeX, SizeX ) )
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "CKF" );
#endif
        this->SetupDesignParametersCubatureBaseSet();
    }

    ///
    /// \brief Конструктор копирования
    /// \param other - экземпляр, с которого делается копия
    ///
    CKalmanCKF( const CKalmanCKF &other ) : CKalmanEKF<SizeX, SizeY>( other )
    {
        this->k_sigma_points_ = other.k_sigma_points_;
        this->weights_mean_ = other.weights_mean_;
        this->weights_covariance_ = other.weights_covariance_;
        this->x_est_sigma_points_ = other.x_est_sigma_points_;
        this->x_pred_sigma_points_ = other.x_pred_sigma_points_;
        this->y_pred_sigma_points_ = other.y_pred_sigma_points_;

        this->dXcal_ = other.dXcal_;
        this->dYcal_ = other.dYcal_;
        this->P_xy_ = other.P_xy_;
        this->sqrt_P_chol_ = other.sqrt_P_chol_;
    }

    ///
    /// \brief Перегрузка оператора присвоения
    /// \param other - экземпляр, с которого делается копия
    /// \return *this
    ///
    CKalmanCKF& operator=( const CKalmanCKF &other )
    {
        CKalmanCKF copy( other );
        swap( *this, copy );
        return *this;
    }

    ///
    /// \brief Конструктор перемещения
    /// \param other - экземпляр, с которого делается копия
    ///
    CKalmanCKF( CKalmanCKF &&other ) noexcept
    {
        swap( *this, other );
    }

    ///
    /// \brief Перегрузка оператора перемещения
    /// \param other - экземпляр, с которого делается копия
    /// \return *this
    ///
    CKalmanCKF& operator=( CKalmanCKF &&other ) noexcept
    {
        swap( *this, other );
        return *this;
    }

    virtual ~CKalmanCKF() = default; ///< Дестркутор

    ///
    /// \brief Метод свапа
    /// \param lhs - left hand side instance
    /// \param rhs - right hand side instance
    ///
    friend void swap( CKalmanCKF<SizeX, SizeY> &lhs, CKalmanCKF<SizeX, SizeY> &rhs ) noexcept
    {
        swap( dynamic_cast< CKalmanEKF<SizeX, SizeY> &>( lhs ), dynamic_cast< CKalmanEKF<SizeX, SizeY> &>( rhs ) ); // Свап предка 1 порядка

        std::swap( lhs.k_sigma_points_, rhs.k_sigma_points_ );
        std::swap( lhs.weights_mean_, rhs.weights_mean_ );
        std::swap( lhs.weights_covariance_, rhs.weights_covariance_ );
        std::swap( lhs.x_est_sigma_points_, rhs.x_est_sigma_points_ );
        std::swap( lhs.x_pred_sigma_points_, rhs.x_pred_sigma_points_ );
        std::swap( lhs.y_pred_sigma_points_, rhs.y_pred_sigma_points_ );

        std::swap( lhs.dXcal_, rhs.dXcal_ );
        std::swap( lhs.dYcal_, rhs.dYcal_ );
        std::swap( lhs.P_xy_, rhs.P_xy_ );
        std::swap( lhs.sqrt_P_chol_, rhs.sqrt_P_chol_ );
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы-сеттеры:

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

    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {        
        PredictionCKF( dt );
    }

    ///
    /// \brief Коррекция
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
        CorrectionCKF( Y_msd );
    }

protected:
    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Прогноз CKF
    /// \param dt - Время прогноза, [с]
    ///
    void PredictionCKF( double dt )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Prediction started, dt = " << dt << std::endl;
#endif
        assert( dt > 0.0 ); // Прогноз имеет смысл только при положительном времени

        // 1. Создание сигма-точек пространства X
        this->sqrt_P_chol_ = arma::chol( this->P_, "lower" ); // Должно быть взято НИЖНЕЕ РАЗЛОЖЕНИЕ!
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Prediction, P before:" );
        ( this->sqrt_P_chol_ ).print( this->filterName_ + " Prediction, sqrt_P_chol:" );
#endif
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->sqrt_P_chol_.col(i) );
            this->x_est_sigma_points_.col(i) = this->X_est_ + add;
            this->x_est_sigma_points_.col(i + SizeX) = this->X_est_ - add;
            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
                this->x_est_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_est_sigma_points_.col(i) );
                this->x_est_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_est_sigma_points_.col(i + SizeX) );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->x_est_sigma_points_ ).print( this->filterName_ + " Prediction, x_est_sigma_points_:" );
#endif
        // 2. Прогноз сигма-точек пространства X на текущий такт
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->x_pred_sigma_points_.col(i) = this->stateTransitionModel_( this->x_est_sigma_points_.col(i), dt );
        }
#ifdef DEBUG_KALMAN
        ( this->x_pred_sigma_points_ ).print( this->filterName_ + " Prediction, x_pred_sigma_points_:" );
#endif
        // 3. Вычисление X_pred по сигма-точкам пространства Х
        if( this->weightedSumStateSigmas_ == nullptr ) {
            this->X_pred_ = this->x_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->X_pred_ = this->weightedSumStateSigmas_( this->weights_mean_, this->x_pred_sigma_points_ );
        }
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }
#ifdef DEBUG_KALMAN
        ( this->X_pred_ ).print( this->filterName_ + " Prediction, X_pred_:" );
#endif
        // 4. Вычисление ковариационной матрицы Р
        this->P_ = ( this->Q_ * dt );
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = this->x_pred_sigma_points_.col(i) - this->X_pred_;
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->P_ += this->weights_covariance_(i) * this->dXcal_.col(i) * arma::trans( this->dXcal_.col(i) );
        }
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->dXcal_ ).print( this->filterName_ + " Prediction, dXcal_:" );
        ( this->P_ ).print( this->filterName_ + " Prediction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );
        // 5. Пересоздание сигма-точек после прогноза по новой P
        this->sqrt_P_chol_ = arma::chol( this->P_, "lower" ); // Должно быть взято НИЖНЕЕ РАЗЛОЖЕНИЕ!
        for( size_t i = 0; i < SizeX; i++ ) {
            arma::vec add = ( this->gamma_ * this->sqrt_P_chol_.col(i) );
            this->x_pred_sigma_points_.col(i) = this->X_pred_ + add;
            this->x_pred_sigma_points_.col(i + SizeX) = this->X_pred_ - add;
            if( this->checkBordersStateAfterPrediction_ != nullptr ) {
                this->x_pred_sigma_points_.col(i) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i) );
                this->x_pred_sigma_points_.col(i + SizeX) = this->checkBordersStateAfterPrediction_( this->x_pred_sigma_points_.col(i + SizeX) );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->x_pred_sigma_points_ ).print( this->filterName_ + " Prediction, x_pred_sigma_points_:" );
#endif
        // 6. Вычисление X_pred заново по пересозданным сигма-точкам пространства Х (в [1] не указано, но это подразумевается)
        if( this->weightedSumStateSigmas_ == nullptr ) {
            this->X_pred_ = this->x_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->X_pred_ = this->weightedSumStateSigmas_( this->weights_mean_, this->x_pred_sigma_points_ );
        }
        if( this->checkBordersStateAfterPrediction_ != nullptr ) {
            this->X_pred_ = this->checkBordersStateAfterPrediction_( this->X_pred_ );
        }

        // 7. Вычисление сигма-точек пространства Y
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }
#ifdef DEBUG_KALMAN
        ( this->y_pred_sigma_points_ ).print( this->filterName_ + " Prediction, y_pred_sigma_points_:" );
#endif
        // 8. Вычисление предсказанного вектора Y_pred на текущий такт по сигма-точкам пространства Х
        if( this->weightedSumMeasurementSigmas_ == nullptr ) {
            this->Y_pred_ = this->y_pred_sigma_points_ * this->weights_mean_; // В матричной форме
        } else {
            this->Y_pred_ = this->weightedSumMeasurementSigmas_( this->weights_mean_, this->y_pred_sigma_points_ );
        }
        if( this->checkBordersMeasurement_ != nullptr ) {
            this->Y_pred_ = this->checkBordersMeasurement_( this->Y_pred_ );
        }

        // 9. Обновление X_est, Y_est
        this->X_est_ = this->X_pred_;
        this->Y_est_ = this->Y_pred_;
#ifdef DEBUG_KALMAN
        ( this->X_est_ ).print( this->filterName_ + " Prediction, X_est:" );
        ( this->Y_est_ ).print( this->filterName_ + " Prediction, Y_est:" );
#endif
        this->prediction_isDone = true; // Выставить признак состоявшегося прогноза
    }

    ///
    /// \brief Коррекция CKF
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    void CorrectionCKF( const arma::vec &Y_msd )
    {
#ifdef DEBUG_KALMAN
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << this->filterName_ + " Correction started" << std::endl;
        Y_msd.print( this->filterName_ + " Correction, Y_msd" );
#endif
        this->SetMeasuredVectorY( Y_msd );

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла
        this->prediction_isDone = false; // Сразу же снять признак

        // 1. Вычисление ковариационной матрицы S (P_yy)
        this->S_ = this->R_;
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dYcal_.col(i) = this->y_pred_sigma_points_.col(i) - this->Y_pred_;
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->dYcal_.col(i) = this->checkDeltaMeasurement_( this->dYcal_.col(i) );
            }
            this->S_ += this->weights_covariance_(i) * this->dYcal_.col(i) * arma::trans( this->dYcal_.col(i) );
        }
        this->fixMatrixMainDiagonalSymmetry( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( this->filterName_ + " Correction, S:" );
#endif
        this->checkMatrixDiagPositive( this->S_ );
        // 2. Вычисление кросс-ковариационной матрицы P_xy
        this->P_xy_.zeros();
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = this->x_pred_sigma_points_.col(i) - this->X_pred_;
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            // dYcal_ - вычислено выше в 1, поэтому заново - не надо
            this->P_xy_ += this->weights_covariance_(i) * this->dXcal_.col(i) * arma::trans( this->dYcal_.col(i) );
        }
#ifdef DEBUG_KALMAN
        ( this->dXcal_ ).print( this->filterName_ + " Prediction, dXcal_:" );
        ( this->dYcal_ ).print( this->filterName_ + " Prediction, dYcal_:" );
        ( this->P_xy_ ).print( this->filterName_ + " Correction, P_xy:" );
#endif
        // 3. Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_xy_ * arma::inv( this->S_ );
#ifdef DEBUG_KALMAN
        ( this->K_ ).print( this->filterName_ + " Correction, K:" );
#endif
        // 4. Вычисление невязки Delta
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
        // 5. Вычисление ковариационной матрицы P
        this->P_ = this->P_ - ( this->K_ * this->S_ * arma::trans( this->K_ ) );
//        this->P_ = this->P_ - P_xy * arma::trans( this->K_ );
        this->fixMatrixMainDiagonalSymmetry( this->P_ );
#ifdef DEBUG_KALMAN
        ( this->P_ ).print( this->filterName_ + " Correction, P after:" );
#endif
        this->checkMatrixDiagPositive( this->P_ );
        // 6. Вычисление X_est, Y_est
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

    //------------------------------------------------------------------------------------------------------------------
    // Параметры, зависящие от SizeX, SizeY:
    int k_sigma_points_;            ///< Число сигма-точек
    arma::vec weights_mean_;        ///< Веса среднего
    arma::vec weights_covariance_;  ///< Веса ковариации
    arma::mat x_est_sigma_points_;  ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве X на текущем такте, размерность [SizeX,k_sigma_points_]
    arma::mat x_pred_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве X, экстраполированный на текущий такт, размерность [SizeX,k_sigma_points_]
    arma::mat y_pred_sigma_points_; ///< Матрица сигма-точек (сигма-точки - столбцы) в пространстве Y, экстраполированный на текущий такт, размерность [SizeX,k_sigma_points_]

    arma::mat dXcal_;       ///< Матрица Х-каллиграфическое (матрица сигма-точек - столбцов)
    arma::mat dYcal_;       ///< Матрица Y-каллиграфическое (матрица сигма-точек - столбцов)
    arma::mat P_xy_;        ///< Матрица кросс-коварации векторов Х и Y, размерность [SizeX * SizeY]
    arma::mat sqrt_P_chol_; ///< Корень из матрицы P
    //------------------------------------------------------------------------------------------------------------------
    // Параметры выбора сигма-точек, определяемые в методах SetDesignParameters*
    double gamma_;  ///< Автоматически вычисляемый (в методах SetDesignParameters*) параметр (множитель при корне из P при создании сигма-точек)

    //------------------------------------------------------------------------------------------------------------------
    // Обертки функций вычисления взвешенной суммы:
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumStateSigmas_; ///< Вычисление взвешенной суммы сигма-точек пространства Х
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumMeasurementSigmas_; ///< Вычисление взвешенной суммы сигма-точек пространства Y

    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Установка кубатурных весов (базовый вариант ансцентного преобразования)
    /// \details Смотри [1] и [2]
    ///
    void SetupDesignParametersCubatureBaseSet()
    {
        double gammaSq = static_cast<double>( SizeX );
        this->gamma_ = std::sqrt( gammaSq );

        // Set the weights for sigma points
        double all_points = 1.0 / ( 2.0 * gammaSq );
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->weights_mean_( i ) = all_points;
            this->weights_covariance_( i ) = all_points;
        }

#ifdef DEBUG_KALMAN
        std::cout << this->filterName_ + " SetDesignParameters_Cubature:" << std::endl;
        std::cout << "weights_mean_:" << std::endl;
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            std::cout << this->weights_mean_( i ) << std::endl;
        }
        std::cout << "weights_covariance_:" << std::endl;
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            std::cout << this->weights_covariance_( i ) << std::endl;
        }
#endif
    }

    //------------------------------------------------------------------------------------------------------------------
    // Запрет доступа к методам родительского класса:
    void SetStateTransitionJacobianF(){} ///< Запрет доступа
    void SetObservationJacobianH(){} ///< Запрет доступа
};

}

#endif // KALMAN_FILTER_CUBATURE_H
/// \}
