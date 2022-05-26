//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_extended_unscented_square_root.h
/// \brief      Шаблонный класс квадратно-корневого расширенного сигма-точечного (ансцентного) фильтра Калмана, КК-РСТФК
///             (Square Root Extended Unscented Kalman Filter, SR-EUKF)
/// \date       11.04.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_EXTENDED_UNSCENTED_SQUARE_ROOT_H
#define KALMAN_FILTER_EXTENDED_UNSCENTED_SQUARE_ROOT_H

// Project includes:
#include <kalman_filter_debug.h>
#include <kalman_filter_extended_square_root.h>
#include <kalman_filter_unscented_square_root.h>
#include <qr_decomposition.h> // JQR разложение матрицы

namespace KalmanFilters /// Фильтры Калмана
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief  Шаблонный класс квадратно-корневого расширенного сигма-точечного (ансцентного) фильтра Калмана, КК-РСТФК
///         (Square Root Extended Unscented Kalman Filter, SR-EUKF)
/// \details Источники:
/// \n[1] A New Extension of the Kalman Filter to Nonlinear Systems, Simon J. Julier, Jeffrey K. Uhlmann,
/// The Robotics Research Group, Department of Engineering Science, The University of Oxford, 1997
/// \n[2] The Unscented Kalman Filter for Nonlinear Estimation, Eric A. Wan and Rudolph van der Merwe
/// Oregon Graduate Institute of Science & Technology 20000 NW Walker Rd, Beaverton, Oregon 97006, 2000
/// ericwan[at]ece.ogi.edu, rvdmerwe[at]ece.ogi.edu
/// \n[3] Sigma-Point Kalman Filters for Probabilistic Inference in Dynamic State-Space Models, Rudolph van der Merwe &
/// Eric Wan, OGI School of Science & Engineering Oregon Health & Science University Beaverton, Oregon, 97006, USA, 2003
/// {rvdmerwe,ericwan}[at]ece.ogi.edu
/// \n[4] THE SQUARE-ROOT UNSCENTED KALMAN FILTER FOR STATE AND PARAMETER-ESTIMATION, Rudolph van der Merwe and
/// Eric A. Wan, Oregon Graduate Institute of Science and Technology 20000 NW Walker Road, Beaverton, Oregon 97006, USA
/// rvdmerwe,ericwan[at]ece.ogi.edu
/// \attention Внимание! В реализации UKF вес нулевой сигма-точки Wcov не может быть отрицательным! (Wmean - может)
/// \attention Фильтр построен по классическому его варианту НИЖНЕтреугольного разложения Холецкого!
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSREUKF : public CKalmanSREKF<SizeX, SizeY>, public CKalmanSRUKF<SizeX, SizeY>
{
public:
    using CKalmanSREKF<SizeX, SizeY>::CKalmanSREKF;
    using CKalmanSRUKF<SizeX, SizeY>::CKalmanSRUKF;
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSREUKF() : CKalmanSREKF<SizeX, SizeY>(), CKalmanSRUKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SREUKF" );
#endif
        createSignMatrices();
    }

    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    ///
    /// \brief Прогноз квадратно-корневого расширенного сигма-точечного фильтра Калмана (КК-РСТФК, SR-EUKF)
    /// \param dt - Время прогноза, [с]
    ///
    virtual void Prediction( double dt )
    {
        this->PredictionSREKF( dt );
    }

    ///
    /// \brief Коррекция квадратно-корневого расширенного сигма-точечного фильтра Калмана (КК-РСТФК, SR-EUKF)
    /// \param Y_msd - вектор измерений, по которым производится коррекция
    ///
    virtual void Correction( const arma::vec &Y_msd )
    {
        // 1a. Создание сигма-точек пространства X                
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
        // 1b. Вычисление сигма-точек пространства Y
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }

        this->CorrectionSRUKF( Y_msd );
    }

protected:
    //------------------------------------------------------------------------------------------------------------------
    // Матрицы знаков
    arma::mat J_; ///< Матрица знаков при Pxy
    arma::mat Jcorrect_; ///< Матрица знаков при коррекции

    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Создание матриц знаков
    ///
    void createSignMatrices()
    {
        int JQR_correct_size = ( 2 * SizeX + 1 ) + SizeY;
        Jcorrect_ = arma::mat( JQR_correct_size, JQR_correct_size, arma::fill::eye );
        if( this->weights_covariance_( this->k_sigma_points_- 1 ) < 0.0 ) {
            Jcorrect_( JQR_correct_size-1, JQR_correct_size-1 ) = -1.0;
        }

        J_ = arma::mat( this->k_sigma_points_, this->k_sigma_points_, arma::fill::eye );
        if( this->weights_covariance_( this->k_sigma_points_- 1 ) < 0.0 ) {
            J_( this->k_sigma_points_- 1, this->k_sigma_points_- 1 ) = -1.0;
        }
    }
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief  Шаблонный класс квадратно-корневого расширенного сигма-точечного (ансцентного) фильтра Калмана
///         (блочная фильтрация), КК-РСТФК (Square Root Extended Unscented Kalman Filter Block, SR-EUKFB)
/// \details Источники:
/// \n[1] A New Extension of the Kalman Filter to Nonlinear Systems, Simon J. Julier, Jeffrey K. Uhlmann,
/// The Robotics Research Group, Department of Engineering Science, The University of Oxford, 1997
/// \n[2] The Unscented Kalman Filter for Nonlinear Estimation, Eric A. Wan and Rudolph van der Merwe
/// Oregon Graduate Institute of Science & Technology 20000 NW Walker Rd, Beaverton, Oregon 97006, 2000
/// ericwan[at]ece.ogi.edu, rvdmerwe[at]ece.ogi.edu
/// \n[3] Sigma-Point Kalman Filters for Probabilistic Inference in Dynamic State-Space Models, Rudolph van der Merwe &
/// Eric Wan, OGI School of Science & Engineering Oregon Health & Science University Beaverton, Oregon, 97006, USA, 2003
/// {rvdmerwe,ericwan}[at]ece.ogi.edu
/// \n[4] THE SQUARE-ROOT UNSCENTED KALMAN FILTER FOR STATE AND PARAMETER-ESTIMATION, Rudolph van der Merwe and
/// Eric A. Wan, Oregon Graduate Institute of Science and Technology 20000 NW Walker Road, Beaverton, Oregon 97006, USA
/// rvdmerwe,ericwan[at]ece.ogi.edu
/// \tparam SizeX - размерность пространства состояния X
/// \tparam SizeY - размерность пространства измерений Y
///
template<size_t SizeX, size_t SizeY>
class CKalmanSREUKFB : public CKalmanSREUKF<SizeX, SizeY>
{
public:
    //------------------------------------------------------------------------------------------------------------------
    // Конструкторы:

    ///
    /// \brief Конструктор по умолчанию
    ///
    CKalmanSREUKFB() : CKalmanSREUKF<SizeX, SizeY>()
    {
#ifdef DEBUG_KALMAN
        this->SetFilterName( "SREUKFB" );
#endif
        this->createSignMatricesBlock();
    }
/*
    ///
    /// \brief Конструктор копирования
    /// \param other - экземпляр, с которого делается копия
    ///
    CKalmanSREUKFB( const CKalmanSREUKFB &other ) : CKalmanEKF<SizeX, SizeY>( other ), CKalmanSREUKF<SizeX, SizeY>( other )
    {}

    ///
    /// \brief Перегрузка оператора присвоения
    /// \param other - экземпляр, с которого делается копия
    /// \return *this
    ///
    CKalmanSREUKFB& operator=( const CKalmanSREUKFB &other )
    {
        CKalmanSREUKFB copy( other );
        swap( *this, copy );
        return *this;
    }

    ///
    /// \brief Конструктор перемещения
    /// \param other - экземпляр, с которого делается копия
    ///
    CKalmanSREUKFB( CKalmanSREUKFB &&other ) noexcept
    {
        swap( *this, other );
    }

    ///
    /// \brief Перегрузка оператора перемещения
    /// \param other - экземпляр, с которого делается копия
    /// \return *this
    ///
    CKalmanSREUKFB& operator=( CKalmanSREUKFB &&other ) noexcept
    {
        swap( *this, other );
        return *this;
    }

    virtual ~CKalmanSREUKFB() = default; ///< Дестркутор

    ///
    /// \brief Метод свапа
    /// \param lhs - left hand side instance
    /// \param rhs - right hand side instance
    ///
    friend void swap( CKalmanSREUKFB<SizeX, SizeY> &lhs, CKalmanSREUKFB<SizeX, SizeY> &rhs ) noexcept
    {
        swap( dynamic_cast< CKalmanSREUKF<SizeX, SizeY> &>( lhs ), dynamic_cast< CKalmanSREUKF<SizeX, SizeY> &>( rhs ) ); // Свап предка 1 порядка

        std::swap( lhs.JcorrectBlock_, rhs.JcorrectBlock_ );
    }
*/
    //------------------------------------------------------------------------------------------------------------------
    // Методы прогноза и коррекции:

    // Prediction как в НЕблочном фильтре

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

        assert( this->prediction_isDone ); // Перед фильтрацией обязательно должен быть выполнен прогноз, иначе не имеет смысла
        this->prediction_isDone = false; // Сразу же снять признак

        // 1a. Создание сигма-точек пространства X
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
        // 1b. Вычисление сигма-точек пространства Y
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->y_pred_sigma_points_.col(i) = this->observationModel_( this->x_pred_sigma_points_.col(i) );
            if( this->checkBordersMeasurement_ != nullptr ) {
                this->y_pred_sigma_points_.col(i) = this->checkBordersMeasurement_( this->y_pred_sigma_points_.col(i) );
            }
        }

        // 1. Составление блочной матрицы:
        // B = [ R    dYcal * WcovSqrtDiag ]
        //     [ 0    dXcal * WcovSqrtDiag ]
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dYcal_.col(i) = ( this->y_pred_sigma_points_.col(i) - this->Y_pred_ );
            if( this->checkDeltaMeasurement_ != nullptr ) {
                this->dYcal_.col(i) = this->checkDeltaMeasurement_( this->dYcal_.col(i) );
            }
            this->dYcal_.col(i) *= std::sqrt( std::abs( this->weights_covariance_(i) ) );
        }
        for( int i = 0; i < this->k_sigma_points_; i++ ) {
            this->dXcal_.col(i) = ( this->x_pred_sigma_points_.col(i) - this->X_pred_ );
            if( this->checkDeltaState_ != nullptr ) {
                this->dXcal_.col(i) = this->checkDeltaState_( this->dXcal_.col(i) );
            }
            this->dXcal_.col(i) *= std::sqrt( std::abs( this->weights_covariance_(i) ) );
        }

        size_t rowsB = ( SizeY + SizeX );
        size_t colsB = ( SizeY + ( 2 * SizeX + 1 ) );
        arma::mat B = arma::mat( rowsB, colsB, arma::fill::zeros );
        for( size_t i = 0; i < rowsB; i++ ) {
            for( size_t j = 0; j < colsB; j++ ) {
                if( ( i < SizeY ) && ( j < SizeY ) ) { // Блок R
                    B( i, j ) = this->R_( i, j );
                } else if( ( i < SizeY ) && ( j >= SizeY ) ) { // Блок dYcal
                    B( i, j ) = this->dYcal_( i, j - SizeY );
                } else if( ( i >= SizeY ) && ( j < SizeY ) ) { // Нулевой блок
                    B( i, j ) = 0.0;
                } else { // Блок dXcal if( ( i >= SizeY ) && ( j >= SizeY ) )
                    B( i, j ) = this->dXcal_( i - SizeY, j - SizeY );
                }
            }
        }
        B = arma::trans( B ); // Транспонировать, т.к. в JQR разложение так надо
#ifdef DEBUG_KALMAN
        B.print( this->filterName_ + " Correction, B' ( block of R_, dYcal_, dXcal_ ):" );
#endif
        // 2. Выполнить JQR разложение блочной матрицы B и считать результат:
        // [ S     0 ]
        // [ Pxy*  P ]
        arma::mat Q_qr_corr, R_qr_corr, J_qr_corr;
        int res_JQR_corr;
        if( this->negativeZeroCovWeight_ ) {
            res_JQR_corr = SPML::QR::J_orthogonal( Q_qr_corr, R_qr_corr, J_qr_corr, B, this->JcorrectBlock_, true );
        } else {
            res_JQR_corr = SPML::QR::ModifiedGramSchmidt( Q_qr_corr, R_qr_corr, B );
        }
        assert( res_JQR_corr == 0 );
        R_qr_corr = arma::trans( R_qr_corr ); // Транспонировать, т.к. в JQR разложение так надо
        for( size_t i = 0; i < rowsB; i++ ) {
            for( size_t j = 0; j < rowsB; j++ ) { //for( int j = 0; j < colsB; j++ ) {
                if( ( i < SizeY ) && ( j < SizeY ) ) { // Матрица S
                   this->S_(i, j) = R_qr_corr(i, j);
                } else if( ( i >= SizeY ) && ( j < SizeY ) ) { // Матрица P_xy*
                    this->P_xy_(i - SizeY, j) = R_qr_corr(i, j);
                } else if( ( i >= SizeY ) && ( ( j >= SizeY ) && ( j < SizeY + SizeX ) ) ) { // Матрица P
                    this->P_(i - SizeY, j - SizeY) = R_qr_corr(i, j);
                }
            }
        }
#ifdef DEBUG_KALMAN
        ( this->S_ ).print( this->filterName_ + " Correction, S_:" );
        ( this->P_xy_ ).print( this->filterName_ + " Correction, P_xy:" );
        ( this->P_ ).print( this->filterName_ + " Correction, P_:" );
#endif
        // 3. Вычисление коэффициента усиления фильтра K
        this->K_ = this->P_xy_ * arma::inv( this->S_ ); // В блочном случае - это ВЕРНО!
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

protected:
    //------------------------------------------------------------------------------------------------------------------
    // Матрицы знаков:
    arma::mat JcorrectBlock_; ///< Матрица знаков для фильтра в блочном виде

    //------------------------------------------------------------------------------------------------------------------
    ///
    /// \brief Создание матриц знаков
    ///
    void createSignMatricesBlock()
    {
        int QRsizeY = ( SizeY + ( 2 * SizeX + 1 ) );
        JcorrectBlock_ = arma::mat( QRsizeY, QRsizeY, arma::fill::eye );
        if( this->weights_covariance_( this->k_sigma_points_- 1 ) < 0.0 ) {
            JcorrectBlock_( QRsizeY - 1, QRsizeY - 1 ) = -1.0;
        }
    }
};

}

#endif // KALMAN_FILTER_EXTENDED_UNSCENTED_SQUARE_ROOT_H
/// \}
