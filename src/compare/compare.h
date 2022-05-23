//----------------------------------------------------------------------------------------------------------------------
///
/// \file       compare.h
/// \brief      Сравнение фильтров Калмана
/// \date       18.03.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef COMPARE_H
#define COMPARE_H

// Включение фильтров:
#define EKF_
#define SREKF_
#define UKF_
#define SRUKF_
#define CKF_
#define SRCKF_
#define ECKF_
#define SRECKF_
#define EUKF_
#define SREUKF_

// System includes:
#include <iostream>
#include <vector>
#include <random>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>

// Project includes:

#ifdef EKF_
    #include <kalman_filter_extended.h>
#endif
#ifdef SREKF_
    #include <kalman_filter_extended_square_root.h>
#endif
#ifdef EUKF_
    #include <kalman_filter_extended_unscented.h>
#endif
#ifdef SREUKF_
    #include <kalman_filter_extended_unscented_square_root.h>
#endif
#ifdef ECKF_
    #include <kalman_filter_extended_cubature.h>
#endif
#ifdef SRECKF_
    #include <kalman_filter_extended_cubature_square_root.h>
#endif
#ifdef UKF_
    #include <kalman_filter_unscented.h>
#endif
#ifdef SRUKF_
    #include <kalman_filter_unscented_square_root.h>
#endif
#ifdef CKF_
    #include <kalman_filter_cubature.h>
#endif
#ifdef SRCKF_
    #include <kalman_filter_cubature_square_root.h>
#endif

#include <matplotlibcpp.h>
#include <timing.h>

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Конвертация float или double в строку с определённоё точностью
///
template <typename T>
std::string to_string_with_precision( const T a_value, const int n = 6 )
{
    static_assert( std::is_same<T, float>::value || std::is_same<T, double>::value, "wrong template class!" );
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

//----------------------------------------------------------------------------------------------------------------------
const size_t SizeX = 5; ///< Размерность вектора состояния: { x, y, V, K, Ka }
const size_t SizeY = 3; ///< Размерность вектора измерений: { R, Az, Vr }

//----------------------------------------------------------------------------------------------------------------------
const double DgToRd = M_PI / 180.0; ///< Градусы в радианы
const double RdToDg = 180.0 / M_PI; ///< Радианы в градусы

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Настройки
///
struct CSettings
{
    bool Debug;                         ///< Признак вывода номера такта и времени в консоль
    uint32_t Seed;                      ///< Зерно ГСЧ (0- рандомное, >0 - какое задано)
    std::vector<std::string> Filters;   ///< Названия фильтров
    std::vector<double> Probabilities;  ///< Вероятности целевой отметки
    double DeltaT;                      ///< Шаг по времени, [с]
    double SimulationTime;              ///< Полное время имитации, [c]
    uint32_t Graphs_0_RMSE_1;           ///< 0 - строить графики состояния, 1 - строить график RMSE
    std::vector<double> Size;           ///< Размеры - ширина/высота
    bool ShowGraphs;                    ///< Показать графики (true - показать, false - сразу в файлыи)
    bool GraphSeparated;                ///< Раздельные графики
    std::string Format;                 ///< Формат файла (png, eps)
    std::vector<double> MatPlotParams;  ///< Параметры полей matplotlib
    uint32_t MCruns;                    ///< Число запусков реализаций ГСЧ

    uint32_t Set;                       ///< Способ установки сигма-точек (0 - Julier, 1 - Merwe)

    std::vector<double> w0;     ///< w0
    std::vector<double> alpha;  ///< alpha
    std::vector<double> beta;   ///< beta
    std::vector<double> kappa;  ///< kappa

    std::vector<double> w0_sr;     ///< w0 sr
    std::vector<double> alpha_sr;  ///< alpha sr
    std::vector<double> beta_sr;   ///< beta sr
    std::vector<double> kappa_sr;  ///< kappa sr

    ///
    /// \brief Конструктор по умолчанию
    ///
    CSettings()
    {
        Debug = false;
        Seed = 1;
        Filters.clear();
        Probabilities.clear();
        DeltaT = 1.0;
        SimulationTime = 100.0 * 60.0;
        Graphs_0_RMSE_1 = 0;
        Size.clear();
        ShowGraphs = false;
        GraphSeparated = false;
        Format = "png";
        MatPlotParams.clear();
        MCruns = 1;

        Set = 1;

        w0.clear();
        alpha.clear();
        beta.clear();
        kappa.clear();

        w0_sr.clear();
        alpha_sr.clear();
        beta_sr.clear();
        kappa_sr.clear();
    }
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Класс тестирования фильтров Калмана
///
class CKalmanFiltersCompare
{
public:
    ///
    /// \brief Основной метод - запуск тестов
    /// \param settings - настройки
    ///
    void RunMain( const CSettings &settings );

    ///
    /// \brief Запуск теста весов SRUKF в сравнении и EKF
    /// \param settings - настройки
    ///
    void Run_RMSE_EKF_SRUKF_SREUKF( const CSettings &settings );

    void print_percent( int cycle, int cycle_max, int &prev_percent );
    
//private:

//    bool ylim_yes = false; ///< Лимит на ось у
    bool ylim_yes = true; ///< Лимит на ось у

    ///
    /// \brief Приведение угла в пределы 0-360 градусов
    /// \param angle - значение приводимого угла, [градусы]
    /// \return Приведенный в пределы 0-360 градусов угол
    ///
    double angleTo360( double angle )
    {
        double _angle = angle;
        double n = std::floor( _angle / 360.0 );
        if( ( _angle >= 360.0 ) || ( _angle < 0.0 ) ) {
            _angle -= ( 360.0 * n );
        }
        return _angle;
    }

    ///
    /// \brief Проверка разности углов
    /// \param deltaAngle - проверяемая разность углов, [градусы]
    /// \return Проверенная разность углов, [градусы]
    ///
    double checkDeltaAngle( double deltaAngle )
    {
        double _deltaAngle = std::fmod( std::abs( deltaAngle ) + 180.0, 360.0 ) - 180.0;
        if( deltaAngle < 0.0 ) {
            _deltaAngle *= ( -1.0 );
        }
        return _deltaAngle;
    }

    ///
    /// \brief Проверка границ вектора состояния Х
    ///
    std::function<arma::vec( const arma::vec &X )> checkBordersState =
        [&]( const arma::vec &X )->arma::vec {
            arma::vec X_ = X;
            if( X_(2) < 0.0 ) {
                X_(2) = std::abs( X_(2) );
            }
            X_(3) = angleTo360( X_(3) ); // Курс в 0-360
            return X_;
        };

    ///
    /// \brief Проверка границ вектора измерений Y
    ///
    std::function<arma::vec( const arma::vec &Y )> checkBordersMeasurement =
        [&]( const arma::vec &Y )->arma::vec {
            arma::vec Y_ = Y;
            Y_(1) = angleTo360( Y_(1) ); // Курс в 0-360
            return Y_;
        };

    ///
    /// \brief Проверка разности векторов состояния Х
    ///
    std::function<arma::vec( const arma::vec &DeltaX )> checkDeltaState =
        [&]( const arma::vec &DeltaX )->arma::vec {
            arma::vec DeltaX_ = DeltaX;
            DeltaX_[3] = checkDeltaAngle( DeltaX_[3] );
            return DeltaX_;
        };

    ///
    /// \brief Проверка разности векторов измерений Y
    ///
    std::function<arma::vec( const arma::vec &DeltaY )> checkDeltaMeasurement =
        [&]( const arma::vec &DeltaY )->arma::vec {
            arma::vec DeltaY_= DeltaY;
            DeltaY_[1] = checkDeltaAngle( DeltaY_[1] );
            return DeltaY;
        };

    ///
    /// \brief Функция прогноза состояния X
    ///
    std::function<arma::vec( const arma::vec &X, double dt )> stateTransitionModel =
        [&]( const arma::vec &X, double dt )->arma::mat {
            arma::vec X_new( SizeX );

            double v = X(2);
            double k = X(3);
            double ka = X(4);
            double arg1 = ( angleTo360( k + ( ka * dt ) ) ) * DgToRd;
//            double arg1 = ( k + ( ka * dt ) ) * DgToRd;
            double arg2 = k * DgToRd;
            double KmToMeters = 1000.0;

            double dX = -( v / ka ) * ( std::cos( arg1 ) - std::cos( arg2 ) );
            double dY =  ( v / ka ) * ( std::sin( arg1 ) - std::sin( arg2 ) );

            X_new(0) = X(0) + ( dX / KmToMeters );
            X_new(1) = X(1) + ( dY / KmToMeters );
            X_new(2) = v;
            X_new(3) = angleTo360( k + ( ka * dt ) );
//            X_new(3) = k + ( ka * dt );
            X_new(4) = ka;
            return X_new;
        };

    ///
    /// \brief Функция перевода состояния X в измерение Y
    ///
    std::function<arma::vec( const arma::vec &X )> observationModel =
        [&]( const arma::vec &X )->arma::mat {
            arma::vec Y_( SizeY );

            // X(0) - X
            // X(1) - Y
            // X(2) - V
            // X(3) - K
            // X(4) - Ka

            Y_(0) = std::sqrt( ( X(0) * X(0) ) + ( X(1) * X(1) ) ); // R
            // 1
            Y_(1) = std::atan2( X(0), X(1) ) * RdToDg; // atan( x / y ) поскольку азимут отсчитывается от вертикальной оси
//            Y_(1) = AngleTo360( Y_(1) ); // Az

            double k = angleTo360( X(3) ); // K
            double az = angleTo360( Y_(1) ); // Az

            Y_(2) = X(2) * std::cos( checkDeltaAngle( k - az ) * DgToRd ); // Vr

//            Y_(2) = X(2) * std::cos( CheckDeltaAngle( X(3) - Y_(1) ) * DgToRd ); // Vr - NEW

//            Y_(2) = X(2) * std::cos( ( X(3) - Y_(1) ) * DgToRd ); // Vr - ORIGINAL

            // 2
//            double tmp = std::atan2( X(0), X(1) );// atan( x / y )
////            Y_(1) = tmp * RdToDg;
//            Y_(1) = AngleTo360( tmp * RdToDg );
//            double cos_K_Az = std::cos( X(3) * DgToRd ) * std::cos( tmp ) + std::sin( X(3) * DgToRd ) * std::sin( tmp );
//            Y_(2) = X(2) * cos_K_Az; // Vr
            return Y_;
        };

    ///
    /// \brief Матрица перехода состояния Х
    ///
    std::function<arma::mat( const arma::vec &X, double dt )> stateTransitionJacobianF =
        [&]( const arma::vec &X, double dt )->arma::mat {
            // Аналитическое дифференцирование
            double v = X(2);
            double k = X(3);
            double ka = X(4);
            double arg1 = ( angleTo360( k + ( ka * dt ) ) ) * DgToRd;
//            double arg1 = ( k + ( ka * dt ) ) * DgToRd;
            double arg2 = k * DgToRd;
            double KmToMeters = 1000.0;

            double f02 = -( std::cos( arg1 ) - std::cos( arg2 ) ) / ( KmToMeters * ka );
            double f03 = -( DgToRd * v * ( std::sin( arg2 ) - std::sin( arg1 ) ) ) / ( KmToMeters * ka );
            double f04 = ( v * ( std::cos( arg1 ) - std::cos( arg2 ) ) / ( KmToMeters * ka * ka ) ) +
                ( ( DgToRd * dt * v * std::sin( arg1 ) ) / ( KmToMeters * ka ) );

            double f12 = ( std::sin( arg1 ) - std::sin( arg2 ) ) / ( KmToMeters * ka );
            double f13 = -( DgToRd * v * ( std::cos( arg2 ) - std::cos( arg1 ) ) ) / ( KmToMeters * ka );
            double f14 = -( v * ( std::sin( arg1 ) - std::sin( arg2 ) ) / ( KmToMeters * ka * ka ) ) +
                ( ( DgToRd * dt * v * std::cos( arg1 ) ) / ( KmToMeters * ka ) );

            arma::mat F = {
                { 1.0, 0.0, f02, f03, f04 },
                { 0.0, 1.0, f12, f13, f14 },
                { 0.0, 0.0, 1.0, 0.0, 0.0 },
                { 0.0, 0.0, 0.0, 1.0, dt  },
                { 0.0, 0.0, 0.0, 0.0, 1.0 }
            };
            return F;

//            // Численное дифференцирование
//            arma::mat F( SizeX, SizeX, arma::fill::zeros );

//            const double delta[SizeX] = {
//                10.0, // X км
//                10.0, // Y км
//                2.7778, // V, [м/с]
//                5.0,//5.0, // K, [градус]
//                0.1 // Ka, [градус/с]
//            };

//            arma::vec Xprev = arma::vec( SizeX, arma::fill::zeros );
//            arma::vec Xnext = arma::vec( SizeX, arma::fill::zeros );

//            for( int j = 0; j < SizeX; j++ ) {
//                Xprev = X;
//                Xnext = X;
//                Xprev(j) -= ( delta[j] * 0.5 );
//                Xnext(j) += ( delta[j] * 0.5 );
//                Xprev = stateTransitionModel( Xprev, dt );
//                Xnext = stateTransitionModel( Xnext, dt );
//                for( int i = 0; i < SizeX; i++ ) {
//                    F(i,j) = ( Xnext(i) - Xprev(i) ) / delta[j];
//                }
//            }
//            return F;
        };

    ///
    /// \brief Матрица перехода состояния X в измерение Y
    ///
    std::function<arma::mat( const arma::vec &X )> observationJacobianH =
        [&]( const arma::vec &X )->arma::mat {
            const double delta[SizeX] = {
                0.001, // X
                0.001, // Y
                0.001, // V, [м/с]
                0.001, // K, [градус]
                0.001 // Ka, [градус/с]
            };
            arma::vec Xprev = arma::vec( SizeX, arma::fill::zeros );
            arma::vec Xnext = arma::vec( SizeX, arma::fill::zeros );
            arma::vec Yprev = arma::vec( SizeY, arma::fill::zeros );
            arma::vec Ynext = arma::vec( SizeY, arma::fill::zeros );
            arma::mat H( SizeY, SizeX, arma::fill::zeros );
//            for( int j = 0; j < SizeX; j++ ) {
//                Xprev = X;
//                Xnext = X;
//                Xprev(j) -= ( delta[j] * 0.5 );
//                Xnext(j) += ( delta[j] * 0.5 );

//                Yprev = observationModel( Xprev );
//                Ynext = observationModel( Xnext );

//                for( int i = 0; i < SizeY; i++ ) {
//                    H(i,j) = ( Ynext(i) - Yprev(i) ) / delta[j];
//                }
//            }

            for( size_t j = 0; j < SizeX; j++ ) {
                Xprev = X;
                Xnext = X;
                Xprev(j) -= ( delta[j] * 0.5 );
                Xnext(j) += ( delta[j] * 0.5 );
//                if( j == 3 ) { // K
//                    Xprev(j) -= ( delta[j] * 0.5 );
//                    Xprev(j) = CheckDeltaAngle( Xprev(j) );
//                    Xnext(j) += ( delta[j] * 0.5 );
//                    Xnext(j) = AngleTo360( Xnext(j) );
//                } else {
//                   Xprev(j) -= ( delta[j] * 0.5 );
//                   Xnext(j) += ( delta[j] * 0.5 );
//                }

                Yprev = observationModel( Xprev );
                Ynext = observationModel( Xnext );

                for( size_t i = 0; i < SizeY; i++ ) {
                    if( i == 1 ) { // Az
                        double tmp = Ynext(i) - Yprev(i);
                        tmp = checkDeltaAngle( tmp );
                        H(i,j) = tmp / delta[j];
                    } else {
                        H(i,j) = ( Ynext(i) - Yprev(i) ) / delta[j];
                    }
                }
            }
            return H;
        };

    ///
    /// \brief Функция вычисления взвешенной суммы векторов состояния Х
    ///
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumStateSigmas =
        [&]( const arma::vec &weights, const arma::mat &sigmaPoints )->arma::vec {
            int k = weights.n_elem; // Число сигма-точек
            // Курс "распадается" на cos и sin составляющую:
            // 0 - X
            // 1 - Y
            // 2 - V
            // 3 - cos(K)
            // 4 - sin(K)
            // 5 - Ka
//            arma::vec sigmaPoint( ( SizeX + 1 ), arma::fill::zeros );
            arma::vec Xsum( ( SizeX + 1 ), arma::fill::zeros );
            arma::vec Xresult( SizeX, arma::fill::zeros );

            for( int i = 0; i < k; i++ ) {
                Xsum( 0 ) += weights( i ) * ( sigmaPoints.col( i ) )[0]; // X
                Xsum( 1 ) += weights( i ) * ( sigmaPoints.col( i ) )[1]; // Y
                Xsum( 2 ) += weights( i ) * ( sigmaPoints.col( i ) )[2]; // V
                Xsum( 3 ) += weights( i ) * std::cos( ( ( sigmaPoints.col( i ) )[3] ) * DgToRd ); // cosK
                Xsum( 4 ) += weights( i ) * std::sin( ( ( sigmaPoints.col( i ) )[3] ) * DgToRd ); // sinK
                Xsum( 5 ) += weights( i ) * ( sigmaPoints.col( i ) )[4]; // Ka
            }
            Xresult( 0 ) = Xsum( 0 ); // X
            Xresult( 1 ) = Xsum( 1 ); // Y
            Xresult( 2 ) = Xsum( 2 ); // V
            Xresult( 3 ) = std::atan2( Xsum( 4 ), Xsum( 3 ) ) * RdToDg ; // K
            Xresult( 4 ) = Xsum( 5 );
            return Xresult;
        };

    ///
    /// \brief Функция вычисления взвешенной суммы векторов измерений Y
    ///
    std::function<arma::vec( const arma::vec &weights, const arma::mat &sigmaPoints )> weightedSumMeasurementSigmas =
        [&]( const arma::vec &weights, const arma::mat &sigmaPoints )->arma::vec {
            int k = weights.n_elem; // Число сигма-точек
            // Азимут "распадается" на cos и sin составляющую:
            // 0 - R
            // 1 - cos(Az)
            // 2 - sin(Az)
            // 3 - Vr
            arma::vec Ysum( ( SizeY + 1 ), arma::fill::zeros );
            arma::vec Yresult( SizeY, arma::fill::zeros );

            for( int i = 0; i < k; i++ ) {
                Ysum( 0 ) += weights( i ) * ( sigmaPoints.col( i ) )[0]; // R
                Ysum( 1 ) += weights( i ) * std::cos( ( ( sigmaPoints.col( i ) )[1] ) * DgToRd ); // cosAz
                Ysum( 2 ) += weights( i ) * std::sin( ( ( sigmaPoints.col( i ) )[1] ) * DgToRd ); // sinAz
                Ysum( 3 ) += weights( i ) * ( sigmaPoints.col( i ) )[2]; // Vr
            }
            Yresult( 0 ) = Ysum( 0 ); // R
            Yresult( 1 ) = std::atan2( Ysum( 2 ), Ysum( 1 ) ) * RdToDg ; // Az
            Yresult( 2 ) = Ysum( 3 ); // Vr
            return Yresult;
        };
};

#endif // COMPARE_H
/// \}
