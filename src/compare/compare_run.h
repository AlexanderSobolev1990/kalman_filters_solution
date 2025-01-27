//----------------------------------------------------------------------------------------------------------------------
///
/// \file       compare_run.h
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
#define SREKFB_
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
#include <array>
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
#include <omp.h> // parallellization

// Project includes:
#ifdef EKF_
    #include <kalman_ekf.h>
#endif
#ifdef SREKF_
    #include <kalman_srekf.h>
#endif
#ifdef UKF_
    #include <kalman_ukf.h>
#endif
#ifdef SRUKF_
    #include <kalman_srukf.h>
#endif
#ifdef CKF_
    #include <kalman_ckf.h>
#endif
#ifdef SRCKF_
    #include <kalman_srckf.h>
#endif
#ifdef EUKF_
    #include <kalman_eukf.h>
#endif
#ifdef SREUKF_
    #include <kalman_sreukf.h>
#endif
#ifdef ECKF_
    #include <kalman_eckf.h>
#endif
#ifdef SRECKF_
    #include <kalman_sreckf.h>
#endif

#include <matplotlibcpp.h>
#include <timing.h>
#include <geodesy.h>
#include <radars.h>
#include <convert.h>

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Конвертация float или double в строку с определённой точностью
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
    std::vector<double> LocLegend;      ///< Параметры расположения легенды
    uint32_t MCruns;                    ///< Число запусков реализаций ГСЧ
    uint32_t MCseed;                    ///< Зерно ГСЧ

    uint32_t Set;                       ///< Способ установки сигма-точек (0 - Julier, 1 - Merwe)

    std::vector<double> w0;     ///< w0
    std::vector<double> alpha;  ///< alpha
    std::vector<double> beta;   ///< beta
    std::vector<double> kappa;  ///< kappa

    std::vector<double> w0_sr;     ///< w0 sr
    std::vector<double> alpha_sr;  ///< alpha sr
    std::vector<double> beta_sr;   ///< beta sr
    std::vector<double> kappa_sr;  ///< kappa sr

    std::vector<double> x_start; ///< Начальное приближение по Х
    std::vector<double> q_koef_ekf; ///< Коэффициент на который умножается матрица Q для EKF-фильтров
    std::vector<double> q_koef_ukf; ///< Коэффициент на который умножается матрица Q для UKF-фильтров и гибридов
    int eksperiment; ///< Номер эксперимента (модели)
    
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
        LocLegend.clear();
        MCruns = 1;
        MCseed = 1;

        Set = 1;

        w0.clear();
        alpha.clear();
        beta.clear();
        kappa.clear();

        w0_sr.clear();
        alpha_sr.clear();
        beta_sr.clear();
        kappa_sr.clear();

        x_start.clear();
        
        q_koef_ekf.clear();
        q_koef_ukf.clear();

        eksperiment = 1;
    }
};

//----------------------------------------------------------------------------------------------------------------------
namespace IndX /// Индексы соответствующих элементов в векторе состояния X фильтра Калмана
{
const int X = 0; ///< X декаратова координата
const int Y = 1; ///< Y декаратова координата
const int V = 2; ///< Полная скорость
const int K = 3; ///< Курс относительно севера
const int Ka = 4; ///< Производная курса
}

//----------------------------------------------------------------------------------------------------------------------
namespace IndY /// Индексы соответствующих элементов в векторе измерения Y фильтра Калмана
{
const int R = 0; ///< Дальность
const int Az = 1; ///< Азимут относительно севера
const int Vf = 2; ///< Радиальная скорость
}

//----------------------------------------------------------------------------------------------------------------------
namespace Consts /// Константы для алгоритма траекторной обработки
{
//
// Единицы измерения - ИЗМЕНЕНИЕ ПРИВЕДЕТ К ОШИБКАМ ПЕРЕСЧЕТА КООРДИНАТ!
const SPML::Units::TRangeUnit UnitRange = SPML::Units::TRangeUnit::RU_Kilometer; ///< Единица измерения дальности
const SPML::Units::TAngleUnit UnitAngle = SPML::Units::TAngleUnit::AU_Degree; ///< Единица измерения угла
//
// Эллипсоиды для решения геодезических задач
const SPML::Geodesy::CEllipsoid Ellipsoid = SPML::Geodesy::Ellipsoids::Sphere6371(); ///< Эллипсоид для геодезических задач
const double Re = Ellipsoid.A(); ///< Радиус Земли для расчетов на сфере, [м]
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Класс тестирования фильтров Калмана
///
class CKalmanFiltersCompare
{
public:
    double RMS_X_X = 0.005; // км
    double RMS_X_Y = 0.005; // км
    double RMS_X_V = 0.1; // м/c
    double RMS_X_K = 0.3; // град
    double RMS_X_Ka = 0.01;//0.01; // град/с

    double resElementR = 1.0; // км
    double resElementAz = 0.57;//0.1; // град
    double resElementVr = 1.0; // м/с

    double coefK = 2.0;
    double weight_dB = 5;//11.0;//  30; //20; //11.0; //11;//9;//
    double weight_times = std::pow( 10.0, ( weight_dB * 0.1 ) ); // Вес в разах
    double RMS_Y_R = coefK * resElementR / std::sqrt( 12.0 * weight_times );
    double RMS_Y_Az = coefK * resElementAz / std::sqrt( 12.0 * weight_times );
    double RMS_Y_Vr = coefK * resElementVr / std::sqrt( 12.0 * weight_times );

    // Диагональ матрицы шумов состояния:
    double k = 1.0;
    arma::vec Q_CV = {
        RMS_X_X * RMS_X_X * k,
        RMS_X_Y * RMS_X_Y * k,
        RMS_X_V * RMS_X_V * k,
        RMS_X_K * RMS_X_K * k,
        RMS_X_Ka * RMS_X_Ka * k
    };

    const double koefR = 1.0;
    // const double koefR = 10.0;
    // Диагональ матрицы шумов измерений:
    arma::vec R = {
        RMS_Y_R * RMS_Y_R * koefR,
        RMS_Y_Az * RMS_Y_Az * koefR,
        RMS_Y_Vr * RMS_Y_Vr * koefR
    };

    double Pg = 0.0;//0.25;// // Коэффициент glint noise - вероятность появления "пиков"
    double glintNoiseCoef = 10.0;//4.0;//

    ///
    /// \brief Основной метод - запуск тестов
    /// \param settings - настройки
    ///
    void RunMain( const CSettings &settings );

    ///
    /// \brief Запуск теста весов SRUKF в сравнении и EKF
    /// \param settings - настройки
    ///
    void Run_RMSE_var_params( const CSettings &settings );

    bool ylim_yes = false; ///< Лимит на ось у

    bool relatEKF = false; ///< строить RMSE относительно EKF
//    bool relatEKF = true; ///< строить RMSE относительно EKF

//    bool recalcQ = false;
    bool recalcQ = true;

    const SPML::Radar::CRadarPosition RadarSector = SPML::Radar::CRadarPosition(
        "radar", 1, 0, 0, 0, 0, 0, Consts::Ellipsoid, Consts::UnitRange, Consts::UnitAngle );

    ///
    /// \brief Проверка границ вектора состояния Х
    ///
    std::function<arma::vec( const arma::vec &X )> checkBordersState =
        [&]( const arma::vec &X )->arma::vec {
            arma::vec X_ = X;
            X_( IndX::K ) = SPML::Convert::AngleTo360( X_( IndX::K ), SPML::Units::AU_Degree ); // Курс в 0-360
            return X_;
        };

    ///
    /// \brief Проверка границ вектора измерений Y
    ///
    std::function<arma::vec( const arma::vec &Y )> checkBordersMeasurement =
        [&]( const arma::vec &Y )->arma::vec {
            arma::vec Y_ = Y;
            Y_( IndY::Az ) = SPML::Convert::AngleTo360( Y_( IndY::Az ), SPML::Units::AU_Degree ); // Азимут в 0-360
            return Y_;
        };

    ///
    /// \brief Проверка разности векторов состояния Х
    ///
    std::function<arma::vec( const arma::vec &DeltaX )> checkDeltaState =
        [&]( const arma::vec &DeltaX )->arma::vec {
            arma::vec DeltaX_ = DeltaX;
            DeltaX_[ IndX::K ] = SPML::Convert::CheckDeltaAngle( DeltaX_[ IndX::K ], SPML::Units::AU_Degree ); // Курс
            return DeltaX_;
        };

    ///
    /// \brief Проверка разности векторов измерений Y
    ///
    std::function<arma::vec( const arma::vec &DeltaY )> checkDeltaMeasurement =
        [&]( const arma::vec &DeltaY )->arma::vec {
            arma::vec DeltaY_= DeltaY;
            DeltaY_[ IndY::Az ] = SPML::Convert::CheckDeltaAngle( DeltaY_[ IndY::Az ], SPML::Units::AU_Degree );
            return DeltaY_;
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
            double arg1 = ( SPML::Convert::AngleTo360( k + ( ka * dt ), SPML::Units::AU_Degree ) ) * SPML::Convert::DgToRd;
            double arg2 = k * SPML::Convert::DgToRd;
            double KmToMeters = 1000.0;

            double dX = -( v / ( ka * SPML::Convert::DgToRd ) ) * ( std::cos( arg1 ) - std::cos( arg2 ) );
            double dY =  ( v / ( ka * SPML::Convert::DgToRd ) ) * ( std::sin( arg1 ) - std::sin( arg2 ) );

            X_new(0) = X(0) + ( dX / KmToMeters );
            X_new(1) = X(1) + ( dY / KmToMeters );
            X_new(2) = v;
            X_new(3) = SPML::Convert::AngleTo360( ( k + ( ka * dt ) ), Consts::UnitAngle );
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
            Y_(1) = SPML::Convert::AngleTo360( ( std::atan2( X(0), X(1) ) * SPML::Convert::RdToDg ), SPML::Units::AU_Degree ); // atan( x / y ) поскольку азимут отсчитывается от вертикальной оси
            double k = SPML::Convert::AngleTo360( X(3), SPML::Units::AU_Degree ); // K
            double az = SPML::Convert::AngleTo360( Y_(1), SPML::Units::AU_Degree ); // Az

            // 2
            Y_(2) = X(2) * std::cos( SPML::Convert::CheckDeltaAngle( ( k - az ), SPML::Units::AU_Degree ) * SPML::Convert::DgToRd ); // Vr            
            return Y_;
        };

    ///
    /// \brief Матрица перехода состояния Х аналитически (analitically)
    ///
    std::function<arma::mat( const arma::vec &X, double dt )> stateTransitionJacobianF =
        [&]( const arma::vec &X, double dt )->arma::mat {
            // Аналитическое дифференцирование
            double v = X(2);
            double k = X(3);
            double ka = X(4);            
            double arg1 = ( SPML::Convert::AngleTo360( k + ( ka * dt ), SPML::Units::AU_Degree ) ) * SPML::Convert::DgToRd;
            double arg2 = k * SPML::Convert::DgToRd;
            double KmToMeters = 1000.0;

            double f02 = -( std::cos( arg1 ) - std::cos( arg2 ) ) / ( KmToMeters * ka * SPML::Convert::DgToRd );
            double f03 = -( v * SPML::Convert::DgToRd * ( std::sin( arg2 ) - std::sin( arg1 ) ) ) / ( KmToMeters * ka * SPML::Convert::DgToRd );
            double f04 = ( ( dt * v * std::sin( arg1 ) ) / ( KmToMeters * ka ) ) +
                ( v * ( std::cos( arg1 ) - std::cos( arg2 ) ) / ( KmToMeters * ka * ka * SPML::Convert::DgToRd ) );

            double f12 = ( std::sin( arg1 ) - std::sin( arg2 ) ) / ( KmToMeters * ka * SPML::Convert::DgToRd );
            double f13 = ( SPML::Convert::DgToRd * v * ( std::cos( arg1 ) - std::cos( arg2 ) ) ) / ( KmToMeters * ka * SPML::Convert::DgToRd );
            double f14 = ( ( dt * v * std::cos( arg1 ) ) / ( KmToMeters * ka ) ) -
                ( v * ( std::sin( arg1 ) - std::sin( arg2 ) ) / ( KmToMeters * ka * ka * SPML::Convert::DgToRd ) );

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
//                0.001, // X
//                0.001, // Y
//                0.001, // V, [м/с]
//                0.001, // K, [градус]
//                0.001 // Ka, [градус/с]
                0.001, // X
                0.001, // Y
                1.0, // V, [м/с]
                1.0, // K, [градус]
                0.1 // Ka, [градус/с]
            };
            arma::vec Xprev = arma::vec( SizeX, arma::fill::zeros );
            arma::vec Xnext = arma::vec( SizeX, arma::fill::zeros );
            arma::vec Yprev = arma::vec( SizeY, arma::fill::zeros );
            arma::vec Ynext = arma::vec( SizeY, arma::fill::zeros );
            arma::mat H( SizeY, SizeX, arma::fill::zeros );
            // for( size_t j = 0; j < SizeX; j++ ) {
            //     Xprev = X;
            //     Xnext = X;
            //     Xprev(j) -= ( delta[j] * 0.5 );
            //     Xnext(j) += ( delta[j] * 0.5 );
            //     // !!!
            //     Xprev = checkBordersState( Xprev );
            //     Xnext = checkBordersState( Xnext );
            //     // !!!
            //     Yprev = observationModel( Xprev );
            //     Ynext = observationModel( Xnext );
            //     for( size_t i = 0; i < SizeY; i++ ) {
            //         H(i,j) = ( Ynext(i) - Yprev(i) ) / delta[j];
            //     }
            // }
            for( size_t j = 0; j < SizeX; j++ ) {
                Xprev = X;
                Xnext = X;
//                Xprev(j) -= ( delta[j] * 0.5 );
//                Xnext(j) += ( delta[j] * 0.5 );
                if( j == 3 ) { // K
                    Xprev(j) -= ( delta[j] * 0.5 );
                    Xprev(j) = SPML::Convert::AngleTo360( Xprev(j), SPML::Units::AU_Degree );
                    Xnext(j) += ( delta[j] * 0.5 );
                    Xnext(j) = SPML::Convert::AngleTo360( Xnext(j), SPML::Units::AU_Degree );
                } else {
                   Xprev(j) -= ( delta[j] * 0.5 );
                   Xnext(j) += ( delta[j] * 0.5 );
                }

                Yprev = observationModel( Xprev );
                Ynext = observationModel( Xnext );

                for( size_t i = 0; i < SizeY; i++ ) {
                    if( i == 1 ) { // Az
                        double tmp = Ynext(i) - Yprev(i);
                        tmp = SPML::Convert::CheckDeltaAngle( tmp, SPML::Units::AU_Degree );
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
            arma::vec Xsum( ( SizeX + 1 ), arma::fill::zeros );
            arma::vec Xresult( SizeX, arma::fill::zeros );

            for( int i = 0; i < k; i++ ) {
                Xsum( 0 ) += weights( i ) * ( sigmaPoints.col( i ) )[0]; // X
                Xsum( 1 ) += weights( i ) * ( sigmaPoints.col( i ) )[1]; // Y
                Xsum( 2 ) += weights( i ) * ( sigmaPoints.col( i ) )[2]; // V
                Xsum( 3 ) += weights( i ) * std::cos( ( ( sigmaPoints.col( i ) )[3] ) * SPML::Convert::DgToRd ); // cosK
                Xsum( 4 ) += weights( i ) * std::sin( ( ( sigmaPoints.col( i ) )[3] ) * SPML::Convert::DgToRd ); // sinK
                Xsum( 5 ) += weights( i ) * ( sigmaPoints.col( i ) )[4]; // Ka
            }
            Xresult( 0 ) = Xsum( 0 ); // X
            Xresult( 1 ) = Xsum( 1 ); // Y
            Xresult( 2 ) = Xsum( 2 ); // V
            Xresult( 3 ) = SPML::Convert::AngleTo360( ( std::atan2( Xsum( 4 ), Xsum( 3 ) ) * SPML::Convert::RdToDg ), SPML::Units::AU_Degree ); // K
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
                Ysum( 1 ) += weights( i ) * std::cos( ( ( sigmaPoints.col( i ) )[1] ) * SPML::Convert::DgToRd ); // cosAz
                Ysum( 2 ) += weights( i ) * std::sin( ( ( sigmaPoints.col( i ) )[1] ) * SPML::Convert::DgToRd ); // sinAz
                Ysum( 3 ) += weights( i ) * ( sigmaPoints.col( i ) )[2]; // Vr
            }
            Yresult( 0 ) = Ysum( 0 ); // R
            Yresult( 1 ) = SPML::Convert::AngleTo360( ( std::atan2( Ysum( 2 ), Ysum( 1 ) ) * SPML::Convert::RdToDg ), SPML::Units::AU_Degree ); // Az
            Yresult( 2 ) = Ysum( 3 ); // Vr
            return Yresult;
        };

    void print_percent( int cycle, int cycle_max, int &prev_percent )
    {
        int cycle_percent = static_cast<int>( std::round( ( static_cast<double>( cycle ) * 100.0 ) / static_cast<double>( cycle_max ) ) );
        if( ( cycle_percent % 10 == 0 ) && ( cycle_percent != prev_percent ) ) {
            prev_percent = cycle_percent;
            std::cout << cycle_percent << "%" << std::endl;
        }
    }
};

#endif // COMPARE_H
/// \}
