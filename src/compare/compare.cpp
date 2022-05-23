//----------------------------------------------------------------------------------------------------------------------
///
/// \file       compare.cpp
/// \brief      Сравнение фильтров Калмана
/// \date       06.04.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#include <compare.h>

/////
///// \brief
/////
//template<typename T>
//std::vector<T> arange( T start, T stop, T step = 1 )
//{
//    std::vector<T> values;
//    for( T value = start; value < stop; value += step ) {
//        values.push_back( value );
//    }
//    return values;
//}


void CKalmanFiltersCompare::print_percent( int cycle, int cycle_max, int &prev_percent )
{
    int cycle_percent = static_cast<int>( std::round( ( static_cast<double>( cycle ) * 100.0 ) / static_cast<double>( cycle_max ) ) );
    if( ( cycle_percent % 10 == 0 ) && ( cycle_percent != prev_percent ) ) {
        prev_percent = cycle_percent;
        std::cout << cycle_percent << "%" << std::endl;
    }
}

void CKalmanFiltersCompare::RunMain( const CSettings &settings )
{
    try{
    if( settings.Filters.size() == 0 ) {
        std::cout << "No filters selected! Read --help" << std::endl;
        return;
    }

    Py_Initialize();

    namespace plt = matplotlibcpp;
    std::map<std::string, double> keywords;
//    keywords.insert( std::make_pair( "left", 0.08 ) );
//    keywords.insert( std::make_pair( "bottom", 0.08 ) );
//    keywords.insert( std::make_pair( "right", 0.92 ) );
//    keywords.insert( std::make_pair( "top", 0.92 ) );
//    keywords.insert( std::make_pair( "wspace", 0.4 ) );
//    keywords.insert( std::make_pair( "hspace", 0.4 ) ); // 0.7
    keywords.insert( std::make_pair( "left", settings.MatPlotParams[0] ) );
    keywords.insert( std::make_pair( "bottom", settings.MatPlotParams[1] ) );
    keywords.insert( std::make_pair( "right", settings.MatPlotParams[2] ) );
    keywords.insert( std::make_pair( "top", settings.MatPlotParams[3] ) );
    keywords.insert( std::make_pair( "wspace", settings.MatPlotParams[4] ) );
    keywords.insert( std::make_pair( "hspace", settings.MatPlotParams[5] ) );

    const double in2mm = 25.4;// mm (fixed)
//    const double pt2mm = 0.3528;// mm (fixed)

    const double dpi = 300;// dpi (variable)
//    double koef = 1.7; // 1.5; //
//    const double width = 195 * pt2mm * koef;// 250 ## mm (variable)
//    const double height =  130 * pt2mm * koef;// 166 ## mm (variable)
    double width = settings.Size[0];
    double height = settings.Size[1];

    const double mm2px = dpi / in2mm;//
    size_t pixels_width = std::round( width * mm2px);//
    size_t pixels_height = std::round( height * mm2px);//
    //------------------------------------------------------------------------------------------------------------------
    std::map<std::string, std::string> true_keywords = { { "color", "black" }, { "linestyle", "--" },
        { "linewidth", "1" }, { "label", "Истинное" } };
    std::map<std::string, std::string> marks_keywords = { { "color", "grey" }, { "linestyle", "" }, { "marker", "." },
        { "label", "Измерение" } };

    std::map<std::string, std::string> graphColors = {
#ifdef EKF_
        { "EKF", "blue" },
#endif
#ifdef UKF_
        { "UKF", "darkgreen" },
#endif
#ifdef SRUKF_
        { "SRUKF", "limegreen" },
        { "SRUKFB", "lime" },
#endif
#ifdef CKF_
        { "CKF", "darkred" },
#endif
#ifdef SRCKF_
        { "SRCKF", "indianred" },
        { "SRCKFB", "red" },
#endif
#ifdef ECKF_
        { "ECKF", "purple" },
#endif
#ifdef SRECKF_
        { "SRECKF", "mediumpurple" },
        { "SRECKFB", "magenta" },
#endif
#ifdef EUKF_
        { "EUKF", "chocolate" },
#endif
#ifdef SREUKF_
        { "SREUKF", "darkorange" },
        { "SREUKFB", "orange" },
#endif
#ifdef SREKF_
        { "SREKF", "cornflowerblue" }
#endif
    };
    std::map<std::string, std::map<std::string, std::string>> estimated_keywords = {
#ifdef EKF_
        { "EKF", { { "color", graphColors.at("EKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "EKF" } } },
#endif
#ifdef UKF_
        { "UKF", { { "color", graphColors.at("UKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "UKF" } } },
#endif
#ifdef SRUKF_
        { "SRUKF", { { "color", graphColors.at("SRUKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRUKF" } } },
        { "SRUKFB", { { "color", graphColors.at("SRUKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRUKFB" } } },
#endif
#ifdef CKF_
        { "CKF", { { "color", graphColors.at("CKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "CKF" } } },
#endif
#ifdef SRCKF_
        { "SRCKF", { { "color", graphColors.at("SRCKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRCKF" } } },
        { "SRCKFB", { { "color", graphColors.at("SRCKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRCKFB" } } },
#endif
#ifdef ECKF_
        { "ECKF", { { "color", graphColors.at("ECKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "ECKF" } } },
#endif
#ifdef SRECKF_
        { "SRECKF", { { "color", graphColors.at("SRECKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRECKF" } } },
        { "SRECKFB", { { "color", graphColors.at("SRECKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRECKFB" } } },
#endif
#ifdef EUKF_
        { "EUKF", { { "color", graphColors.at("EUKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "EUKF" } } },
#endif
#ifdef SREUKF_
        { "SREUKF", { { "color", graphColors.at("SREUKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREUKF" } } },
        { "SREUKFB", { { "color", graphColors.at("SREUKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREUKFB" } } },
#endif
#ifdef SREKF_
        { "SREKF", { { "color", graphColors.at("SREKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREKF" } } }
#endif
    };
    std::map<std::string, std::map<std::string, std::string>> estimated_keywords_time = {
#ifdef EKF_
        { "EKF", { { "color", graphColors.at("EKF") }, { "label", "EKF" } } },
#endif
#ifdef UKF_
        { "UKF", { { "color", graphColors.at("UKF") }, { "label", "UKF" } } },
#endif
#ifdef SRUKF_
        { "SRUKF", { { "color", graphColors.at("SRUKF") }, { "label", "SRUKF" } } },
        { "SRUKFB", { { "color", graphColors.at("SRUKFB") }, { "label", "SRUKFB" } } },
#endif
#ifdef CKF_
        { "CKF", { { "color", graphColors.at("CKF") }, { "label", "CKF" } } },
#endif
#ifdef SRCKF_
        { "SRCKF", { { "color", graphColors.at("SRCKF") }, { "label", "SRCKF" } } },
        { "SRCKFB", { { "color", graphColors.at("SRCKFB") }, { "label", "SRCKFB" } } },
#endif
#ifdef ECKF_
        { "ECKF", { { "color", graphColors.at("ECKF") }, { "label", "ECKF" } } },
#endif
#ifdef SRECKF_
        { "SRECKF", { { "color", graphColors.at("SRECKF") }, { "label", "SRECKF" } } },
        { "SRECKFB", { { "color", graphColors.at("SRECKFB") }, { "label", "SRECKFB" } } },
#endif
#ifdef EUKF_
        { "EUKF", { { "color", graphColors.at("EUKF") }, { "label", "EUKF" } } },
#endif
#ifdef SREUKF_
        { "SREUKF", { { "color", graphColors.at("SREUKF") }, { "label", "SREUKF" } } },
        { "SREUKFB", { { "color", graphColors.at("SREUKFB") }, { "label", "SREUKFB" } } },
#endif
#ifdef SREKF_
        { "SREKF", { { "color", graphColors.at("SREKF") }, { "label", "SREKF" } } },
#endif
    };

    std::string filters_names;
    for( unsigned long long i = 0; i < settings.Filters.size(); i++ ) {
        filters_names += settings.Filters[i];
        if( i < settings.Filters.size() - 1 ) {
            filters_names += "_";
        }
    }
    std::map<std::string, std::string>legend_loc { { "loc", "upper right" } };
    std::map<std::string, std::string>legend_loc_time { { "loc", "upper left" } };
    std::string loc_time = "upper left";

    //------------------------------------------------------------------------------------------------------------------
    // Буферы
    int N = static_cast<int>( settings.SimulationTime / settings.DeltaT ); // Количество моментов времени

    std::vector<double> time(N);
    arma::mat true_X( SizeX, N, arma::fill::zeros ); // Состояние Х (столбец - вектор в i-ый момент времени)
    arma::mat true_Y( SizeY, N, arma::fill::zeros );
    arma::mat measured_Y( SizeY, N, arma::fill::zeros );

    std::vector<double> template_vector_N(N);
    arma::mat template_mat_X( SizeX, N, arma::fill::zeros );
    arma::mat template_mat_Y( SizeY, N, arma::fill::zeros );

    std::map<std::string, arma::mat> template_map_X;
    std::map<std::string, arma::mat> template_map_Y;
    std::map<std::string, std::vector<double>> template_map_N;

    std::map<std::string, SPML::Timing::CTimeKeeper> timerPrediction;
    std::map<std::string, SPML::Timing::CTimeKeeper> timerCorrection;
    std::map<std::string, double> timerPredictionAverage;
    std::map<std::string, double> timerCorrectionAverage;

    for( auto &filter : settings.Filters ) {
        template_map_X.insert( std::make_pair( filter, template_mat_X ) );
        template_map_Y.insert( std::make_pair( filter, template_mat_Y ) );
        template_map_N.insert( std::make_pair( filter, template_vector_N ) );
        timerPrediction.insert( std::make_pair( filter, SPML::Timing::CTimeKeeper() ) );
        timerCorrection.insert( std::make_pair( filter, SPML::Timing::CTimeKeeper() ) );
        timerPredictionAverage.insert( std::make_pair( filter, 0.0 ) );
        timerCorrectionAverage.insert( std::make_pair( filter, 0.0 ) );
    }
    std::map<std::string, arma::mat> estimated_X = template_map_X;
    std::map<std::string, arma::mat> estimated_Y = template_map_Y;
    std::map<std::string, arma::mat> estimated_Pdiag = template_map_X;
    std::map<std::string, arma::mat> estimated_Sdiag = template_map_Y;
    std::map<std::string, arma::mat> delta_Y = template_map_Y;
    std::map<std::string, std::vector<double>> mahalanobis = template_map_N;
    std::map<std::string, std::vector<double>> SDCM = template_map_N;
    std::map<std::string, arma::mat> RMSE_X = template_map_X;

    std::map<std::string, std::vector<double>> timeMeasuredPrediction = template_map_N;
    std::map<std::string, std::vector<double>> timeMeasuredCorrection = template_map_N;
    std::map<std::string, std::vector<double>> timeMeasuredSumm = template_map_N;    
    std::string name_tmp;

    std::mt19937 generator; // Генератор псевдослучайных чисел Mersenne Twister

    double RMS_X_X = 0.001; // км
    double RMS_X_Y = 0.001; // км
    double RMS_X_V = 0.1;  // м/c
    double RMS_X_K = 0.3; // град
    double RMS_X_Ka = 0.04; // град/с

    double resElementR = 0.30; // км
    double resElementAz = 0.05;//0.1; // град
    double resElementVr = 0.1; // м/с

    double coefK = 2.0;
    double weight_dB = 11.0; //11;//9;//
    double weight_times = std::pow( 10.0, ( weight_dB * 0.1 ) ); // Вес в разах
    double RMS_Y_R = coefK * resElementR / std::sqrt( 12.0 * weight_times );
    double RMS_Y_Az = coefK * resElementAz / std::sqrt( 12.0 * weight_times );
    double RMS_Y_Vr = coefK * resElementVr / std::sqrt( 12.0 * weight_times );

    std::normal_distribution<double> noiseY_R( 0.0, RMS_Y_R );
    std::normal_distribution<double> noiseY_Az( 0.0, RMS_Y_Az );
    std::normal_distribution<double> noiseY_Vr( 0.0, RMS_Y_Vr );

    std::uniform_real_distribution<double> random_0_1( 0.0, 1.0 ); // Вещественное случайное число от 0 до 1 с равномерной плотностью вероятности

    // Диагональ матрицы шумов состояния:
    arma::vec Q = {
        RMS_X_X * RMS_X_X,
        RMS_X_Y * RMS_X_Y,
        RMS_X_V * RMS_X_V,
        RMS_X_K * RMS_X_K,
        RMS_X_Ka * RMS_X_Ka
    };
    // Диагональ матрицы шумов измерений:
    arma::vec R = {
        RMS_Y_R * RMS_Y_R,
        RMS_Y_Az * RMS_Y_Az,
        RMS_Y_Vr * RMS_Y_Vr,
    };

#ifdef EKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanEKF<SizeX, SizeY> EKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "EKF" ) ) {
        EKF.SetStateTransitionModel( stateTransitionModel );
        EKF.SetObservationModel( observationModel );
        EKF.SetStateTransitionJacobianF( stateTransitionJacobianF );        
        EKF.SetObservationJacobianH( observationJacobianH );

        EKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        EKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        EKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        EKF.SetCheckDeltaState( checkDeltaState );
        EKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        EKF.SetProcessCovarianceMatrixQdiag( Q );
        EKF.SetObservationCovarianceMatrixRdiag( R );
    }
#endif
#ifdef UKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanUKF<SizeX, SizeY> UKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "UKF" ) ) {

        if( settings.Set == 0 ) { // Julier            
            UKF.SetupDesignParametersMeanSet( settings.w0[0] );
        } else if( settings.Set == 1 ) { // Merwe            
            UKF.SetupDesignParametersScaledSet( settings.alpha[0], settings.beta[0], settings.kappa[0] );
        } else {
            assert( false );
        }

        UKF.SetStateTransitionModel( stateTransitionModel );
        UKF.SetObservationModel( observationModel );

        UKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        UKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        UKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        UKF.SetCheckDeltaState( checkDeltaState );
        UKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        UKF.SetWeightedSumStateSigmas( weightedSumStateSigmas );
        UKF.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

        UKF.SetProcessCovarianceMatrixQdiag( Q );
        UKF.SetObservationCovarianceMatrixRdiag( R );
    }
#endif
#ifdef SRUKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanSRUKF<SizeX, SizeY> SRUKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {

        if( settings.Set == 0 ) { // Julier            
            SRUKF.SetupDesignParametersMeanSet( settings.w0_sr[0] );
        } else if( settings.Set == 1 ) { // Merwe            
            SRUKF.SetupDesignParametersScaledSet( settings.alpha_sr[0], settings.beta_sr[0], settings.kappa_sr[0] );
        } else {
            assert( false );
        }

        SRUKF.SetStateTransitionModel( stateTransitionModel );
        SRUKF.SetObservationModel( observationModel );

        SRUKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        SRUKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        SRUKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        SRUKF.SetCheckDeltaState( checkDeltaState );
        SRUKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        SRUKF.SetWeightedSumStateSigmas( weightedSumStateSigmas );
        SRUKF.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

        SRUKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
        SRUKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    }
    KalmanFilters::CKalmanSRUKFB<SizeX, SizeY> SRUKFB;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKFB" ) ) {

        if( settings.Set == 0 ) { // Julier            
            SRUKFB.SetupDesignParametersMeanSet( settings.w0_sr[0] );
        } else if( settings.Set == 1 ) { // Merwe            
            SRUKFB.SetupDesignParametersScaledSet( settings.alpha_sr[0], settings.beta_sr[0], settings.kappa_sr[0] );
        } else {
            assert( false );
        }

        SRUKFB.SetStateTransitionModel( stateTransitionModel );
        SRUKFB.SetObservationModel( observationModel );

        SRUKFB.SetCheckBordersStateAfterPrediction( checkBordersState );
        SRUKFB.SetCheckBordersStateAfterCorrection( checkBordersState );
        SRUKFB.SetCheckBordersMeasurement( checkBordersMeasurement );
        SRUKFB.SetCheckDeltaState( checkDeltaState );
        SRUKFB.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        SRUKFB.SetWeightedSumStateSigmas( weightedSumStateSigmas );
        SRUKFB.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

        SRUKFB.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
        SRUKFB.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    }
#endif
#ifdef CKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanCKF<SizeX, SizeY> CKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "CKF" ) ) {
        CKF.SetStateTransitionModel( stateTransitionModel );
        CKF.SetObservationModel( observationModel );

        CKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        CKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        CKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        CKF.SetCheckDeltaState( checkDeltaState );
        CKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        CKF.SetWeightedSumStateSigmas( weightedSumStateSigmas );
        CKF.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

        CKF.SetProcessCovarianceMatrixQdiag( Q );
        CKF.SetObservationCovarianceMatrixRdiag( R );
    }
#endif
#ifdef SRCKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanSRCKF<SizeX, SizeY> SRCKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
        SRCKF.SetStateTransitionModel( stateTransitionModel );
        SRCKF.SetObservationModel( observationModel );

        SRCKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        SRCKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        SRCKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        SRCKF.SetCheckDeltaState( checkDeltaState );
        SRCKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        SRCKF.SetWeightedSumStateSigmas( weightedSumStateSigmas );
        SRCKF.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

        SRCKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
        SRCKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    }
    KalmanFilters::CKalmanSRCKFB<SizeX, SizeY> SRCKFB;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKFB" ) ) {
        SRCKFB.SetStateTransitionModel( stateTransitionModel );
        SRCKFB.SetObservationModel( observationModel );

        SRCKFB.SetCheckBordersStateAfterPrediction( checkBordersState );
        SRCKFB.SetCheckBordersStateAfterCorrection( checkBordersState );
        SRCKFB.SetCheckBordersMeasurement( checkBordersMeasurement );
        SRCKFB.SetCheckDeltaState( checkDeltaState );
        SRCKFB.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        SRCKFB.SetWeightedSumStateSigmas( weightedSumStateSigmas );
        SRCKFB.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

        SRCKFB.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
        SRCKFB.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    }
#endif
#ifdef ECKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanECKF<SizeX, SizeY> ECKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "ECKF" ) ) {
        ECKF.SetStateTransitionModel( stateTransitionModel );
        ECKF.SetStateTransitionJacobianF( stateTransitionJacobianF );
        ECKF.SetObservationModel( observationModel );
        ECKF.SetObservationJacobianH( observationJacobianH );

        ECKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        ECKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        ECKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        ECKF.SetCheckDeltaState( checkDeltaState );
        ECKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        ECKF.SetProcessCovarianceMatrixQdiag( Q );
        ECKF.SetObservationCovarianceMatrixRdiag( R );
    }
#endif
#ifdef SRECKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanSRECKF<SizeX, SizeY> SRECKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
        SRECKF.SetStateTransitionModel( stateTransitionModel );
        SRECKF.SetStateTransitionJacobianF( stateTransitionJacobianF );
        SRECKF.SetObservationModel( observationModel );
        SRECKF.SetObservationJacobianH( observationJacobianH );

        SRECKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        SRECKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        SRECKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        SRECKF.SetCheckDeltaState( checkDeltaState );
        SRECKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        SRECKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
        SRECKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    }
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanSRECKFB<SizeX, SizeY> SRECKFB;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKFB" ) ) {
        SRECKFB.SetStateTransitionModel( stateTransitionModel );
        SRECKFB.SetStateTransitionJacobianF( stateTransitionJacobianF );
        SRECKFB.SetObservationModel( observationModel );
        SRECKFB.SetObservationJacobianH( observationJacobianH );

        SRECKFB.SetCheckBordersStateAfterPrediction( checkBordersState );
        SRECKFB.SetCheckBordersStateAfterCorrection( checkBordersState );
        SRECKFB.SetCheckBordersMeasurement( checkBordersMeasurement );
        SRECKFB.SetCheckDeltaState( checkDeltaState );
        SRECKFB.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        SRECKFB.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
        SRECKFB.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    }
#endif
#ifdef EUKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanEUKF<SizeX, SizeY> EUKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "EUKF" ) ) {

        if( settings.Set == 0 ) { // Julier            
            EUKF.SetupDesignParametersMeanSet( settings.w0[0] );
        } else if( settings.Set == 1 ) { // Merwe            
            EUKF.SetupDesignParametersScaledSet( settings.alpha[0], settings.beta[0], settings.kappa[0] );
        } else {
            assert( false );
        }

        EUKF.SetStateTransitionModel( stateTransitionModel );
        EUKF.SetStateTransitionJacobianF( stateTransitionJacobianF );
        EUKF.SetObservationModel( observationModel );
        EUKF.SetObservationJacobianH( observationJacobianH );

        EUKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        EUKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        EUKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        EUKF.SetCheckDeltaState( checkDeltaState );
        EUKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        EUKF.SetProcessCovarianceMatrixQdiag( Q );
        EUKF.SetObservationCovarianceMatrixRdiag( R );
    }
#endif
#ifdef SREUKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanSREUKF<SizeX, SizeY> SREUKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {

        if( settings.Set == 0 ) { // Julier            
            SREUKF.SetupDesignParametersMeanSet( settings.w0_sr[0] );
        } else if( settings.Set == 1 ) { // Merwe            
            SREUKF.SetupDesignParametersScaledSet( settings.alpha_sr[0], settings.beta_sr[0], settings.kappa_sr[0] );
        } else {
            assert( false );
        }

        SREUKF.SetStateTransitionModel( stateTransitionModel );
        SREUKF.SetStateTransitionJacobianF( stateTransitionJacobianF );
        SREUKF.SetObservationModel( observationModel );
        SREUKF.SetObservationJacobianH( observationJacobianH );

        SREUKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        SREUKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        SREUKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        SREUKF.SetCheckDeltaState( checkDeltaState );
        SREUKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        SREUKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
        SREUKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    }
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanSREUKFB<SizeX, SizeY> SREUKFB;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKFB" ) ) {

        if( settings.Set == 0 ) { // Julier            
            SREUKFB.SetupDesignParametersMeanSet( settings.w0_sr[0] );
        } else if( settings.Set == 1 ) { // Merwe            
            SREUKFB.SetupDesignParametersScaledSet( settings.alpha_sr[0], settings.beta_sr[0], settings.kappa_sr[0] );
        } else {
            assert( false );
        }

        SREUKFB.SetStateTransitionModel( stateTransitionModel );
        SREUKFB.SetStateTransitionJacobianF( stateTransitionJacobianF );
        SREUKFB.SetObservationModel( observationModel );
        SREUKFB.SetObservationJacobianH( observationJacobianH );

        SREUKFB.SetCheckBordersStateAfterPrediction( checkBordersState );
        SREUKFB.SetCheckBordersStateAfterCorrection( checkBordersState );
        SREUKFB.SetCheckBordersMeasurement( checkBordersMeasurement );
        SREUKFB.SetCheckDeltaState( checkDeltaState );
        SREUKFB.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        SREUKFB.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
        SREUKFB.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    }
#endif
#ifdef SREKF_
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanSREKF<SizeX, SizeY> SREKF;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREKF" ) ) {
        SREKF.SetStateTransitionModel( stateTransitionModel );
        SREKF.SetStateTransitionJacobianF( stateTransitionJacobianF );
        SREKF.SetObservationModel( observationModel );
        SREKF.SetObservationJacobianH( observationJacobianH );

        SREKF.SetCheckBordersStateAfterPrediction( checkBordersState );
        SREKF.SetCheckBordersStateAfterCorrection( checkBordersState );
        SREKF.SetCheckBordersMeasurement( checkBordersMeasurement );
        SREKF.SetCheckDeltaState( checkDeltaState );
        SREKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

        SREKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
        SREKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    }
#endif

//    //------------------------------------------------------------------------------------------------------------------
//    // Определим сколько итераций по параметрам настройки UKF
//    int iter = 0;
//    int iter_sr = 0;
//    if( settings.Set == 0 ) { // Julier
//        iter = settings.w0.size();
//        iter_sr = settings.w0_sr.size();
//    } else if( settings.Set == 1 ) { // Merwe
//        iter = settings.alpha.size() * settings.beta.size() * settings.kappa.size();
//        iter_sr = settings.alpha_sr.size() * settings.beta_sr.size() * settings.kappa_sr.size();
//    } else {
//        assert( false );
//    }
    //------------------------------------------------------------------------------------------------------------------
    // Начало рабочих циклов
    //------------------------------------------------------------------------------------------------------------------

    int cycle_max = settings.MCruns * N * settings.Filters.size();
    int cycle = 0;
    int prev_percent = -1;

    for( auto &p_ : settings.Probabilities ) {
        std::string name_p = "_P_" + to_string_with_precision( p_, 1 ) + "_dt_" + std::to_string( settings.DeltaT );
        //--------------------------------------------------------------------------------------------------------------
        for( uint32_t seed_ = 1; seed_ <= settings.MCruns; seed_++ ) {
//            if( seed_ != 3 ) { //1 //3
//                continue;
//            }
            generator.seed( seed_ ); // Выставить зерно ГСЧ
//            std::string name_seed_p = "_seed_" + std::to_string( seed_ ) + "_P_" + to_string_with_precision( p_, 1 ) + "_dt_" + std::to_string( settings.DeltaT );
            std::string name_seed_p = "_seed_" + std::to_string( seed_ ) + "_" + name_p;
            //----------------------------------------------------------------------------------------------------------
            // Начальное положение X, Y, Delta, P
            true_X.col(0) = arma::vec{ 100.0, 200.0, 100.0, 45.0, 1.0e-5 };
            true_Y.col(0) = observationModel( true_X.col(0) );
            measured_Y.col(0) = true_Y.col(0) + arma::vec{ noiseY_R( generator ), noiseY_Az( generator ), noiseY_Vr( generator ) };
            ///
            double v0 = std::abs( ( measured_Y.col(0) )(2) ); // Vr
            double NV = 2.0;// Более или равен 1 (эмпирически)
            double Vmax = 300.0; // м/с
            double startV = ( Vmax / NV ) + ( ( ( NV - 1.0 ) / ( NV * Vmax ) ) * v0 * v0 );
            int sign = 1;
            if( ( measured_Y.col(0) )(2) < 0 ) { // Vr < 0
                sign = -1;
            }
            double startK = ( measured_Y.col(0) )(1) * sign; // Az
            ///
            arma::vec startX {
                ( measured_Y.col(0) )(0) * std::sin( ( measured_Y.col(0) )(1) * DgToRd ),
                ( measured_Y.col(0) )(0) * std::cos( ( measured_Y.col(0) )(1) * DgToRd ),
                startV,
                startK,
                1.0e-5
            };
            arma::vec startY = observationModel( startX );
            arma::vec deltaY = measured_Y.col(0) - startY;
            deltaY = checkDeltaMeasurement( deltaY );

            double dispR = resElementR * resElementR / 12.0;
            double sigmaV = ( Vmax - ( measured_Y.col(0) )(2) ); // ( Vmax - measured_Y_Vr[0] )
            double sigmaK = 2.0 * std::acos( v0 / Vmax ) * RdToDg;
            double dispV = sigmaV * sigmaV / 12.0;
            double dispK = sigmaK * sigmaK / 12.0;

            arma::vec startP {
                dispR * 1.0,
                dispR * 1.0,
                dispV * 0.25,
                dispK * 0.25,
                ( RMS_X_Ka * RMS_X_Ka ) * 1.0
            };
#ifdef EKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "EKF" ) ) {
                EKF.SetEstimatedVectorX( startX );
                EKF.SetEstimatedVectorY( startY );
                EKF.SetMeasuredVectorY( measured_Y.col(0) );
                EKF.SetDeltaY( deltaY );

                arma::mat PdenseEKF = arma::mat( SizeX, SizeX );
                PdenseEKF.fill( 1.0e-9 );
                PdenseEKF.diag() = startP;
                EKF.SetEstimateCovarianceMatrixP( PdenseEKF );
                if( settings.Debug ) {
                    PdenseEKF.print("PdenseEKF:");
                }
                estimated_X.at("EKF").col(0) = startX;
                estimated_Y.at("EKF").col(0) = startY;
                estimated_Pdiag.at("EKF").col(0) = arma::sqrt( PdenseEKF.diag() );
                SDCM.at("EKF")[0] = std::sqrt( arma::trace( PdenseEKF ) );
                RMSE_X.at("EKF").col(0) += arma::square( checkDeltaState( estimated_X.at("EKF").col(0) - true_X.col(0) ) );
            }
#endif
#ifdef UKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "UKF" ) ) {
                UKF.SetEstimatedVectorX( startX );
                UKF.SetEstimatedVectorY( startY );
                UKF.SetMeasuredVectorY( measured_Y.col(0) );
                UKF.SetDeltaY( deltaY );

                arma::mat PdenseUKF = arma::mat( SizeX, SizeX );
                PdenseUKF.fill( 1.0e-9 );
                PdenseUKF.diag() = startP ;
                UKF.SetEstimateCovarianceMatrixP( PdenseUKF );
                if( settings.Debug ) {
                    PdenseUKF.print("PdenseUKF:");
                }
                estimated_X.at("UKF").col(0) = startX;
                estimated_Y.at("UKF").col(0) = startY;
                estimated_Pdiag.at("UKF").col(0) = arma::sqrt( PdenseUKF.diag() );
                SDCM.at("UKF")[0] = std::sqrt( arma::trace( PdenseUKF ) );
                RMSE_X.at("UKF").col(0) += arma::square( checkDeltaState( estimated_X.at("UKF").col(0) - true_X.col(0) ) );
            }
#endif
#ifdef SRUKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                SRUKF.SetEstimatedVectorX( startX );
                SRUKF.SetEstimatedVectorY( startY );
                SRUKF.SetMeasuredVectorY( measured_Y.col(0) );
                SRUKF.SetDeltaY( deltaY );

                arma::mat PdenseSRUKF = arma::mat( SizeX, SizeX );
                PdenseSRUKF.fill( 1.0e-9 );
                PdenseSRUKF.diag() = startP;
                PdenseSRUKF = arma::chol( PdenseSRUKF, "lower" );
                SRUKF.SetEstimateCovarianceMatrixP( PdenseSRUKF );
                if( settings.Debug ) {
                    PdenseSRUKF.print("PdenseSRUKF:");
                }
                estimated_X.at("SRUKF").col(0) = startX;
                estimated_Y.at("SRUKF").col(0) = startY;
                estimated_Pdiag.at("SRUKF").col(0) = PdenseSRUKF.diag();
                SDCM.at("SRUKF")[0] = std::sqrt( arma::trace( PdenseSRUKF * arma::trans( PdenseSRUKF ) ) );
                RMSE_X.at("SRUKF").col(0) += arma::square( checkDeltaState( estimated_X.at("SRUKF").col(0) - true_X.col(0) ) );
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKFB" ) ) {
                SRUKFB.SetEstimatedVectorX( startX );
                SRUKFB.SetEstimatedVectorY( startY );
                SRUKFB.SetMeasuredVectorY( measured_Y.col(0) );
                SRUKFB.SetDeltaY( deltaY );

                arma::mat PdenseSRUKFB = arma::mat( SizeX, SizeX );
                PdenseSRUKFB.fill( 1.0e-9 );
                PdenseSRUKFB.diag() = startP;
                PdenseSRUKFB = arma::chol( PdenseSRUKFB, "lower" );
                SRUKFB.SetEstimateCovarianceMatrixP( PdenseSRUKFB );
                if( settings.Debug ) {
                    PdenseSRUKFB.print("PdenseSRUKFB:");
                }
                estimated_X.at("SRUKFB").col(0) = startX;
                estimated_Y.at("SRUKFB").col(0) = startY;
                estimated_Pdiag.at("SRUKFB").col(0) = PdenseSRUKFB.diag();
                SDCM.at("SRUKFB")[0] = std::sqrt( arma::trace( PdenseSRUKFB * arma::trans( PdenseSRUKFB ) ) );
                RMSE_X.at("SRUKFB").col(0) += arma::square( checkDeltaState( estimated_X.at("SRUKFB").col(0) - true_X.col(0) ) );
            }
#endif
#ifdef CKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "CKF" ) ) {
                CKF.SetEstimatedVectorX( startX );
                CKF.SetEstimatedVectorY( startY );
                CKF.SetMeasuredVectorY( measured_Y.col(0) );
                CKF.SetDeltaY( deltaY );

                arma::mat PdenseCKF = arma::mat( SizeX, SizeX );
                PdenseCKF.fill( 1.0e-9 );
                PdenseCKF.diag() = startP;
                CKF.SetEstimateCovarianceMatrixP( PdenseCKF );
                if( settings.Debug ) {
                    PdenseCKF.print("PdenseCKF:");
                }
                estimated_X.at("CKF").col(0) = startX;
                estimated_Y.at("CKF").col(0) = startY;
                estimated_Pdiag.at("CKF").col(0) = arma::sqrt( PdenseCKF.diag() );
                SDCM.at("CKF")[0] = std::sqrt( arma::trace( PdenseCKF ) );
                RMSE_X.at("CKF").col(0) += arma::square( checkDeltaState( estimated_X.at("CKF").col(0) - true_X.col(0) ) );
            }
#endif
#ifdef SRCKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
                SRCKF.SetEstimatedVectorX( startX );
                SRCKF.SetEstimatedVectorY( startY );
                SRCKF.SetMeasuredVectorY( measured_Y.col(0) );
                SRCKF.SetDeltaY( deltaY );

                arma::mat PdenseSRCKF = arma::mat( SizeX, SizeX );
                PdenseSRCKF.fill( 1.0e-9 );
                PdenseSRCKF.diag() = startP;
                PdenseSRCKF = arma::chol( PdenseSRCKF, "lower" );
                SRCKF.SetEstimateCovarianceMatrixP( PdenseSRCKF );
                if( settings.Debug ) {
                    PdenseSRCKF.print("PdenseSRCKF:");
                }
                estimated_X.at("SRCKF").col(0) = startX;
                estimated_Y.at("SRCKF").col(0) = startY;
                estimated_Pdiag.at("SRCKF").col(0) = PdenseSRCKF.diag();
                SDCM.at("SRCKF")[0] = std::sqrt( arma::trace( PdenseSRCKF * arma::trans( PdenseSRCKF ) ) );
                RMSE_X.at("SRCKF").col(0) += arma::square( checkDeltaState( estimated_X.at("SRCKF").col(0) - true_X.col(0) ) );
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKFB" ) ) {
                SRCKFB.SetEstimatedVectorX( startX );
                SRCKFB.SetEstimatedVectorY( startY );
                SRCKFB.SetMeasuredVectorY( measured_Y.col(0) );
                SRCKFB.SetDeltaY( deltaY );

                arma::mat PdenseSRCKFB = arma::mat( SizeX, SizeX );
                PdenseSRCKFB.fill( 1.0e-9 );
                PdenseSRCKFB.diag() = startP;
                PdenseSRCKFB = arma::chol( PdenseSRCKFB, "lower" );
                SRCKFB.SetEstimateCovarianceMatrixP( PdenseSRCKFB );
                if( settings.Debug ) {
                    PdenseSRCKFB.print("PdenseSRCKFB:");
                }
                estimated_X.at("SRCKFB").col(0) = startX;
                estimated_Y.at("SRCKFB").col(0) = startY;
                estimated_Pdiag.at("SRCKFB").col(0) = PdenseSRCKFB.diag();
                SDCM.at("SRCKFB")[0] = std::sqrt( arma::trace( PdenseSRCKFB * arma::trans( PdenseSRCKFB ) ) );
                RMSE_X.at("SRCKFB").col(0) += arma::square( checkDeltaState( estimated_X.at("SRCKFB").col(0) - true_X.col(0) ) );
            }
#endif
#ifdef ECKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "ECKF" ) ) {
                ECKF.SetEstimatedVectorX( startX );
                ECKF.SetEstimatedVectorY( startY );
                ECKF.SetMeasuredVectorY( measured_Y.col(0) );
                ECKF.SetDeltaY( deltaY );

                arma::mat PdenseECKF = arma::mat( SizeX, SizeX );
                PdenseECKF.fill( 1.0e-9 );
                PdenseECKF.diag() = startP;
                ECKF.SetEstimateCovarianceMatrixP( PdenseECKF );
                if( settings.Debug ) {
                    PdenseECKF.print("PdenseECKF:");
                }
                estimated_X.at("ECKF").col(0) = startX;
                estimated_Y.at("ECKF").col(0) = startY;
                estimated_Pdiag.at("ECKF").col(0) = arma::sqrt( PdenseECKF.diag() );
                SDCM.at("ECKF")[0] = std::sqrt( arma::trace( PdenseECKF ) );
                RMSE_X.at("ECKF").col(0) += arma::square( checkDeltaState( estimated_X.at("ECKF").col(0) - true_X.col(0) ) );
            }
#endif
#ifdef SRECKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
                SRECKF.SetEstimatedVectorX( startX );
                SRECKF.SetEstimatedVectorY( startY );
                SRECKF.SetMeasuredVectorY( measured_Y.col(0) );
                SRECKF.SetDeltaY( deltaY );

                arma::mat PdenseSRECKF = arma::mat( SizeX, SizeX );
                PdenseSRECKF.fill( 1.0e-9 );
                PdenseSRECKF.diag() = startP;
                PdenseSRECKF = arma::chol( PdenseSRECKF, "lower" );
                SRECKF.SetEstimateCovarianceMatrixP( PdenseSRECKF );
                if( settings.Debug ) {
                    PdenseSRECKF.print("PdenseSRECKF:");
                }
                estimated_X.at("SRECKF").col(0) = startX;
                estimated_Y.at("SRECKF").col(0) = startY;
                estimated_Pdiag.at("SRECKF").col(0) = PdenseSRECKF.diag();
                SDCM.at("SRECKF")[0] = std::sqrt( arma::trace( PdenseSRECKF * arma::trans( PdenseSRECKF ) ) );
                RMSE_X.at("SRECKF").col(0) += arma::square( checkDeltaState( estimated_X.at("SRECKF").col(0) - true_X.col(0) ) );
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKFB" ) ) {
                SRECKFB.SetEstimatedVectorX( startX );
                SRECKFB.SetEstimatedVectorY( startY );
                SRECKFB.SetMeasuredVectorY( measured_Y.col(0) );
                SRECKFB.SetDeltaY( deltaY );

                arma::mat PdenseSRECKFB = arma::mat( SizeX, SizeX );
                PdenseSRECKFB.fill( 1.0e-9 );
                PdenseSRECKFB.diag() = startP;
                PdenseSRECKFB = arma::chol( PdenseSRECKFB, "lower" );
                SRECKFB.SetEstimateCovarianceMatrixP( PdenseSRECKFB );
                if( settings.Debug ) {
                    PdenseSRECKFB.print("PdenseSRECKFB:");
                }
                estimated_X.at("SRECKFB").col(0) = startX;
                estimated_Y.at("SRECKFB").col(0) = startY;
                estimated_Pdiag.at("SRECKFB").col(0) = PdenseSRECKFB.diag();
                SDCM.at("SRECKFB")[0] = std::sqrt( arma::trace( PdenseSRECKFB * arma::trans( PdenseSRECKFB ) ) );
                RMSE_X.at("SRECKFB").col(0) += arma::square( checkDeltaState( estimated_X.at("SRECKFB").col(0) - true_X.col(0) ) );
            }
#endif
#ifdef EUKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "EUKF" ) ) {
                EUKF.SetEstimatedVectorX( startX );
                EUKF.SetEstimatedVectorY( startY );
                EUKF.SetMeasuredVectorY( measured_Y.col(0) );
                EUKF.SetDeltaY( deltaY );

                arma::mat PdenseEUKF = arma::mat( SizeX, SizeX );
                PdenseEUKF.fill( 1.0e-9 );
                PdenseEUKF.diag() = startP;
                EUKF.SetEstimateCovarianceMatrixP( PdenseEUKF );
                if( settings.Debug ) {
                    PdenseEUKF.print("PdenseEUKF:");
                }
                estimated_X.at("EUKF").col(0) = startX;
                estimated_Y.at("EUKF").col(0) = startY;
                estimated_Pdiag.at("EUKF").col(0) = arma::sqrt( PdenseEUKF.diag() );
                SDCM.at("EUKF")[0] = std::sqrt( arma::trace( PdenseEUKF ) );
                RMSE_X.at("EUKF").col(0) += arma::square( checkDeltaState( estimated_X.at("EUKF").col(0) - true_X.col(0) ) );
            }
#endif
#ifdef SREUKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                SREUKF.SetEstimatedVectorX( startX );
                SREUKF.SetEstimatedVectorY( startY );
                SREUKF.SetMeasuredVectorY( measured_Y.col(0) );
                SREUKF.SetDeltaY( deltaY );

                arma::mat PdenseSREUKF = arma::mat( SizeX, SizeX );
                PdenseSREUKF.fill( 1.0e-9 );
                PdenseSREUKF.diag() = startP;
                PdenseSREUKF = arma::chol( PdenseSREUKF, "lower" );
                SREUKF.SetEstimateCovarianceMatrixP( PdenseSREUKF );
                if( settings.Debug ) {
                    PdenseSREUKF.print("PdenseSREUKF:");
                }
                estimated_X.at("SREUKF").col(0) = startX;
                estimated_Y.at("SREUKF").col(0) = startY;
                estimated_Pdiag.at("SREUKF").col(0) = PdenseSREUKF.diag();
                SDCM.at("SREUKF")[0] = std::sqrt( arma::trace( PdenseSREUKF * arma::trans( PdenseSREUKF ) ) );
                RMSE_X.at("SREUKF").col(0) += arma::square( checkDeltaState( estimated_X.at("SREUKF").col(0) - true_X.col(0) ) );
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKFB" ) ) {
                SREUKFB.SetEstimatedVectorX( startX );
                SREUKFB.SetEstimatedVectorY( startY );
                SREUKFB.SetMeasuredVectorY( measured_Y.col(0) );
                SREUKFB.SetDeltaY( deltaY );

                arma::mat PdenseSREUKFB = arma::mat( SizeX, SizeX );
                PdenseSREUKFB.fill( 1.0e-9 );
                PdenseSREUKFB.diag() = startP;
                PdenseSREUKFB = arma::chol( PdenseSREUKFB, "lower" );
                SREUKFB.SetEstimateCovarianceMatrixP( PdenseSREUKFB );
                if( settings.Debug ) {
                    PdenseSREUKFB.print("PdenseSREUKFB:");
                }
                estimated_X.at("SREUKFB").col(0) = startX;
                estimated_Y.at("SREUKFB").col(0) = startY;
                estimated_Pdiag.at("SREUKFB").col(0) = PdenseSREUKFB.diag();
                SDCM.at("SREUKFB")[0] = std::sqrt( arma::trace( PdenseSREUKFB * arma::trans( PdenseSREUKFB ) ) );
                RMSE_X.at("SREUKFB").col(0) += arma::square( checkDeltaState( estimated_X.at("SREUKFB").col(0) - true_X.col(0) ) );
            }
#endif
#ifdef SREKF_
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREKF" ) ) {
                SREKF.SetEstimatedVectorX( startX );
                SREKF.SetEstimatedVectorY( startY );
                SREKF.SetMeasuredVectorY( measured_Y.col(0) );
                SREKF.SetDeltaY( deltaY );

                arma::mat PdenseSREKF = arma::mat( SizeX, SizeX );
                PdenseSREKF.fill( 1.0e-9 );
                PdenseSREKF.diag() = startP;
                PdenseSREKF = arma::chol( PdenseSREKF, "lower" );
                SREKF.SetEstimateCovarianceMatrixP( PdenseSREKF );
                if( settings.Debug ) {
                    PdenseSREKF.print("PdenseSREKF:");
                }
                estimated_X.at("SREKF").col(0) = startX;
                estimated_Y.at("SREKF").col(0) = startY;
                estimated_Pdiag.at("SREKF").col(0) = PdenseSREKF.diag();
                SDCM.at("SREKF")[0] = std::sqrt( arma::trace( PdenseSREKF * arma::trans( PdenseSREKF ) ) );
                RMSE_X.at("SREKF").col(0) += arma::square( checkDeltaState( estimated_X.at("SREKF").col(0) - true_X.col(0) ) );
            }
#endif
            //----------------------------------------------------------------------------------------------------------
            // Симуляция
            for( int i = 1; i < N; i++ ) { // Цикл по тактам времени
                time[i] = i * settings.DeltaT;
                if( settings.Debug ) {
                    std::cout << "Takt = " << i << " Time = " << time[i] << std::endl;
                }

                // Прогноз
#ifdef EKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "EKF" ) ) {
                    timerPrediction.at("EKF").StartTimer();
                    EKF.Prediction( settings.DeltaT );
                    timerPrediction.at("EKF").EndTimer();
                    timeMeasuredPrediction.at("EKF")[i] = timerPrediction.at("EKF").TimeCur() * 1.0e6; // мкс
                }
#endif
#ifdef UKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "UKF" ) ) {
                    timerPrediction.at("UKF").StartTimer();
                    UKF.Prediction( settings.DeltaT );
                    timerPrediction.at("UKF").EndTimer();
                    timeMeasuredPrediction.at("UKF")[i] = timerPrediction.at("UKF").TimeCur() * 1.0e6; // мкс
                }
#endif
#ifdef SRUKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                    timerPrediction.at("SRUKF").StartTimer();
                    SRUKF.Prediction( settings.DeltaT );
                    timerPrediction.at("SRUKF").EndTimer();
                    timeMeasuredPrediction.at("SRUKF")[i] = timerPrediction.at("SRUKF").TimeCur() * 1.0e6; // мкс
                }
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKFB" ) ) {
                    timerPrediction.at("SRUKFB").StartTimer();
                    SRUKFB.Prediction( settings.DeltaT );
                    timerPrediction.at("SRUKFB").EndTimer();
                    timeMeasuredPrediction.at("SRUKFB")[i] = timerPrediction.at("SRUKFB").TimeCur() * 1.0e6; // мкс
                }
#endif
#ifdef CKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "CKF" ) ) {
                    timerPrediction.at("CKF").StartTimer();
                    CKF.Prediction( settings.DeltaT );
                    timerPrediction.at("CKF").EndTimer();
                    timeMeasuredPrediction.at("CKF")[i] = timerPrediction.at("CKF").TimeCur() * 1.0e6; // мкс
                }
#endif
#ifdef SRCKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
                    timerPrediction.at("SRCKF").StartTimer();
                    SRCKF.Prediction( settings.DeltaT );
                    timerPrediction.at("SRCKF").EndTimer();
                    timeMeasuredPrediction.at("SRCKF")[i] = timerPrediction.at("SRCKF").TimeCur() * 1.0e6; // мкс
                }
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKFB" ) ) {
                    timerPrediction.at("SRCKFB").StartTimer();
                    SRCKFB.Prediction( settings.DeltaT );
                    timerPrediction.at("SRCKFB").EndTimer();
                    timeMeasuredPrediction.at("SRCKFB")[i] = timerPrediction.at("SRCKFB").TimeCur() * 1.0e6; // мкс
                }
#endif
#ifdef ECKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "ECKF" ) ) {
                    timerPrediction.at("ECKF").StartTimer();
                    ECKF.Prediction( settings.DeltaT );
                    timerPrediction.at("ECKF").EndTimer();
                    timeMeasuredPrediction.at("ECKF")[i] = timerPrediction.at("ECKF").TimeCur() * 1.0e6; // мкс
                }
#endif
#ifdef SRECKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
                    timerPrediction.at("SRECKF").StartTimer();
                    SRECKF.Prediction( settings.DeltaT );
                    timerPrediction.at("SRECKF").EndTimer();
                    timeMeasuredPrediction.at("SRECKF")[i] = timerPrediction.at("SRECKF").TimeCur() * 1.0e6; // мкс
                }
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKFB" ) ) {
                    timerPrediction.at("SRECKFB").StartTimer();
                    SRECKFB.Prediction( settings.DeltaT );
                    timerPrediction.at("SRECKFB").EndTimer();
                    timeMeasuredPrediction.at("SRECKFB")[i] = timerPrediction.at("SRECKFB").TimeCur() * 1.0e6; // мкс
                }
#endif
#ifdef EUKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "EUKF" ) ) {
                    timerPrediction.at("EUKF").StartTimer();
                    EUKF.Prediction( settings.DeltaT );
                    timerPrediction.at("EUKF").EndTimer();
                    timeMeasuredPrediction.at("EUKF")[i] = timerPrediction.at("EUKF").TimeCur() * 1.0e6; // мкс
                }
#endif
#ifdef SREUKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                    timerPrediction.at("SREUKF").StartTimer();
                    SREUKF.Prediction( settings.DeltaT );
                    timerPrediction.at("SREUKF").EndTimer();
                    timeMeasuredPrediction.at("SREUKF")[i] = timerPrediction.at("SREUKF").TimeCur() * 1.0e6; // мкс
                }
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKFB" ) ) {
                    timerPrediction.at("SREUKFB").StartTimer();
                    SREUKFB.Prediction( settings.DeltaT );
                    timerPrediction.at("SREUKFB").EndTimer();
                    timeMeasuredPrediction.at("SREUKFB")[i] = timerPrediction.at("SREUKFB").TimeCur() * 1.0e6; // мкс
                }
#endif
#ifdef SREKF_
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREKF" ) ) {
                    timerPrediction.at("SREKF").StartTimer();
                    SREKF.Prediction( settings.DeltaT );
                    timerPrediction.at("SREKF").EndTimer();
                    timeMeasuredPrediction.at("SREKF")[i] = timerPrediction.at("SREKF").TimeCur() * 1.0e6; // мкс
                }
#endif
                // Создание измерений текущего такта:
                true_X.col(i) = stateTransitionModel( true_X.col(i - 1), settings.DeltaT );

                // Маневр:
                double start_K_man = 2000.0;
                double end_K_man = start_K_man + 60 * 5;
//                double start_K_man = 4000.0;
//                double end_K_man = start_K_man + 60 * 15;
                if( time[i] > start_K_man && time[i] < end_K_man ) // Если манёвр
                {
                    ( true_X.col(i) )(4) = 0.6;//1.0;// // Ka
                }
                if( time[i] >= end_K_man ) { // Манёвр кончился
                    ( true_X.col(i) )(4) = 1.0e-5;//0.0;//
                }
                // Measurement Y
                true_Y.col(i) = observationModel( true_X.col(i) );
                measured_Y.col(i) = true_Y.col(i) +
                    arma::vec{ noiseY_R( generator ), noiseY_Az( generator ), noiseY_Vr( generator ) };

                // Коррекция
                double randomNum = random_0_1( generator ); // Случайное число от 0 до 1
                if( ( p_ > randomNum ) ) { // || ( i < 2 )
#ifdef EKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "EKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("EKF").col(i) = EKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_EKF = measured_Y.col(i) - estimated_Y.at("EKF").col(i);
                        deltaY_EKF = checkDeltaMeasurement( deltaY_EKF );

                        // Коррекция
                        timerCorrection.at("EKF").StartTimer();
                        EKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("EKF").EndTimer();
                        timeMeasuredCorrection.at("EKF")[i] = timerCorrection.at("EKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_EKF = EKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("EKF").col(i) = arma::sqrt( arma::abs( S_EKF.diag() ) );
                        arma::mat Sinv_EKF = arma::inv( S_EKF );
                        arma::mat lam_EKF = arma::trans( deltaY_EKF ) * Sinv_EKF * deltaY_EKF;
                        double md_EKF = std::sqrt( lam_EKF[0] );
                        mahalanobis.at("EKF")[i] = md_EKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
#ifdef UKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "UKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("UKF").col(i) = UKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_UKF = measured_Y.col(i) - estimated_Y.at("UKF").col(i);
                        deltaY_UKF = checkDeltaMeasurement( deltaY_UKF );

                        // Коррекция
                        timerCorrection.at("UKF").StartTimer();
                        UKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("UKF").EndTimer();
                        timeMeasuredCorrection.at("UKF")[i] = timerCorrection.at("UKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_UKF = UKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("UKF").col(i) = arma::sqrt( arma::abs( S_UKF.diag() ) );
                        arma::mat Sinv_UKF = arma::inv( S_UKF );
                        arma::mat lam_UKF = arma::trans( deltaY_UKF ) * Sinv_UKF * deltaY_UKF;
                        double md_UKF = std::sqrt( lam_UKF[0] );
                        mahalanobis.at("UKF")[i] = md_UKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
#ifdef SRUKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("SRUKF").col(i) = SRUKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_SRUKF = measured_Y.col(i) - estimated_Y.at("SRUKF").col(i);
                        deltaY_SRUKF = checkDeltaMeasurement( deltaY_SRUKF );

                        // Коррекция
                        timerCorrection.at("SRUKF").StartTimer();
                        SRUKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("SRUKF").EndTimer();
                        timeMeasuredCorrection.at("SRUKF")[i] = timerCorrection.at("SRUKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_SRUKF = SRUKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("SRUKF").col(i) = S_SRUKF.diag();
                        arma::mat Sinv_SRUKF = arma::inv( S_SRUKF * arma::trans( S_SRUKF ) );
                        arma::mat lam_SRUKF = arma::trans( deltaY_SRUKF ) * Sinv_SRUKF * deltaY_SRUKF;
                        double md_SRUKF = std::sqrt( lam_SRUKF[0] );
                        mahalanobis.at("SRUKF")[i] = md_SRUKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKFB" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("SRUKFB").col(i) = SRUKFB.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_SRUKFB = measured_Y.col(i) - estimated_Y.at("SRUKFB").col(i);
                        deltaY_SRUKFB = checkDeltaMeasurement( deltaY_SRUKFB );

                        // Коррекция
                        timerCorrection.at("SRUKFB").StartTimer();
                        SRUKFB.Correction( measured_Y.col(i) );
                        timerCorrection.at("SRUKFB").EndTimer();
                        timeMeasuredCorrection.at("SRUKFB")[i] = timerCorrection.at("SRUKFB").TimeCur() * 1.0e6; // мкс

                        arma::mat S_SRUKFB = SRUKFB.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("SRUKFB").col(i) = S_SRUKFB.diag();
                        arma::mat Sinv_SRUKFB = arma::inv( S_SRUKFB * arma::trans( S_SRUKFB ) );
                        arma::mat lam_SRUKFB = arma::trans( deltaY_SRUKFB ) * Sinv_SRUKFB * deltaY_SRUKFB;
                        double md_SRUKFB = std::sqrt( lam_SRUKFB[0] );
                        mahalanobis.at("SRUKFB")[i] = md_SRUKFB;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
#ifdef CKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "CKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("CKF").col(i) = CKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_CKF = measured_Y.col(i) - estimated_Y.at("CKF").col(i);
                        deltaY_CKF = checkDeltaMeasurement( deltaY_CKF );

                        // Коррекция
                        timerCorrection.at("CKF").StartTimer();
                        CKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("CKF").EndTimer();
                        timeMeasuredCorrection.at("CKF")[i] = timerCorrection.at("CKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_CKF = CKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("CKF").col(i) = arma::sqrt( arma::abs( S_CKF.diag() ) );
                        arma::mat Sinv_CKF = arma::inv( S_CKF );
                        arma::mat lam_CKF = arma::trans( deltaY_CKF ) * Sinv_CKF * deltaY_CKF;
                        double md_CKF = std::sqrt( lam_CKF[0] );
                        mahalanobis.at("CKF")[i] = md_CKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
#ifdef SRCKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("SRCKF").col(i) = SRCKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_SRCKF = measured_Y.col(i) - estimated_Y.at("SRCKF").col(i);
                        deltaY_SRCKF = checkDeltaMeasurement( deltaY_SRCKF );

                        // Коррекция
                        timerCorrection.at("SRCKF").StartTimer();
                        SRCKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("SRCKF").EndTimer();
                        timeMeasuredCorrection.at("SRCKF")[i] = timerCorrection.at("SRCKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_SRCKF = SRCKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("SRCKF").col(i) = S_SRCKF.diag();
                        arma::mat Sinv_SRCKF = arma::inv( S_SRCKF * arma::trans( S_SRCKF ) );
                        arma::mat lam_SRCKF = arma::trans( deltaY_SRCKF ) * Sinv_SRCKF * deltaY_SRCKF;
                        double md_SRCKF = std::sqrt( lam_SRCKF[0] );
                        mahalanobis.at("SRCKF")[i] = md_SRCKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKFB" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("SRCKFB").col(i) = SRCKFB.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_SRCKFB = measured_Y.col(i) - estimated_Y.at("SRCKFB").col(i);
                        deltaY_SRCKFB = checkDeltaMeasurement( deltaY_SRCKFB );

                        // Коррекция
                        timerCorrection.at("SRCKFB").StartTimer();
                        SRCKFB.Correction( measured_Y.col(i) );
                        timerCorrection.at("SRCKFB").EndTimer();
                        timeMeasuredCorrection.at("SRCKFB")[i] = timerCorrection.at("SRCKFB").TimeCur() * 1.0e6; // мкс

                        arma::mat S_SRCKFB = SRCKFB.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("SRCKFB").col(i) = S_SRCKFB.diag();
                        arma::mat Sinv_SRCKFB = arma::inv( S_SRCKFB * arma::trans( S_SRCKFB ) );
                        arma::mat lam_SRCKFB = arma::trans( deltaY_SRCKFB ) * Sinv_SRCKFB * deltaY_SRCKFB;
                        double md_SRCKFB = std::sqrt( lam_SRCKFB[0] );
                        mahalanobis.at("SRCKFB")[i] = md_SRCKFB;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
#ifdef ECKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "ECKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("ECKF").col(i) = ECKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_ECKF = measured_Y.col(i) - estimated_Y.at("ECKF").col(i);
                        deltaY_ECKF = checkDeltaMeasurement( deltaY_ECKF );

                        // Коррекция
                        timerCorrection.at("ECKF").StartTimer();
                        ECKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("ECKF").EndTimer();
                        timeMeasuredCorrection.at("ECKF")[i] = timerCorrection.at("ECKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_ECKF = ECKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("ECKF").col(i) = arma::sqrt( arma::abs( S_ECKF.diag() ) );
                        arma::mat Sinv_ECKF = arma::inv( S_ECKF );
                        arma::mat lam_ECKF = arma::trans( deltaY_ECKF ) * Sinv_ECKF * deltaY_ECKF;
                        double md_ECKF = std::sqrt( lam_ECKF[0] );
                        mahalanobis.at("ECKF")[i] = md_ECKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
#ifdef SRECKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("SRECKF").col(i) = SRECKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_SRECKF = measured_Y.col(i) - estimated_Y.at("SRECKF").col(i);
                        deltaY_SRECKF = checkDeltaMeasurement( deltaY_SRECKF );

                        // Коррекция
                        timerCorrection.at("SRECKF").StartTimer();
                        SRECKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("SRECKF").EndTimer();
                        timeMeasuredCorrection.at("SRECKF")[i] = timerCorrection.at("SRECKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_SRECKF = SRECKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("SRECKF").col(i) = S_SRECKF.diag();
                        arma::mat Sinv_SRECKF = arma::inv( S_SRECKF * arma::trans( S_SRECKF ) );
                        arma::mat lam_SRECKF = arma::trans( deltaY_SRECKF ) * Sinv_SRECKF * deltaY_SRECKF;
                        double md_SRECKF = std::sqrt( lam_SRECKF[0] );
                        mahalanobis.at("SRECKF")[i] = md_SRECKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKFB" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("SRECKFB").col(i) = SRECKFB.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_SRECKFB = measured_Y.col(i) - estimated_Y.at("SRECKFB").col(i);
                        deltaY_SRECKFB = checkDeltaMeasurement( deltaY_SRECKFB );

                        // Коррекция
                        timerCorrection.at("SRECKFB").StartTimer();
                        SRECKFB.Correction( measured_Y.col(i) );
                        timerCorrection.at("SRECKFB").EndTimer();
                        timeMeasuredCorrection.at("SRECKFB")[i] = timerCorrection.at("SRECKFB").TimeCur() * 1.0e6; // мкс

                        arma::mat S_SRECKFB = SRECKFB.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("SRECKFB").col(i) = S_SRECKFB.diag();
                        arma::mat Sinv_SRECKFB = arma::inv( S_SRECKFB * arma::trans( S_SRECKFB ) );
                        arma::mat lam_SRECKFB = arma::trans( deltaY_SRECKFB ) * Sinv_SRECKFB * deltaY_SRECKFB;
                        double md_SRECKFB = std::sqrt( lam_SRECKFB[0] );
                        mahalanobis.at("SRECKFB")[i] = md_SRECKFB;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
#ifdef EUKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "EUKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("EUKF").col(i) = EUKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_EUKF = measured_Y.col(i) - estimated_Y.at("EUKF").col(i);
                        deltaY_EUKF = checkDeltaMeasurement( deltaY_EUKF );

                        // Коррекция
                        timerCorrection.at("EUKF").StartTimer();
                        EUKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("EUKF").EndTimer();
                        timeMeasuredCorrection.at("EUKF")[i] = timerCorrection.at("EUKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_EUKF = EUKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("EUKF").col(i) = arma::sqrt( arma::abs( S_EUKF.diag() ) );
                        arma::mat Sinv_EUKF = arma::inv( S_EUKF );
                        arma::mat lam_EUKF = arma::trans( deltaY_EUKF ) * Sinv_EUKF * deltaY_EUKF;
                        double md_EUKF = std::sqrt( lam_EUKF[0] );
                        mahalanobis.at("EUKF")[i] = md_EUKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
#ifdef SREUKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("SREUKF").col(i) = SREUKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_SREUKF = measured_Y.col(i) - estimated_Y.at("SREUKF").col(i);
                        deltaY_SREUKF = checkDeltaMeasurement( deltaY_SREUKF );

                        // Коррекция
                        timerCorrection.at("SREUKF").StartTimer();
                        SREUKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("SREUKF").EndTimer();
                        timeMeasuredCorrection.at("SREUKF")[i] = timerCorrection.at("SREUKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_SREUKF = SREUKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("SREUKF").col(i) = S_SREUKF.diag();
                        arma::mat Sinv_SREUKF = arma::inv( S_SREUKF * arma::trans( S_SREUKF ) );
                        arma::mat lam_SREUKF = arma::trans( deltaY_SREUKF ) * Sinv_SREUKF * deltaY_SREUKF;
                        double md_SREUKF = std::sqrt( lam_SREUKF[0] );
                        mahalanobis.at("SREUKF")[i] = md_SREUKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKFB" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("SREUKFB").col(i) = SREUKFB.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_SREUKFB = measured_Y.col(i) - estimated_Y.at("SREUKFB").col(i);
                        deltaY_SREUKFB = checkDeltaMeasurement( deltaY_SREUKFB );

                        // Коррекция
                        timerCorrection.at("SREUKFB").StartTimer();
                        SREUKFB.Correction( measured_Y.col(i) );
                        timerCorrection.at("SREUKFB").EndTimer();
                        timeMeasuredCorrection.at("SREUKFB")[i] = timerCorrection.at("SREUKFB").TimeCur() * 1.0e6; // мкс

                        arma::mat S_SREUKFB = SREUKFB.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("SREUKFB").col(i) = S_SREUKFB.diag();
                        arma::mat Sinv_SREUKFB = arma::inv( S_SREUKFB * arma::trans( S_SREUKFB ) );
                        arma::mat lam_SREUKFB = arma::trans( deltaY_SREUKFB ) * Sinv_SREUKFB * deltaY_SREUKFB;
                        double md_SREUKFB = std::sqrt( lam_SREUKFB[0] );
                        mahalanobis.at("SREUKFB")[i] = md_SREUKFB;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
#ifdef SREKF_
                    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREKF" ) ) {
                        //----------------------------------------------------------------------------------------------
                        estimated_Y.at("SREKF").col(i) = SREKF.GetEstimatedVectorY(); // Прогнозный Y
                        arma::vec deltaY_SREKF = measured_Y.col(i) - estimated_Y.at("SREKF").col(i);
                        deltaY_SREKF = checkDeltaMeasurement( deltaY_SREKF );

                        // Коррекция
                        timerCorrection.at("SREKF").StartTimer();
                        SREKF.Correction( measured_Y.col(i) );
                        timerCorrection.at("SREKF").EndTimer();
                        timeMeasuredCorrection.at("SREKF")[i] = timerCorrection.at("SREKF").TimeCur() * 1.0e6; // мкс

                        arma::mat S_SREKF = SREKF.GetInnovationCovarianceMatrixS();
                        estimated_Sdiag.at("SREKF").col(i) = S_SREKF.diag();
                        arma::mat Sinv_SREKF = arma::inv( S_SREKF * arma::trans( S_SREKF ) );
                        arma::mat lam_SREKF = arma::trans( deltaY_SREKF ) * Sinv_SREKF * deltaY_SREKF;
                        double md_SREKF = std::sqrt( lam_SREKF[0] );
                        mahalanobis.at("SREKF")[i] = md_SREKF;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
#endif
                } // end if( ( p_ > randomNum ) ) - Correction
#ifdef EKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "EKF" ) ) {
                    timeMeasuredSumm.at("EKF")[i] = timeMeasuredPrediction.at("EKF")[i] + timeMeasuredCorrection.at("EKF")[i];
                    estimated_X.at("EKF").col(i) = EKF.GetEstimatedVectorX();
                    estimated_Y.at("EKF").col(i) = EKF.GetEstimatedVectorY();

                    arma::mat estimatedP_EKF = EKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("EKF").col(i) = arma::sqrt( estimatedP_EKF.diag() );
                    SDCM.at("EKF")[i] = std::sqrt( arma::trace( estimatedP_EKF ) );

                    delta_Y.at("EKF").col(i) = EKF.GetDeltaY();
                    RMSE_X.at("EKF").col(i) += arma::square( checkDeltaState( estimated_X.at("EKF").col(i) - true_X.col(i) ) );
                }
#endif
#ifdef UKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "UKF" ) ) {
                    timeMeasuredSumm.at("UKF")[i] = timeMeasuredPrediction.at("UKF")[i] + timeMeasuredCorrection.at("UKF")[i];
                    estimated_X.at("UKF").col(i) = UKF.GetEstimatedVectorX();
                    estimated_Y.at("UKF").col(i) = UKF.GetEstimatedVectorY();

                    arma::mat estimatedP_UKF = UKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("UKF").col(i) = arma::sqrt( estimatedP_UKF.diag() );
                    SDCM.at("UKF")[i] = std::sqrt( arma::trace( estimatedP_UKF ) );

                    delta_Y.at("UKF").col(i) = UKF.GetDeltaY();
                    RMSE_X.at("UKF").col(i) += arma::square( checkDeltaState( estimated_X.at("UKF").col(i) - true_X.col(i) ) );
                }
#endif
#ifdef SRUKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                    timeMeasuredSumm.at("SRUKF")[i] = timeMeasuredPrediction.at("SRUKF")[i] + timeMeasuredCorrection.at("SRUKF")[i];
                    estimated_X.at("SRUKF").col(i) = SRUKF.GetEstimatedVectorX();
                    estimated_Y.at("SRUKF").col(i) = SRUKF.GetEstimatedVectorY();

                    arma::mat estimatedP_SRUKF = SRUKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("SRUKF").col(i) = estimatedP_SRUKF.diag();
                    SDCM.at("SRUKF")[i] = std::sqrt( arma::trace( estimatedP_SRUKF * arma::trans( estimatedP_SRUKF ) ) );

                    delta_Y.at("SRUKF").col(i) = SRUKF.GetDeltaY();
                    RMSE_X.at("SRUKF").col(i) += arma::square( checkDeltaState( estimated_X.at("SRUKF").col(i) - true_X.col(i) ) );
                }
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKFB" ) ) {
                    timeMeasuredSumm.at("SRUKFB")[i] = timeMeasuredPrediction.at("SRUKFB")[i] + timeMeasuredCorrection.at("SRUKFB")[i];
                    estimated_X.at("SRUKFB").col(i) = SRUKFB.GetEstimatedVectorX();
                    estimated_Y.at("SRUKFB").col(i) = SRUKFB.GetEstimatedVectorY();

                    arma::mat estimatedP_SRUKFB = SRUKFB.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("SRUKFB").col(i) = estimatedP_SRUKFB.diag();
                    SDCM.at("SRUKFB")[i] = std::sqrt( arma::trace( estimatedP_SRUKFB * arma::trans( estimatedP_SRUKFB ) ) );

                    delta_Y.at("SRUKFB").col(i) = SRUKFB.GetDeltaY();
                    RMSE_X.at("SRUKFB").col(i) += arma::square( checkDeltaState( estimated_X.at("SRUKFB").col(i) - true_X.col(i) ) );
                }
#endif
#ifdef CKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "CKF" ) ) {
                    timeMeasuredSumm.at("CKF")[i] = timeMeasuredPrediction.at("CKF")[i] + timeMeasuredCorrection.at("CKF")[i];
                    estimated_X.at("CKF").col(i) = CKF.GetEstimatedVectorX();
                    estimated_Y.at("CKF").col(i) = CKF.GetEstimatedVectorY();

                    arma::mat estimatedP_CKF = CKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("CKF").col(i) = arma::sqrt( estimatedP_CKF.diag() );
                    SDCM.at("CKF")[i] = std::sqrt( arma::trace( estimatedP_CKF ) );

                    delta_Y.at("CKF").col(i) = CKF.GetDeltaY();
                    RMSE_X.at("CKF").col(i) += arma::square( checkDeltaState( estimated_X.at("CKF").col(i) - true_X.col(i) ) );
                }
#endif
#ifdef SRCKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
                    timeMeasuredSumm.at("SRCKF")[i] = timeMeasuredPrediction.at("SRCKF")[i] + timeMeasuredCorrection.at("SRCKF")[i];
                    estimated_X.at("SRCKF").col(i) = SRCKF.GetEstimatedVectorX();
                    estimated_Y.at("SRCKF").col(i) = SRCKF.GetEstimatedVectorY();

                    arma::mat estimatedP_SRCKF = SRCKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("SRCKF").col(i) = estimatedP_SRCKF.diag();
                    SDCM.at("SRCKF")[i] = std::sqrt( arma::trace( estimatedP_SRCKF * arma::trans( estimatedP_SRCKF ) ) );

                    delta_Y.at("SRCKF").col(i) = SRCKF.GetDeltaY();
                    RMSE_X.at("SRCKF").col(i) += arma::square( checkDeltaState( estimated_X.at("SRCKF").col(i) - true_X.col(i) ) );
                }
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKFB" ) ) {
                    timeMeasuredSumm.at("SRCKFB")[i] = timeMeasuredPrediction.at("SRCKFB")[i] + timeMeasuredCorrection.at("SRCKFB")[i];
                    estimated_X.at("SRCKFB").col(i) = SRCKFB.GetEstimatedVectorX();
                    estimated_Y.at("SRCKFB").col(i) = SRCKFB.GetEstimatedVectorY();

                    arma::mat estimatedP_SRCKFB = SRCKFB.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("SRCKFB").col(i) = estimatedP_SRCKFB.diag();
                    SDCM.at("SRCKFB")[i] = std::sqrt( arma::trace( estimatedP_SRCKFB * arma::trans( estimatedP_SRCKFB ) ) );

                    delta_Y.at("SRCKFB").col(i) = SRCKFB.GetDeltaY();
                    RMSE_X.at("SRCKFB").col(i) += arma::square( checkDeltaState( estimated_X.at("SRCKFB").col(i) - true_X.col(i) ) );
                }
#endif
#ifdef ECKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "ECKF" ) ) {
                    timeMeasuredSumm.at("ECKF")[i] = timeMeasuredPrediction.at("ECKF")[i] + timeMeasuredCorrection.at("ECKF")[i];
                    estimated_X.at("ECKF").col(i) = ECKF.GetEstimatedVectorX();
                    estimated_Y.at("ECKF").col(i) = ECKF.GetEstimatedVectorY();

                    arma::mat estimatedP_ECKF = ECKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("ECKF").col(i) = arma::sqrt( estimatedP_ECKF.diag() );
                    SDCM.at("ECKF")[i] = std::sqrt( arma::trace( estimatedP_ECKF ) );

                    delta_Y.at("ECKF").col(i) = ECKF.GetDeltaY();
                    RMSE_X.at("ECKF").col(i) += arma::square( checkDeltaState( estimated_X.at("ECKF").col(i) - true_X.col(i) ) );
                }
#endif
#ifdef SRECKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
                    timeMeasuredSumm.at("SRECKF")[i] = timeMeasuredPrediction.at("SRECKF")[i] + timeMeasuredCorrection.at("SRECKF")[i];
                    estimated_X.at("SRECKF").col(i) = SRECKF.GetEstimatedVectorX();
                    estimated_Y.at("SRECKF").col(i) = SRECKF.GetEstimatedVectorY();

                    arma::mat estimatedP_SRECKF = SRECKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("SRECKF").col(i) = estimatedP_SRECKF.diag();
                    SDCM.at("SRECKF")[i] = std::sqrt( arma::trace( estimatedP_SRECKF * arma::trans( estimatedP_SRECKF ) ) );

                    delta_Y.at("SRECKF").col(i) = SRECKF.GetDeltaY();
                    RMSE_X.at("SRECKF").col(i) += arma::square( checkDeltaState( estimated_X.at("SRECKF").col(i) - true_X.col(i) ) );
                }
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKFB" ) ) {
                    timeMeasuredSumm.at("SRECKFB")[i] = timeMeasuredPrediction.at("SRECKFB")[i] + timeMeasuredCorrection.at("SRECKFB")[i];
                    estimated_X.at("SRECKFB").col(i) = SRECKFB.GetEstimatedVectorX();
                    estimated_Y.at("SRECKFB").col(i) = SRECKFB.GetEstimatedVectorY();

                    arma::mat estimatedP_SRECKFB = SRECKFB.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("SRECKFB").col(i) = estimatedP_SRECKFB.diag();
                    SDCM.at("SRECKFB")[i] = std::sqrt( arma::trace( estimatedP_SRECKFB * arma::trans( estimatedP_SRECKFB ) ) );

                    delta_Y.at("SRECKFB").col(i) = SRECKFB.GetDeltaY();
                    RMSE_X.at("SRECKFB").col(i) += arma::square( checkDeltaState( estimated_X.at("SRECKFB").col(i) - true_X.col(i) ) );
                }
#endif
#ifdef EUKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "EUKF" ) ) {
                    timeMeasuredSumm.at("EUKF")[i] = timeMeasuredPrediction.at("EUKF")[i] + timeMeasuredCorrection.at("EUKF")[i];
                    estimated_X.at("EUKF").col(i) = EUKF.GetEstimatedVectorX();
                    estimated_Y.at("EUKF").col(i) = EUKF.GetEstimatedVectorY();

                    arma::mat estimatedP_EUKF = EUKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("EUKF").col(i) = arma::sqrt( estimatedP_EUKF.diag() );
                    SDCM.at("EUKF")[i] = std::sqrt( arma::trace( estimatedP_EUKF ) );

                    delta_Y.at("EUKF").col(i) = EUKF.GetDeltaY();
                    RMSE_X.at("EUKF").col(i) += arma::square( checkDeltaState( estimated_X.at("EUKF").col(i) - true_X.col(i) ) );
                }
#endif
#ifdef SREUKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                    timeMeasuredSumm.at("SREUKF")[i] = timeMeasuredPrediction.at("SREUKF")[i] + timeMeasuredCorrection.at("SREUKF")[i];
                    estimated_X.at("SREUKF").col(i) = SREUKF.GetEstimatedVectorX();
                    estimated_Y.at("SREUKF").col(i) = SREUKF.GetEstimatedVectorY();

                    arma::mat estimatedP_SREUKF = SREUKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("SREUKF").col(i) = estimatedP_SREUKF.diag();
                    SDCM.at("SREUKF")[i] = std::sqrt( arma::trace( estimatedP_SREUKF * arma::trans( estimatedP_SREUKF ) ) );

                    delta_Y.at("SREUKF").col(i) = SREUKF.GetDeltaY();
                    RMSE_X.at("SREUKF").col(i) += arma::square( checkDeltaState( estimated_X.at("SREUKF").col(i) - true_X.col(i) ) );
                }
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKFB" ) ) {
                    timeMeasuredSumm.at("SREUKFB")[i] = timeMeasuredPrediction.at("SREUKFB")[i] + timeMeasuredCorrection.at("SREUKFB")[i];
                    estimated_X.at("SREUKFB").col(i) = SREUKFB.GetEstimatedVectorX();
                    estimated_Y.at("SREUKFB").col(i) = SREUKFB.GetEstimatedVectorY();

                    arma::mat estimatedP_SREUKFB = SREUKFB.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("SREUKFB").col(i) = estimatedP_SREUKFB.diag();
                    SDCM.at("SREUKFB")[i] = std::sqrt( arma::trace( estimatedP_SREUKFB * arma::trans( estimatedP_SREUKFB ) ) );

                    delta_Y.at("SREUKFB").col(i) = SREUKFB.GetDeltaY();
                    RMSE_X.at("SREUKFB").col(i) += arma::square( checkDeltaState( estimated_X.at("SREUKFB").col(i) - true_X.col(i) ) );
                }
#endif
#ifdef SREKF_
                //------------------------------------------------------------------------------------------------------
                if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREKF" ) ) {
                    timeMeasuredSumm.at("SREKF")[i] = timeMeasuredPrediction.at("SREKF")[i] + timeMeasuredCorrection.at("SREKF")[i];
                    estimated_X.at("SREKF").col(i) = SREKF.GetEstimatedVectorX();
                    estimated_Y.at("SREKF").col(i) = SREKF.GetEstimatedVectorY();

                    arma::mat estimatedP_SREKF = SREKF.GetEstimatedCovarianceMatrixP();
                    estimated_Pdiag.at("SREKF").col(i) = estimatedP_SREKF.diag();
                    SDCM.at("SREKF")[i] = std::sqrt( arma::trace( estimatedP_SREKF * arma::trans( estimatedP_SREKF ) ) );

                    delta_Y.at("SREKF").col(i) = SREKF.GetDeltaY();
                    RMSE_X.at("SREKF").col(i) += arma::square( checkDeltaState( estimated_X.at("SREKF").col(i) - true_X.col(i) ) );
                }
#endif
            } // end Цикл по тактам времени
            //----------------------------------------------------------------------------------------------------------
            // Суммирование времени для усреднения
            for( auto &filter : settings.Filters ) {
                timerPredictionAverage.at( filter ) += timerPrediction.at( filter ).TimePerOp() * 1.0e6; // Прогноз
                timerCorrectionAverage.at( filter ) += timerCorrection.at( filter ).TimePerOp() * 1.0e6; // Коррекция
            }
            //----------------------------------------------------------------------------------------------------------
            // Построение графиков
            //----------------------------------------------------------------------------------------------------------
            if( settings.MCruns == 1 ) {
                std::cout << "\nПостроение графиков..." << std::endl;
            }
            if( settings.Graphs_0_RMSE_1 == 0 ) {
            if( settings.Debug ) {
                std::cout << "seed = " << seed_ << "/" << settings.MCruns << " P = " << p_ << std::endl;
            }
            //----------------------------------------------------------------------------------------------------------
            // X-Y, R-Az            
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::figure_size( pixels_width, pixels_height );
                plt::subplot( 1, 2, 1 );
            }
            plt::title( "Координаты X-Y" );
            plt::xlabel( "X, [км]" );
            plt::ylabel( "Y, [км]" );
            plt::plot(
                arma::conv_to< std::vector<double> >::from( true_X.row(0) ), // X
                arma::conv_to< std::vector<double> >::from( true_X.row(1) ), // Y
                true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot(
                    arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row(0) ), // X
                    arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row(1) ), // Y
                    estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
            plt::xlim( 98.0, 103.0 );
            plt::ylim( 197.0, 203.0 );
        //    plt::xlim( 99.5, 100.4 );
        //    plt::ylim( 199.9, 200.4 );
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_X_" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {                
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 1, 2, 2 );
            }
            plt::title( "Координаты R-Az" );
            plt::xlabel( "Az, [градусы]" );
            plt::ylabel( "R, [км]" );
            plt::plot(
                arma::conv_to< std::vector<double> >::from( measured_Y.row(1) ), // Az
                arma::conv_to< std::vector<double> >::from( measured_Y.row(0) ), // R
                marks_keywords );
            plt::plot(
                arma::conv_to< std::vector<double> >::from( true_Y.row(1) ), // Az
                arma::conv_to< std::vector<double> >::from( true_Y.row(0) ), // R
                true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot(
                    arma::conv_to< std::vector<double> >::from( ( estimated_Y.at( filter ) ).row(1) ), // Az
                        arma::conv_to< std::vector<double> >::from( ( estimated_Y.at( filter ) ).row(0) ), // R
                        estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
            plt::xlim( 26.4, 27.0 );
            plt::ylim( 220.0, 227.5 );
        //    plt::xlim( 26.4, 27.0 );
        //    plt::ylim( 220.0, 227.5 );
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_Y_" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( !settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_1_X-Y" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //----------------------------------------------------------------------------------------------------------
            // mahalanobis
            plt::figure_size( pixels_width, pixels_height );

            plt::title( "Расстояние Махаланобиса" );
            plt::xlabel( "Время, [с]" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, mahalanobis.at( filter ), { { "color", graphColors.at( filter ) }, { "linestyle", "-" },
                    { "linewidth", "2" }, { "label", filter } } );
            }
            plt::axhline( 3.5, 0, settings.SimulationTime - settings.DeltaT, { { "color", "black" }, { "linestyle", "--" },
                { "linewidth", "2" }, { "label", "Порог по уровню вероятности ~0.99" } } );
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.0, 6.0 ); // 7.0
            }
            plt::subplots_adjust( keywords );
            name_tmp = "./" + filters_names + "_2_mahalanobis" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );

//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //----------------------------------------------------------------------------------------------------------
            // State            
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::figure_size( pixels_width, pixels_height );
                plt::subplot( 3, 2, 1 );
            }
            plt::title( "Оценка координаты X" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "X, [км]");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row(0) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row(0) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
        //        plt::ylim( 98.0, 103.0 );
            if( ylim_yes ) {
                plt::ylim( 99.7, 101.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State_X" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 3 );
            }
            plt::title( "Оценка координаты Y");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "Y, [км]");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row(1) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row(1) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
        //        plt::ylim( 197.0, 203.0 );
            if( ylim_yes ) {
                plt::ylim( 199.5, 201.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State_Y" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 2 );
            }
            plt::title( "Оценка полной скорости");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "V, [м/с]");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row(2) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row(2) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
        //        plt::ylim( 94.0, 130.0 );
            if( ylim_yes ) {
                plt::ylim( 94.0, 120.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State_V" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 4 );
            }
            plt::title( "Оценка курса");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "К, [градусы]");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row(3) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row(3) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( -1.0, 60.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State_K" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 6 );
            }
            plt::title( "Оценка скорости изменения курса" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "dK/dt, [градусы/с]");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row(4) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row(4) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true) ;
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( -0.2, 0.3 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State_Ka" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( !settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //----------------------------------------------------------------------------------------------------------
            // P            
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::figure_size( pixels_width, pixels_height );
                plt::subplot( 3, 2, 1 );
            }
            plt::title( "СКО координаты X");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "СКО X, [км]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row(0) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.01, 0.03 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P_X" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 3 );
            }
            plt::title( "СКО координаты Y");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "СКО Y, [км]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row(1) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.01, 0.03 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P_Y" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 2 );
            }
            plt::title( "СКО полной скорости");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "СКО V, [м/с]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row(2) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.0, 4.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P_V" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 4 );
            }
            plt::title( "СКО курса");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "СКО K, [градусы]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row(3) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.0, 12.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P_K" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 6 );
            }
            plt::title("СКО скорости изменения курса");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "СКО dK/dt, [градусы/с]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row(4) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.1, 0.3 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P_Ka" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( !settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //----------------------------------------------------------------------------------------------------------
            // Measurement
            plt::figure_size( pixels_width, pixels_height );

            plt::subplot( 3, 1, 1 );
            plt::title( "Дальность отметки" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "R, [км]" );
            plt::plot( time, arma::conv_to< std::vector<double> >::from( measured_Y.row(0) ), marks_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Y.at( filter ).row(0) ), estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );

            plt::subplot( 3, 1, 2 );
            plt::title( "Азимут отметки" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "Az, [градусы]");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( measured_Y.row(1) ), marks_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Y.at( filter ).row(1) ), estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );

            plt::subplot( 3, 1, 3 );
            plt::title( "Радиальная скорость отметки");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "Vr, [м/с]");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( measured_Y.row(2) ), marks_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Y.at( filter ).row(2) ), estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );

            plt::subplots_adjust( keywords );
            name_tmp = "./" + filters_names + "_5_Measurement" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );

//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //----------------------------------------------------------------------------------------------------------
            // S matrix
            plt::figure_size( pixels_width, pixels_height );

            plt::subplot( 3, 1, 1 );
            plt::title( "S, корень диаг. элемент дальности" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "dR, [км]" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Sdiag.at( filter ).row(0) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );

            plt::subplot( 3, 1, 2 );
            plt::title( "S, корень диаг. элемент азимута" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "dAz, [градусы]" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Sdiag.at( filter ).row(1) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );

            plt::subplot( 3, 1, 3 );
            plt::title( "S, корень диаг. элемент рад.скорости" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "dVr, [м/с]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Sdiag.at( filter ).row(2) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );

            plt::subplots_adjust( keywords );
            name_tmp = "./" + filters_names + "_6_sqrtS" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );
//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //--------------------------------------------------------
            // deltaY
            plt::figure_size( pixels_width, pixels_height );

            plt::subplot( 3, 1, 1 );
            plt::title( "Невязка дальности" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "dR, [км]" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( delta_Y.at( filter ).row(0) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );

            plt::subplot( 3, 1, 2 );
            plt::title( "Невязка азимута" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "dAz, [градусы]" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( delta_Y.at( filter ).row(1) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );

            plt::subplot( 3, 1, 3 );
            plt::title( "Невязка радиальной скорости" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "dVr, [м/с]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( delta_Y.at( filter ).row(2) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );

            plt::subplots_adjust( keywords );
            name_tmp = "./" + filters_names + "_7_deltaY" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );
//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //--------------------------------------------------------
            // SDCM
            plt::figure_size( pixels_width, pixels_height );

            plt::title("SDCM");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "SDCM" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, SDCM.at( filter ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            name_tmp = "./" + filters_names + "_8_SDCM" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );
//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();

            //--------------------------------------------------------
            // Time
            plt::figure_size( pixels_width, pixels_height );
            plt::title("Время выполнения прогноза/коррекции/суммарно");
            plt::ylabel( "Время, [мкс]" );

            std::map<std::string, std::vector<double>> y;
            std::map<std::string, std::vector<double>> x;
            std::vector<double> x_ticks;

            double dy = settings.Filters.size() + 1;
            double fi = 0.0;
            for( auto &filter : settings.Filters ) {
                y.insert( std::make_pair( filter, std::vector<double>{
                    timerPrediction.at( filter ).TimePerOp() * 1.0e6, // Прогноз
                    timerCorrection.at( filter ).TimePerOp() * 1.0e6, // Коррекция
                    ( timerPrediction.at( filter ).TimePerOp() + timerCorrection.at( filter ).TimePerOp() ) * 1.0e6 // Суммарно
                } ) );
                x.insert( std::make_pair( filter, std::vector<double>{
                    ( 0.0 + fi ),       // Прогноз
                    ( dy + fi ),        // Коррекция
                    ( 2.0 * dy ) + fi   // Суммарно
                } ) );
                plt::bar( x.at( filter ), y.at( filter ), "black", "-", 1.0, estimated_keywords_time.at( filter ) );
                fi += 1.0;
            }
            double first_tick = ( settings.Filters.size() - 1.0 ) / 2.0;
            for( int t = 0; t < 3; t++ ) { // 3 группы - Прогноз, Коррекция, Суммарно
                x_ticks.push_back( first_tick + ( t * dy ) );
            }
            plt::xticks( x_ticks, { "Экстраполяция", "Коррекция", "Суммарно" } );

            plt::legend( legend_loc_time );
            plt::grid( true );
            name_tmp = "./" + filters_names + "_9_time" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );

            if( settings.ShowGraphs ) {
                plt::show();
            }
            plt::close();

            } //if( settings.Graphs_0_RMSE_1 == 0 )

        } // seed
        //--------------------------------------------------------------------------------------------------------------
        // RMSE
        if( settings.Graphs_0_RMSE_1 == 1 ) {

            std::cout << "\nПостроение графиков..." << std::endl;

            if( settings.Debug ) {
                std::cout << "P = " << p_ << std::endl;
            }

            for( auto &filter : settings.Filters ) {
                RMSE_X.at( filter ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
                RMSE_X.at( filter ) = arma::sqrt( RMSE_X.at( filter ) );
            }            
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::figure_size( pixels_width, pixels_height );
                plt::subplot( 3, 2, 1 );
            }
            plt::title( "RMSE координаты X" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "RMSE X, [км]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row(0) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.0, 0.03 ); // 0.05
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 3 );
            }
            plt::title( "RMSE координаты Y");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "RMSE Y, [км]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row(1) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.0, 0.03 ); // 0.05
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_Y_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 2 );
            }
            plt::title( "RMSE полной скорости");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "RMSE V, [м/с]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row(2) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.0, 10.0 ); // 20.0
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_V_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 4 );
            }
            plt::title( "RMSE курса");
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "RMSE К, [градусы]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row(3) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.0, 25.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_K_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 6 );
            }
            plt::title( "RMSE скорости изменения курса" );
            plt::xlabel( "Время, [с]" );
            plt::ylabel( "RMSE dK/dt, [градусы/с]");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row(4) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true) ;
            plt::xlim( 0.0, settings.SimulationTime );
            if( ylim_yes ) {
                plt::ylim( 0.0, 0.18 ); // 0.7
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_Ka_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( !settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //----------------------------------------------------------------------------------------------------------
            // Время усреднённое
            for( auto &filter : settings.Filters ) {
                timerPredictionAverage.at( filter ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
                timerCorrectionAverage.at( filter ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
            }
            plt::figure_size( pixels_width, pixels_height );
            std::string my_title = "Время выполнения прогноза/коррекции/суммарно усреднённое по " + std::to_string( settings.MCruns ) + " реализациям";
//            plt::title( my_title );
            plt::ylabel( "Время, [мкс]" );

            std::map<std::string, std::vector<double>> y;
            std::map<std::string, std::vector<double>> x;
            std::vector<double> x_ticks;

            double dy = settings.Filters.size() + 1;
            double fi = 0.0;
            for( auto &filter : settings.Filters ) {
                y.insert( std::make_pair( filter, std::vector<double>{
                    timerPredictionAverage.at( filter ), // Прогноз
                    timerCorrectionAverage.at( filter ), // Коррекция
                    ( timerPredictionAverage.at( filter ) + timerCorrectionAverage.at( filter ) ) // Суммарно
                } ) );
                x.insert( std::make_pair( filter, std::vector<double>{
                    ( 0.0 + fi ),       // Прогноз
                    ( dy + fi ),        // Коррекция
                    ( 2.0 * dy ) + fi   // Суммарно
                } ) );
                plt::bar( x.at( filter ), y.at( filter ), "black", "-", 1.0, estimated_keywords_time.at( filter ) );
                fi += 1.0;
            }
            double first_tick = ( settings.Filters.size() - 1.0 ) / 2.0;
            for( int t = 0; t < 3; t++ ) { // 3 группы - Прогноз, Коррекция, Суммарно
                x_ticks.push_back( first_tick + ( t * dy ) );
            }
            plt::xticks( x_ticks, { "Экстраполяция", "Коррекция", "Суммарно" } );

//            plt::legend( legend_loc_time );
            plt::legend( loc_time, { 1, 1 } );
            plt::grid( true );
            plt::subplots_adjust( keywords );
            name_tmp = "./" + filters_names + "_time_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
            plt::save( name_tmp, dpi );

//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //----------------------------------------------------------------------------------------------------------
            // Время усреднённое относительно EKF
            plt::figure_size( pixels_width, pixels_height );
            my_title = "Время выполнения прогноза/коррекции/суммарно (относительно EKF) усреднённое по " + std::to_string( settings.MCruns ) + " реализациям";
            name_tmp = "./" + filters_names + "_timeRelative_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
            std::string name_tmp_txt = "./" + filters_names + "_timeRelative_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + ".xls";
//            plt::title( my_title );
            plt::ylabel( "разы" );

            x.clear();
            y.clear();
            x_ticks.clear();

            dy = settings.Filters.size() + 1;
            fi = 0.0;

            std::ofstream ofs;
            ofs.open( name_tmp_txt, std::ofstream::out ); // txt
            ofs <<
                "Фильтр" << "\t" <<
                "Время прогноза относительно EKF" << "\t" <<
                "Время коррекции относительно EKF" << "\t" <<
                "Суммарное время коррекции относительно EKF" << "\t" <<
                std::endl;

            for( auto &filter : settings.Filters ) {
                if( filter == "EKF" ) {
                    continue;
                }
                double pred = ( timerPredictionAverage.at( filter ) / timerPredictionAverage.at( "EKF" ) );
                double corr = ( timerCorrectionAverage.at( filter ) / timerCorrectionAverage.at( "EKF" ) );
                double summ = ( timerPredictionAverage.at( filter ) + timerCorrectionAverage.at( filter ) ) /
                    ( timerPredictionAverage.at( "EKF" ) + timerCorrectionAverage.at( "EKF" ) );

                ofs <<
                    filter << "\t" <<
                    pred << "\t" <<
                    corr << "\t" <<
                    summ << "\t" <<
                    std::endl;

                y.insert( std::make_pair( filter, std::vector<double>{
                    pred, // Прогноз
                    corr, // Коррекция
                    summ // Суммарно
                } ) );
                x.insert( std::make_pair( filter, std::vector<double>{
                    ( 0.0 + fi ),       // Прогноз
                    ( dy + fi ),        // Коррекция
                    ( 2.0 * dy ) + fi   // Суммарно
                } ) );
                plt::bar( x.at( filter ), y.at( filter ), "black", "-", 1.0, estimated_keywords_time.at( filter ) );
                fi += 1.0;
            }

            plt::axhline( 1.0, 0.0, 1.0, { { "color", "black" }, { "linestyle", "-" },
                { "linewidth", "2" }, { "label", "EKF" } } );
    //        plt::ylim( 0.0, 10.0 );

            if( ofs.is_open() ) {
                ofs.close();
            }

            first_tick = ( settings.Filters.size() - 1.0 ) / 2.0;
            for( int t = 0; t < 3; t++ ) { // 3 группы - Прогноз, Коррекция, Суммарно
                x_ticks.push_back( first_tick + ( t * dy ) );
            }
            plt::xticks( x_ticks, { "Экстраполяция", "Коррекция", "Суммарно" } );

//            plt::legend( legend_loc_time );
            plt::legend( loc_time, { 1, 1 } );
            plt::grid( true );
            plt::subplots_adjust( keywords );
            plt::save( name_tmp, dpi );

            if( settings.ShowGraphs ) {
                plt::show();
            }
            plt::close();
        } // end if( settings.Graphs_0_RMSE_1 == 1 ) {
    } // p - по вероятностям
    //------------------------------------------------------------------------------------------------------------------
    plt::clf();
    plt::cla();
    Py_Finalize();

    }
    catch( std::exception &ex )
    {
        std::cout << "Exception occured: " << ex.what() << std::endl;
    }
    std::cout << "Завершено." << std::endl;
}

void CKalmanFiltersCompare::Run_RMSE_EKF_SRUKF_SREUKF( const CSettings &settings )
{
    try{

    Py_Initialize();

    namespace plt = matplotlibcpp;
    std::map<std::string, double> keywords;
//    keywords.insert( std::make_pair( "left", 0.08 ) );
//    keywords.insert( std::make_pair( "bottom", 0.08 ) );
//    keywords.insert( std::make_pair( "right", 0.92 ) );
//    keywords.insert( std::make_pair( "top", 0.92 ) );
//    keywords.insert( std::make_pair( "wspace", 0.4 ) );
//    keywords.insert( std::make_pair( "hspace", 0.4 ) ); // 0.7
    keywords.insert( std::make_pair( "left", settings.MatPlotParams[0] ) );
    keywords.insert( std::make_pair( "bottom", settings.MatPlotParams[1] ) );
    keywords.insert( std::make_pair( "right", settings.MatPlotParams[2] ) );
    keywords.insert( std::make_pair( "top", settings.MatPlotParams[3] ) );
    keywords.insert( std::make_pair( "wspace", settings.MatPlotParams[4] ) );
    keywords.insert( std::make_pair( "hspace", settings.MatPlotParams[5] ) );

    const double in2mm = 25.4;// mm (fixed)
//    const double pt2mm = 0.3528;// mm (fixed)

    const double dpi = 300;// dpi (variable)
//    double koef = 1.7; // 1.5; //
//    const double width = 195 * pt2mm * koef;// 250 ## mm (variable)
//    const double height =  130 * pt2mm * koef;// 166 ## mm (variable)
    double width = settings.Size[0];
    double height = settings.Size[1];

    const double mm2px = dpi / in2mm;//
    size_t pixels_width = std::round( width * mm2px);//
    size_t pixels_height = std::round( height * mm2px);//
    //------------------------------------------------------------------------------------------------------------------
//    std::map<std::string, std::string> true_keywords = { { "color", "black" }, { "linestyle", "--" },
//        { "linewidth", "1" }, { "label", "Истинное" } };
//    std::map<std::string, std::string> marks_keywords = { { "color", "grey" }, { "linestyle", "" }, { "marker", "." },
//        { "label", "Измерение" } };
/*
    std::map<std::string, std::string> graphColors = {
#ifdef EKF_
        { "EKF", "blue" },
#endif
#ifdef UKF_
        { "UKF", "darkgreen" },
#endif
#ifdef SRUKF_
        { "SRUKF", "limegreen" },
        { "SRUKFB", "lime" },
#endif
#ifdef CKF_
        { "CKF", "darkred" },
#endif
#ifdef SRCKF_
        { "SRCKF", "indianred" },
        { "SRCKFB", "red" },
#endif
#ifdef ECKF_
        { "ECKF", "purple" },
#endif
#ifdef SRECKF_
        { "SRECKF", "mediumpurple" },
        { "SRECKFB", "magenta" },
#endif
#ifdef EUKF_
        { "EUKF", "chocolate" },
#endif
#ifdef SREUKF_
        { "SREUKF", "darkorange" },
        { "SREUKFB", "orange" },
#endif
#ifdef SREKF_
        { "SREKF", "cornflowerblue" }
#endif
    };
    std::map<std::string, std::map<std::string, std::string>> estimated_keywords = {
#ifdef EKF_
        { "EKF", { { "color", graphColors.at("EKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "EKF" } } },
#endif
#ifdef UKF_
        { "UKF", { { "color", graphColors.at("UKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "UKF" } } },
#endif
#ifdef SRUKF_
        { "SRUKF", { { "color", graphColors.at("SRUKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRUKF" } } },
        { "SRUKFB", { { "color", graphColors.at("SRUKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRUKFB" } } },
#endif
#ifdef CKF_
        { "CKF", { { "color", graphColors.at("CKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "CKF" } } },
#endif
#ifdef SRCKF_
        { "SRCKF", { { "color", graphColors.at("SRCKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRCKF" } } },
        { "SRCKFB", { { "color", graphColors.at("SRCKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRCKFB" } } },
#endif
#ifdef ECKF_
        { "ECKF", { { "color", graphColors.at("ECKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "ECKF" } } },
#endif
#ifdef SRECKF_
        { "SRECKF", { { "color", graphColors.at("SRECKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRECKF" } } },
        { "SRECKFB", { { "color", graphColors.at("SRECKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRECKFB" } } },
#endif
#ifdef EUKF_
        { "EUKF", { { "color", graphColors.at("EUKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "EUKF" } } },
#endif
#ifdef SREUKF_
        { "SREUKF", { { "color", graphColors.at("SREUKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREUKF" } } },
        { "SREUKFB", { { "color", graphColors.at("SREUKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREUKFB" } } },
#endif
#ifdef SREKF_
        { "SREKF", { { "color", graphColors.at("SREKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREKF" } } }
#endif
    };
    std::map<std::string, std::map<std::string, std::string>> estimated_keywords_time = {
#ifdef EKF_
        { "EKF", { { "color", graphColors.at("EKF") }, { "label", "EKF" } } },
#endif
#ifdef UKF_
        { "UKF", { { "color", graphColors.at("UKF") }, { "label", "UKF" } } },
#endif
#ifdef SRUKF_
        { "SRUKF", { { "color", graphColors.at("SRUKF") }, { "label", "SRUKF" } } },
        { "SRUKFB", { { "color", graphColors.at("SRUKFB") }, { "label", "SRUKFB" } } },
#endif
#ifdef CKF_
        { "CKF", { { "color", graphColors.at("CKF") }, { "label", "CKF" } } },
#endif
#ifdef SRCKF_
        { "SRCKF", { { "color", graphColors.at("SRCKF") }, { "label", "SRCKF" } } },
        { "SRCKFB", { { "color", graphColors.at("SRCKFB") }, { "label", "SRCKFB" } } },
#endif
#ifdef ECKF_
        { "ECKF", { { "color", graphColors.at("ECKF") }, { "label", "ECKF" } } },
#endif
#ifdef SRECKF_
        { "SRECKF", { { "color", graphColors.at("SRECKF") }, { "label", "SRECKF" } } },
        { "SRECKFB", { { "color", graphColors.at("SRECKFB") }, { "label", "SRECKFB" } } },
#endif
#ifdef EUKF_
        { "EUKF", { { "color", graphColors.at("EUKF") }, { "label", "EUKF" } } },
#endif
#ifdef SREUKF_
        { "SREUKF", { { "color", graphColors.at("SREUKF") }, { "label", "SREUKF" } } },
        { "SREUKFB", { { "color", graphColors.at("SREUKFB") }, { "label", "SREUKFB" } } },
#endif
#ifdef SREKF_
        { "SREKF", { { "color", graphColors.at("SREKF") }, { "label", "SREKF" } } },
#endif
    };
*/
    std::string filters_names = "COMPARE_EKF";

    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        filters_names += "_SRUKF";
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        filters_names += "_SREUKF";
    }

    //std::string filters_names = "COMPARE_EKF_SRUKF";

//    for( int i = 0; i < settings.Filters.size(); i++ ) {
//        filters_names += settings.Filters[i];
//        if( i < settings.Filters.size() - 1 ) {
//            filters_names += "_";
//        }
//    }
    std::map<std::string, std::string>legend_loc { { "loc", "upper right" } };
    //std::map<std::string, std::string>legend_loc_time { { "loc", "upper left" } };

    //------------------------------------------------------------------------------------------------------------------
    // Буферы
    int N = static_cast<int>( settings.SimulationTime / settings.DeltaT ); // Количество моментов времени

    std::vector<double> time(N);
    arma::mat true_X( SizeX, N, arma::fill::zeros ); // Состояние Х (столбец - вектор в i-ый момент времени)
    arma::mat true_Y( SizeY, N, arma::fill::zeros );
    arma::mat measured_Y( SizeY, N, arma::fill::zeros );

    std::vector<double> template_vector_N(N);
    arma::mat template_mat_X( SizeX, N, arma::fill::zeros );
    arma::mat template_mat_Y( SizeY, N, arma::fill::zeros );

//    std::map<std::string, arma::mat> template_map_X;
//    std::map<std::string, arma::mat> template_map_Y;
//    std::map<std::string, std::vector<double>> template_map_N;

//    std::map<std::string, SPML::Timing::CTimeKeeper> timerPrediction;
//    std::map<std::string, SPML::Timing::CTimeKeeper> timerCorrection;
//    std::map<std::string, double> timerPredictionAverage;
//    std::map<std::string, double> timerCorrectionAverage;

//    for( auto &filter : settings.Filters ) {
//        template_map_X.insert( std::make_pair( filter, template_mat_X ) );
//        template_map_Y.insert( std::make_pair( filter, template_mat_Y ) );
//        template_map_N.insert( std::make_pair( filter, template_vector_N ) );
//        timerPrediction.insert( std::make_pair( filter, SPML::Timing::CTimeKeeper() ) );
//        timerCorrection.insert( std::make_pair( filter, SPML::Timing::CTimeKeeper() ) );
//        timerPredictionAverage.insert( std::make_pair( filter, 0.0 ) );
//        timerCorrectionAverage.insert( std::make_pair( filter, 0.0 ) );
//    }
//    std::map<std::string, arma::mat> estimated_X = template_map_X;
//    std::map<std::string, arma::mat> estimated_Y = template_map_Y;
//    std::map<std::string, arma::mat> estimated_Pdiag = template_map_X;
//    std::map<std::string, arma::mat> estimated_Sdiag = template_map_Y;
//    std::map<std::string, arma::mat> delta_Y = template_map_Y;
//    std::map<std::string, std::vector<double>> mahalanobis = template_map_N;
//    std::map<std::string, std::vector<double>> SDCM = template_map_N;

//    std::map<std::string, arma::mat> RMSE_X = template_map_X;
//    std::map<std::string, arma::mat> RMSE_EKF = template_map_X;

    arma::mat RMSE_X_EKF = template_mat_X;
    std::map<std::string, arma::mat> RMSE_X_SRUKF;
    std::map<std::string, arma::mat> RMSE_X_SREUKF;

//    std::map<std::string, std::vector<double>> timeMeasuredPrediction = template_map_N;
//    std::map<std::string, std::vector<double>> timeMeasuredCorrection = template_map_N;
//    std::map<std::string, std::vector<double>> timeMeasuredSumm = template_map_N;
    std::string name_tmp;

    std::mt19937 generator; // Генератор псевдослучайных чисел Mersenne Twister

    double RMS_X_X = 0.001; // км
    double RMS_X_Y = 0.001; // км
    double RMS_X_V = 0.1;  // м/c
    double RMS_X_K = 0.3; // град
    double RMS_X_Ka = 0.04; // град/с

    double resElementR = 0.30; // км
    double resElementAz = 0.05;//0.1; // град
    double resElementVr = 0.1; // м/с

    double coefK = 2.0;
    double weight_dB = 11.0; //11;//9;//
    double weight_times = std::pow( 10.0, ( weight_dB * 0.1 ) ); // Вес в разах
    double RMS_Y_R = coefK * resElementR / std::sqrt( 12.0 * weight_times );
    double RMS_Y_Az = coefK * resElementAz / std::sqrt( 12.0 * weight_times );
    double RMS_Y_Vr = coefK * resElementVr / std::sqrt( 12.0 * weight_times );

    std::normal_distribution<double> noiseY_R( 0.0, RMS_Y_R );
    std::normal_distribution<double> noiseY_Az( 0.0, RMS_Y_Az );
    std::normal_distribution<double> noiseY_Vr( 0.0, RMS_Y_Vr );

    std::uniform_real_distribution<double> random_0_1( 0.0, 1.0 ); // Вещественное случайное число от 0 до 1 с равномерной плотностью вероятности

    // Диагональ матрицы шумов состояния:
    arma::vec Q = {
        RMS_X_X * RMS_X_X,
        RMS_X_Y * RMS_X_Y,
        RMS_X_V * RMS_X_V,
        RMS_X_K * RMS_X_K,
        RMS_X_Ka * RMS_X_Ka
    };
    // Диагональ матрицы шумов измерений:
    arma::vec R = {
        RMS_Y_R * RMS_Y_R,
        RMS_Y_Az * RMS_Y_Az,
        RMS_Y_Vr * RMS_Y_Vr,
    };
    //------------------------------------------------------------------------------------------------------------------
    KalmanFilters::CKalmanEKF<SizeX, SizeY> EKF;
    EKF.SetStateTransitionModel( stateTransitionModel );
    EKF.SetObservationModel( observationModel );
    EKF.SetStateTransitionJacobianF( stateTransitionJacobianF );
    EKF.SetObservationJacobianH( observationJacobianH );

    EKF.SetCheckBordersStateAfterPrediction( checkBordersState );
    EKF.SetCheckBordersStateAfterCorrection( checkBordersState );
    EKF.SetCheckBordersMeasurement( checkBordersMeasurement );
    EKF.SetCheckDeltaState( checkDeltaState );
    EKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

    EKF.SetProcessCovarianceMatrixQdiag( Q );
    EKF.SetObservationCovarianceMatrixRdiag( R );
    //------------------------------------------------------------------------------------------------------------------
    std::map<std::string, KalmanFilters::CKalmanSRUKF<SizeX, SizeY>> filtersSRUKF;
    KalmanFilters::CKalmanSRUKF<SizeX, SizeY> SRUKF;
    SRUKF.SetStateTransitionModel( stateTransitionModel );
    SRUKF.SetObservationModel( observationModel );

    SRUKF.SetCheckBordersStateAfterPrediction( checkBordersState );
    SRUKF.SetCheckBordersStateAfterCorrection( checkBordersState );
    SRUKF.SetCheckBordersMeasurement( checkBordersMeasurement );
    SRUKF.SetCheckDeltaState( checkDeltaState );
    SRUKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

    SRUKF.SetWeightedSumStateSigmas( weightedSumStateSigmas );
    SRUKF.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

    SRUKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
    SRUKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    //------------------------------------------------------------------------------------------------------------------
    std::map<std::string, KalmanFilters::CKalmanSREUKF<SizeX, SizeY>> filtersSREUKF;
    KalmanFilters::CKalmanSREUKF<SizeX, SizeY> SREUKF;
    SREUKF.SetStateTransitionModel( stateTransitionModel );
    SREUKF.SetObservationModel( observationModel );
    SREUKF.SetStateTransitionJacobianF( stateTransitionJacobianF );
    SREUKF.SetObservationJacobianH( observationJacobianH );

    SREUKF.SetCheckBordersStateAfterPrediction( checkBordersState );
    SREUKF.SetCheckBordersStateAfterCorrection( checkBordersState );
    SREUKF.SetCheckBordersMeasurement( checkBordersMeasurement );
    SREUKF.SetCheckDeltaState( checkDeltaState );
    SREUKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

    SREUKF.SetWeightedSumStateSigmas( weightedSumStateSigmas );
    SREUKF.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

    SREUKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q ) );
    SREUKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    //------------------------------------------------------------------------------------------------------------------
    if( settings.Set == 0 ) { // Julier
        filters_names += "_Julier";
        for( unsigned long long i = 0; i < settings.w0_sr.size(); i++ ) {
            SRUKF.SetupDesignParametersMeanSet( settings.w0_sr[i] );
            SREUKF.SetupDesignParametersMeanSet( settings.w0_sr[i] );

            std::string key = "w0=" + std::to_string( settings.w0_sr[i] );
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                filtersSRUKF.insert( std::make_pair( key, SRUKF ) );
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                filtersSREUKF.insert( std::make_pair( key, SREUKF ) );
            }
        }
    } else if( settings.Set == 1 ) { // Merwe
        filters_names += "_Merwe";
        for( unsigned long long a = 0; a < settings.alpha_sr.size(); a++ )
        for( unsigned long long b = 0; b < settings.beta_sr.size(); b++ )
        for( unsigned long long k = 0; k < settings.kappa_sr.size(); k++ ) {
            SRUKF.SetupDesignParametersScaledSet( settings.alpha_sr[a], settings.beta_sr[b], settings.kappa_sr[k] );
            SREUKF.SetupDesignParametersScaledSet( settings.alpha_sr[a], settings.beta_sr[b], settings.kappa_sr[k] );

            std::string key = "a_b_k=" +
                std::to_string( settings.alpha_sr[a] ) + "_" +
                std::to_string( settings.beta_sr[b] ) + "_" +
                std::to_string( settings.kappa_sr[k] );
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                filtersSRUKF.insert( std::make_pair( key, SRUKF ) );
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                filtersSREUKF.insert( std::make_pair( key, SREUKF ) );
            }
        }
    } else {
        assert( false );
    }
    //------------------------------------------------------------------------------------------------------------------
    // Начало рабочих циклов
    //------------------------------------------------------------------------------------------------------------------

    std::string name_dt = "_dt_" + std::to_string( settings.DeltaT );

    int cycle_max = settings.MCruns * N * ( 1.0 + filtersSRUKF.size() + filtersSREUKF.size() );
    int cycle = 0;
    int prev_percent = -1;

    //------------------------------------------------------------------------------------------------------------------
    for( uint32_t seed_ = 1; seed_ <= settings.MCruns; seed_++ ) {
        generator.seed( seed_ ); // Выставить зерно ГСЧ
//            std::string name_seed_p = "_seed_" + std::to_string( seed_ ) + "_P_" + to_string_with_precision( p_, 1 ) + "_dt_" + std::to_string( settings.DeltaT );
        std::string name_seed_p = "_seed_" + std::to_string( seed_ ) + "_" + name_dt;
        //--------------------------------------------------------------------------------------------------------------
        // Начальное положение X, Y, Delta, P
        true_X.col(0) = arma::vec{ 100.0, 200.0, 100.0, 45.0, 1.0e-5 };
        true_Y.col(0) = observationModel( true_X.col(0) );
        measured_Y.col(0) = true_Y.col(0) + arma::vec{ noiseY_R( generator ), noiseY_Az( generator ), noiseY_Vr( generator ) };
        ///
        double v0 = std::abs( ( measured_Y.col(0) )(2) ); // Vr
        double NV = 2.0;// Более или равен 1 (эмпирически)
        double Vmax = 300.0; // м/с
        double startV = ( Vmax / NV ) + ( ( ( NV - 1.0 ) / ( NV * Vmax ) ) * v0 * v0 );
        int sign = 1;
        if( ( measured_Y.col(0) )(2) < 0 ) { // Vr < 0
            sign = -1;
        }
        double startK = ( measured_Y.col(0) )(1) * sign; // Az
        ///
        arma::vec startX {
            ( measured_Y.col(0) )(0) * std::sin( ( measured_Y.col(0) )(1) * DgToRd ),
            ( measured_Y.col(0) )(0) * std::cos( ( measured_Y.col(0) )(1) * DgToRd ),
            startV,
            startK,
            1.0e-5
        };
        arma::vec startY = observationModel( startX );
        arma::vec deltaY = measured_Y.col(0) - startY;
        deltaY = checkDeltaMeasurement( deltaY );

        double dispR = resElementR * resElementR / 12.0;
        double sigmaV = ( Vmax - ( measured_Y.col(0) )(2) ); // ( Vmax - measured_Y_Vr[0] )
        double sigmaK = 2.0 * std::acos( v0 / Vmax ) * RdToDg;
        double dispV = sigmaV * sigmaV / 12.0;
        double dispK = sigmaK * sigmaK / 12.0;

        arma::vec startP {
            dispR * 1.0,
            dispR * 1.0,
            dispV * 0.25,
            dispK * 0.25,
            ( RMS_X_Ka * RMS_X_Ka ) * 1.0
        };

        //--------------------------------------------------------------------------------------------------------------
        // EKF
        EKF.SetEstimatedVectorX( startX );
        EKF.SetEstimatedVectorY( startY );
        EKF.SetMeasuredVectorY( measured_Y.col(0) );
        EKF.SetDeltaY( deltaY );

        arma::mat PdenseEKF = arma::mat( SizeX, SizeX );
        PdenseEKF.fill( 1.0e-9 );
        PdenseEKF.diag() = startP;
        EKF.SetEstimateCovarianceMatrixP( PdenseEKF );
        if( settings.Debug ) {
            PdenseEKF.print("PdenseEKF:");
        }

//        RMSE_X.at("EKF").col(0) += arma::square( checkDeltaState( startX - true_X.col(0) ) );
        RMSE_X_EKF.col(0) += arma::square( checkDeltaState( startX - true_X.col(0) ) );
        //--------------------------------------------------------------------------------------------------------------
        // MAP OF SRUKF
        if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
            for( auto &item : filtersSRUKF ) {
                auto &key = item.first;
                auto &filter = item.second;

                filter.SetEstimatedVectorX( startX );
                filter.SetEstimatedVectorY( startY );
                filter.SetMeasuredVectorY( measured_Y.col(0) );
                filter.SetDeltaY( deltaY );

                arma::mat PdenseSRUKF = arma::mat( SizeX, SizeX );
                PdenseSRUKF.fill( 1.0e-9 );
                PdenseSRUKF.diag() = startP;
                PdenseSRUKF = arma::chol( PdenseSRUKF, "lower" );
                filter.SetEstimateCovarianceMatrixP( PdenseSRUKF );
                if( settings.Debug ) {
                    PdenseSRUKF.print("PdenseSRUKF:");
                }

                RMSE_X_SRUKF.insert( std::make_pair( key, template_mat_X ) );

                RMSE_X_SRUKF.at( key ).col(0) += arma::square( checkDeltaState( startX - true_X.col(0) ) );
            }
        }
        //--------------------------------------------------------------------------------------------------------------
        // MAP OF SREUKF
        if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
            for( auto &item : filtersSREUKF ) {
                auto &key = item.first;
                auto &filter = item.second;

                filter.SetEstimatedVectorX( startX );
                filter.SetEstimatedVectorY( startY );
                filter.SetMeasuredVectorY( measured_Y.col(0) );
                filter.SetDeltaY( deltaY );

                arma::mat PdenseSREUKF = arma::mat( SizeX, SizeX );
                PdenseSREUKF.fill( 1.0e-9 );
                PdenseSREUKF.diag() = startP;
                PdenseSREUKF = arma::chol( PdenseSREUKF, "lower" );
                filter.SetEstimateCovarianceMatrixP( PdenseSREUKF );
                if( settings.Debug ) {
                    PdenseSREUKF.print("PdenseSREUKF:");
                }

                RMSE_X_SREUKF.insert( std::make_pair( key, template_mat_X ) );

                RMSE_X_SREUKF.at( key ).col(0) += arma::square( checkDeltaState( startX - true_X.col(0) ) );
            }
        }
        //----------------------------------------------------------------------------------------------------------
        // Симуляция
        for( int i = 1; i < N; i++ ) { // Цикл по тактам времени
            time[i] = i * settings.DeltaT;
            if( settings.Debug ) {
                std::cout << "Takt = " << i << " Time = " << time[i] << std::endl;
            }

            // Прогноз
            EKF.Prediction( settings.DeltaT );

            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                for( auto &item : filtersSRUKF ) {
//                    auto &key = item.first;
                    auto &filter = item.second;
                    filter.Prediction( settings.DeltaT );
                }
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                for( auto &item : filtersSREUKF ) {
//                    auto &key = item.first;
                    auto &filter = item.second;
                    filter.Prediction( settings.DeltaT );
                }
            }

            // Создание измерений текущего такта:
            true_X.col(i) = stateTransitionModel( true_X.col(i - 1), settings.DeltaT );

            // Маневр:
            double start_K_man = 2000.0;
            double end_K_man = start_K_man + 60 * 5;
//                double start_K_man = 4000.0;
//                double end_K_man = start_K_man + 60 * 15;
            if( time[i] > start_K_man && time[i] < end_K_man ) // Если манёвр
            {
                ( true_X.col(i) )(4) = 0.6;//1.0;// // Ka
            }
            if( time[i] >= end_K_man ) { // Манёвр кончился
                ( true_X.col(i) )(4) = 1.0e-5;//0.0;//
            }
            // Measurement Y
            true_Y.col(i) = observationModel( true_X.col(i) );
            measured_Y.col(i) = true_Y.col(i) +
                arma::vec{ noiseY_R( generator ), noiseY_Az( generator ), noiseY_Vr( generator ) };

            // Коррекция
            EKF.Correction( measured_Y.col(i) );
            cycle++;
            print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения

            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                for( auto &item : filtersSRUKF ) {
//                    auto &key = item.first;
                    auto &filter = item.second;
                    filter.Correction( measured_Y.col(i) );
                    cycle++;
                    print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                }
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                for( auto &item : filtersSREUKF ) {
//                    auto &key = item.first;
                    auto &filter = item.second;
                    filter.Correction( measured_Y.col(i) );
                    cycle++;
                    print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                }
            }

            // RMSE
            RMSE_X_EKF.col(i) += arma::square( checkDeltaState( EKF.GetEstimatedVectorX() - true_X.col(i) ) );
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                for( auto &item : filtersSRUKF ) {
                    auto &key = item.first;
                    auto &filter = item.second;
                    RMSE_X_SRUKF.at( key ).col(i) += arma::square( checkDeltaState( filter.GetEstimatedVectorX() - true_X.col(i) ) );
                }
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                for( auto &item : filtersSREUKF ) {
                    auto &key = item.first;
                    auto &filter = item.second;
                    RMSE_X_SREUKF.at( key ).col(i) += arma::square( checkDeltaState( filter.GetEstimatedVectorX() - true_X.col(i) ) );
                }
            }
        } // end Цикл по тактам времени
    } // seed
    //--------------------------------------------------------------------------------------------------------------
    // RMSE
    int color_SRUKF_start = 70;//100;
    int color_SRUKF_end = 255;
    int color_SRUKF_step = 0;
    if( filtersSRUKF.size() > 0 ) {
        color_SRUKF_step = ( color_SRUKF_end - color_SRUKF_start ) / filtersSRUKF.size();
    }
    int color_SREUKF_start = 70;//100;
    int color_SREUKF_end = 255;
    int color_SREUKF_step = 0;
    if( filtersSREUKF.size() > 0 ) {
        color_SREUKF_step = ( color_SREUKF_end - color_SREUKF_start ) / filtersSREUKF.size();
    }

    RMSE_X_EKF *= ( 1.0 / static_cast<double>( settings.MCruns ) );
    RMSE_X_EKF = arma::sqrt( RMSE_X_EKF );

    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;
            RMSE_X_SRUKF.at( key ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
            RMSE_X_SRUKF.at( key ) = arma::sqrt( RMSE_X_SRUKF.at( key ) );
        }
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;
            RMSE_X_SREUKF.at( key ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
            RMSE_X_SREUKF.at( key ) = arma::sqrt( RMSE_X_SREUKF.at( key ) );
        }
    }
    if( settings.GraphSeparated ) {
        plt::figure_size( pixels_width, pixels_height );
    } else {
        plt::figure_size( pixels_width, pixels_height );
        plt::subplot( 3, 2, 1 );
    }
    plt::title( "RMSE координаты X" );
    plt::xlabel( "Время, [с]" );
    plt::ylabel( "RMSE X, [км]");
    int i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row(0) ),
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#" + stream.str() + "0000"; // RED

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row(0) ), //estimated_keywords.at( key ) );
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row(0) ), //estimated_keywords.at( filter ) );
        { { "color", "blue" }, { "linestyle", "-" }, { "linewidth", "1" } } );

//    plt::legend( legend_loc );
    plt::grid( true );
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 0.03 );
    }
    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_x_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }

    if( settings.GraphSeparated ) {
        plt::figure_size( pixels_width, pixels_height );
    } else {
        plt::subplot( 3, 2, 3 );
    }
    plt::title( "RMSE координаты Y");
    plt::xlabel( "Время, [с]" );
    plt::ylabel( "RMSE Y, [км]");
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row(1) ), //estimated_keywords.at( key ) );
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#" + stream.str() + "0000"; // RED

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row(1) ), //estimated_keywords.at( key ) );
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row(1) ), //estimated_keywords.at( filter ) );
        { { "color", "blue" }, { "linestyle", "-" }, { "linewidth", "1" } } );
//    plt::legend( legend_loc );
    plt::grid( true );
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 0.03 );
    }
    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_y_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }

    if( settings.GraphSeparated ) {
        plt::figure_size( pixels_width, pixels_height );
    } else {
        plt::subplot( 3, 2, 2 );
    }
    plt::title( "RMSE полной скорости");
    plt::xlabel( "Время, [с]" );
    plt::ylabel( "RMSE V, [м/с]");
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row(2) ), //estimated_keywords.at( key ) );
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#" + stream.str() + "0000"; // RED

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row(2) ), //estimated_keywords.at( key ) );
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row(2) ), //estimated_keywords.at( filter ) );
        { { "color", "blue" }, { "linestyle", "-" }, { "linewidth", "1" } } );
//    plt::legend( legend_loc );
    plt::grid( true );
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 10.0 ); // 15
    }
    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_v_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }

    if( settings.GraphSeparated ) {
        plt::figure_size( pixels_width, pixels_height );
    } else {
        plt::subplot( 3, 2, 4 );
    }
    plt::title( "RMSE курса");
    plt::xlabel( "Время, [с]" );
    plt::ylabel( "RMSE К, [градусы]");
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row(3) ), //estimated_keywords.at( key ) );
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#" + stream.str() + "0000"; // RED

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row(3) ), //estimated_keywords.at( key ) );
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row(3) ), //estimated_keywords.at( filter ) );
        { { "color", "blue" }, { "linestyle", "-" }, { "linewidth", "1" } } );
//    plt::legend( legend_loc );
    plt::grid( true );
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 25.0 );
    }
    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_k_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }

    if( settings.GraphSeparated ) {
        plt::figure_size( pixels_width, pixels_height );
    } else {
        plt::subplot( 3, 2, 6 );
    }
    plt::title( "RMSE скорости изменения курса" );
    plt::xlabel( "Время, [с]" );
    plt::ylabel( "RMSE dK/dt, [градусы/с]");
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row(4) ), //estimated_keywords.at( key ) );
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
//            auto &filter = item.second;

            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#" + stream.str() + "0000"; // RED

            plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row(4) ), //estimated_keywords.at( key ) );
                { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row(4) ), //estimated_keywords.at( filter ) );
        { { "color", "blue" }, { "linestyle", "-" }, { "linewidth", "1" } } );
//    plt::legend( legend_loc );
    plt::grid( true) ;
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 0.18 );
    }
    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_ka_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }

    if( !settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }

    if( settings.ShowGraphs ) {
        plt::show();
    }
    plt::close();
    //------------------------------------------------------------------------------------------------------------------
    plt::clf();
    plt::cla();
    Py_Finalize();

    }
    catch( std::exception &ex )
    {
        std::cout << "Exception occured: " << ex.what() << std::endl;
    }
    std::cout << "Завершено." << std::endl;
}

/// \}

