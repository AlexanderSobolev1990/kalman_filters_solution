//------------------------------------------------------------------------------
///
/// \file       compare_run_main.cpp
/// \brief      Сравнение фильтров Калмана
/// \date       26.12.24 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#include <compare_run.h>

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

    keywords.insert( std::make_pair( "left", settings.MatPlotParams[0] ) );
    keywords.insert( std::make_pair( "bottom", settings.MatPlotParams[1] ) );
    keywords.insert( std::make_pair( "right", settings.MatPlotParams[2] ) );
    keywords.insert( std::make_pair( "top", settings.MatPlotParams[3] ) );
    keywords.insert( std::make_pair( "wspace", settings.MatPlotParams[4] ) );
    keywords.insert( std::make_pair( "hspace", settings.MatPlotParams[5] ) );

    const double in2mm = 25.4;// mm (fixed)

    const double dpi = 300;// dpi (variable)
    double width = settings.Size[0];
    double height = settings.Size[1];

    const double mm2px = dpi / in2mm;//
    size_t pixels_width = std::round( width * mm2px);//
    size_t pixels_height = std::round( height * mm2px);//
    //--------------------------------------------------------------------------
    std::map<std::string, std::string> true_start = { { "color", "black" }, { "linestyle", "" }, { "marker", "." },
        { "markersize", "12" }, { "label", "Истинное начальное положение" } };
    std::map<std::string, std::string> estimated_start = { { "color", "black" }, { "marker", "o" },
        { "label", "Оцененное начальное положение" } };
    std::map<std::string, std::string> true_keywords = { { "color", "black" }, { "linestyle", "--" }, //{ "marker", "." },
        { "linewidth", "1" }, { "label", "Истинная траектория" } };
    std::map<std::string, std::string> marks_keywords = { { "color", "grey" }, { "linestyle", "" }, { "marker", "." },
        { "label", "Измерение" } };

    std::map<std::string, std::string> graphColors = {        
        { "EKF", "darkblue" }, 
        { "SREKF", "cornflowerblue" },
        { "SREKFB", "blue" },
        { "UKF", "green" },
        { "SRUKF", "limegreen" },
        { "SRUKFB", "lime" },
        { "CKF", "darkred" },
        { "SRCKF", "indianred" },
        { "SRCKFB", "red" },        
        { "EUKF", "chocolate" },
        { "SREUKF", "orange" },
        { "SREUKFB", "gold" },
        { "ECKF", "purple" },
        { "SRECKF", "mediumpurple" },
        { "SRECKFB", "magenta" }
    };

    std::map<std::string, std::map<std::string, std::string>> estimated_keywords = {
        { "EKF", { { "color", graphColors.at("EKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "EKF" } } },
        { "SREKF", { { "color", graphColors.at("SREKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREKF" } } },
        { "SREKFB", { { "color", graphColors.at("SREKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREKFB" } } },
        { "UKF", { { "color", graphColors.at("UKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "UKF" } } },
        { "SRUKF", { { "color", graphColors.at("SRUKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRUKF" } } },
        { "SRUKFB", { { "color", graphColors.at("SRUKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRUKFB" } } },
        { "CKF", { { "color", graphColors.at("CKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "CKF" } } },
        { "SRCKF", { { "color", graphColors.at("SRCKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRCKF" } } },
        { "SRCKFB", { { "color", graphColors.at("SRCKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRCKFB" } } },        
        { "EUKF", { { "color", graphColors.at("EUKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "EUKF" } } },
        { "SREUKF", { { "color", graphColors.at("SREUKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREUKF" } } },
        { "SREUKFB", { { "color", graphColors.at("SREUKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREUKFB" } } },
        { "ECKF", { { "color", graphColors.at("ECKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "ECKF" } } },
        { "SRECKF", { { "color", graphColors.at("SRECKF") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRECKF" } } },
        { "SRECKFB", { { "color", graphColors.at("SRECKFB") }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SRECKFB" } } },        
    };
    
    std::map<std::string, std::map<std::string, std::string>> estimated_keywords_time = {
        { "EKF", { { "color", graphColors.at("EKF") }, { "label", "EKF" } } },
        { "SREKF", { { "color", graphColors.at("SREKF") }, { "label", "SREKF" } } },
        { "SREKFB", { { "color", graphColors.at("SREKFB") }, { "label", "SREKFB" } } },
        { "UKF", { { "color", graphColors.at("UKF") }, { "label", "UKF" } } },
        { "SRUKF", { { "color", graphColors.at("SRUKF") }, { "label", "SRUKF" } } },
        { "SRUKFB", { { "color", graphColors.at("SRUKFB") }, { "label", "SRUKFB" } } },
        { "CKF", { { "color", graphColors.at("CKF") }, { "label", "CKF" } } },
        { "SRCKF", { { "color", graphColors.at("SRCKF") }, { "label", "SRCKF" } } },
        { "SRCKFB", { { "color", graphColors.at("SRCKFB") }, { "label", "SRCKFB" } } },        
        { "EUKF", { { "color", graphColors.at("EUKF") }, { "label", "EUKF" } } },
        { "SREUKF", { { "color", graphColors.at("SREUKF") }, { "label", "SREUKF" } } },
        { "SREUKFB", { { "color", graphColors.at("SREUKFB") }, { "label", "SREUKFB" } } },
        { "ECKF", { { "color", graphColors.at("ECKF") }, { "label", "ECKF" } } },
        { "SRECKF", { { "color", graphColors.at("SRECKF") }, { "label", "SRECKF" } } },
        { "SRECKFB", { { "color", graphColors.at("SRECKFB") }, { "label", "SRECKFB" } } }
    };

    std::string filters_names;
    for( std::size_t i = 0; i < settings.Filters.size(); i++ ) {
        filters_names += settings.Filters[i];
        if( i < settings.Filters.size() - 1 ) {
            filters_names += "_";
        }
    }
    std::map<std::string, std::string>legend_loc{ { "loc", "upper right" } };
    std::map<std::string, std::string>legend_loc_time { { "loc", "upper left" } };
    std::string loc_time = "upper left";
    std::string loc = "center right";

    //--------------------------------------------------------------------------
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
    std::map<std::string, arma::mat> template_map_mu;

    std::map< std::string, KalmanFilters::CKalmanBase<SizeX, SizeY>* > filters;

    for( auto &filter : settings.Filters ) {
        template_map_X.insert( std::make_pair( filter, template_mat_X ) );
        template_map_Y.insert( std::make_pair( filter, template_mat_Y ) );
        template_map_N.insert( std::make_pair( filter, template_vector_N ) );
        timerPrediction.insert( std::make_pair( filter, SPML::Timing::CTimeKeeper() ) );
        timerCorrection.insert( std::make_pair( filter, SPML::Timing::CTimeKeeper() ) );
        timerPredictionAverage.insert( std::make_pair( filter, 0.0 ) );
        timerCorrectionAverage.insert( std::make_pair( filter, 0.0 ) );
        
        if( filter == "EKF" ) {
            filters.insert( std::make_pair( "EKF", new KalmanFilters::CKalmanEKF<SizeX, SizeY>() ) );
        } else if( filter == "SREKF" ) {
            filters.insert( std::make_pair( "SREKF", new KalmanFilters::CKalmanSREKF<SizeX, SizeY>() ) );
        } else if( filter == "SREKFB" ) {
            filters.insert( std::make_pair( "SREKFB", new KalmanFilters::CKalmanSREKFB<SizeX, SizeY>() ) );
        } else if( filter == "UKF" ) {
            filters.insert( std::make_pair( "UKF", new KalmanFilters::CKalmanUKF<SizeX, SizeY>() ) );
        } else if( filter == "SRUKF" ) {
            filters.insert( std::make_pair( "SRUKF", new KalmanFilters::CKalmanSRUKF<SizeX, SizeY>() ) );
        } else if( filter == "SRUKFB" ) {
            filters.insert( std::make_pair( "SRUKFB", new KalmanFilters::CKalmanSRUKFB<SizeX, SizeY>() ) );
        } else if( filter == "CKF" ) {
            filters.insert( std::make_pair( "CKF", new KalmanFilters::CKalmanCKF<SizeX, SizeY>() ) );
        } else if( filter == "SRCKF" ) {
            filters.insert( std::make_pair( "SRCKF", new KalmanFilters::CKalmanSRCKF<SizeX, SizeY>() ) );
        } else if( filter == "SRCKFB" ) {
            filters.insert( std::make_pair( "SRCKFB", new KalmanFilters::CKalmanSRCKFB<SizeX, SizeY>() ) );
        } else if( filter == "EUKF" ) {
            filters.insert( std::make_pair( "EUKF", new KalmanFilters::CKalmanEUKF<SizeX, SizeY>() ) );
        } else if( filter == "SREUKF" ) {
            filters.insert( std::make_pair( "SREUKF", new KalmanFilters::CKalmanSREUKF<SizeX, SizeY>() ) );
        } else if( filter == "SREUKFB" ) {
            filters.insert( std::make_pair( "SREUKFB", new KalmanFilters::CKalmanSREUKFB<SizeX, SizeY>() ) );
        } else if( filter == "ECKF" ) {
            filters.insert( std::make_pair( "ECKF", new KalmanFilters::CKalmanECKF<SizeX, SizeY>() ) );        
        } else if( filter == "SRECKF" ) {
            filters.insert( std::make_pair( "SRECKF", new KalmanFilters::CKalmanSRECKF<SizeX, SizeY>() ) );
        } else if( filter == "SRECKFB" ) {
            filters.insert( std::make_pair( "SRECKFB", new KalmanFilters::CKalmanSRECKFB<SizeX, SizeY>() ) );
        }
    }
    std::map<std::string, arma::mat> estimated_X = template_map_X;
    std::map<std::string, arma::mat> estimated_Y = template_map_Y;
    std::map<std::string, arma::mat> estimated_Pdiag = template_map_X;
    std::map<std::string, arma::mat> estimated_Sdiag = template_map_Y;
    std::map<std::string, arma::mat> delta_Y = template_map_Y;
    std::map<std::string, std::vector<double>> mahalanobis = template_map_N;
    std::map<std::string, std::vector<double>> SDCM = template_map_N;
    std::map<std::string, arma::mat> RMSE_X = template_map_X;
    std::map<std::string, arma::mat> mean_Pdiag = template_map_X;
    std::map<std::string, arma::mat> estimated_mu = template_map_mu;

    std::map<std::string, std::vector<double>> timeMeasuredPrediction = template_map_N;
    std::map<std::string, std::vector<double>> timeMeasuredCorrection = template_map_N;
    std::map<std::string, std::vector<double>> timeMeasuredSumm = template_map_N;    
    std::string name_tmp;

    std::random_device rd;
    std::mt19937 generator; // Генератор псевдослучайных чисел Mersenne Twister

    std::normal_distribution<double> noise_Rm_Y_R( 0.0, RMS_Y_R );
    std::normal_distribution<double> noise_Rm_Y_Az( 0.0, RMS_Y_Az );
    std::normal_distribution<double> noise_Rm_Y_Vr( 0.0, RMS_Y_Vr );
    
    std::normal_distribution<double> noise_Rg_Y_R( 0.0, glintNoiseCoef * RMS_Y_R );
    std::normal_distribution<double> noise_Rg_Y_Az( 0.0, glintNoiseCoef * RMS_Y_Az );
    std::normal_distribution<double> noise_Rg_Y_Vr( 0.0, glintNoiseCoef * RMS_Y_Vr );

    std::uniform_real_distribution<double> random_0_1( 0.0, 1.0 ); // Вещественное случайное число от 0 до 1 с равномерной плотностью вероятности
        
    omp_set_num_threads( NUM_OF_THREADS );

    // #pragma omp parallel for
    // for( auto &filterName : settings.Filters  ) {        
    //    auto filter = filters.at( filterName );
    for( int f = 0; f < settings.Filters.size(); f++ ) {            
        auto filter = filters.at( settings.Filters[f] );
        filter->SetStateTransitionModel( stateTransitionModel );
        filter->SetObservationModel( observationModel );
        filter->SetStateTransitionJacobianF( stateTransitionJacobianF );
        filter->SetObservationJacobianH( observationJacobianH );

        filter->SetCheckBordersStateAfterPrediction( checkBordersState );
        filter->SetCheckBordersStateAfterCorrection( checkBordersState );
        filter->SetCheckBordersMeasurement( checkBordersMeasurement );
        filter->SetCheckDeltaState( checkDeltaState );
        filter->SetCheckDeltaMeasurement( checkDeltaMeasurement );

        // Установить матрицы Q и R
        arma::vec Q;
        switch( settings.eksperiment ) {
            case 0: {
                Q = Q_CV;                    
                break;
            }
//            case 1:
//            case 2:
//            case 3: {
//                Q = Q_CT;
//                break;
//            }
//            case 4: {
//                Q = Q_CA;
//                break;
//            }
            default:
                assert( false );
        }

        if( filter->GetCovType() == KalmanFilters::TCovType::CT_FullMatrices ) {
            if( filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Unscented ||
                filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Cubature ) 
            {
                filter->SetProcessCovarianceMatrixQdiag( Q % arma::vec( settings.q_koef_ukf ) ); // Поэлементное умножение
            } else {
                filter->SetProcessCovarianceMatrixQdiag( Q % arma::vec( settings.q_koef_ekf ) ); // Поэлементное умножение
            }
            filter->SetObservationCovarianceMatrixRdiag( R );
        } else if( filter->GetCovType() == KalmanFilters::TCovType::CT_SquareRootMatrices ) {                       
            if( filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Unscented ||
                filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Cubature ) 
            {
                filter->SetProcessCovarianceMatrixQdiag( arma::sqrt( Q % arma::vec( settings.q_koef_ukf ) ) ); // Поэлементное умножение
            } else {
                filter->SetProcessCovarianceMatrixQdiag( arma::sqrt( Q % arma::vec( settings.q_koef_ekf ) ) ); // Поэлементное умножение
            }
            filter->SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
        } else {
            assert( false );
        }

        // Функции взвешивания для сигма-точечных фильтров
        if( filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Unscented ||
            filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Cubature ) 
        {
            filter->SetWeightedSumStateSigmas( weightedSumStateSigmas );
            filter->SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );
        }

        if( filter->GetCovType() == KalmanFilters::TCovType::CT_FullMatrices ) {
            if( filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Unscented ) {
                if( settings.Set == 0 ) { // Julier            
                    filter->SetDesignParametersMeanSet( settings.w0[0] );
                } else if( settings.Set == 1 ) { // Merwe            
                    filter->SetDesignParametersScaledSet( settings.alpha[0], settings.beta[0], settings.kappa[0] );
                } else {
                    assert( false );
                }
            } else if( filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Cubature ) {
                filter->SetDesignParametersCubatureBaseSet();
            }
        } else if( filter->GetCovType() == KalmanFilters::TCovType::CT_SquareRootMatrices ) {
            if( filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Unscented ) {
                if( settings.Set == 0 ) { // Julier            
                    filter->SetDesignParametersMeanSet( settings.w0_sr[0] );
                } else if( settings.Set == 1 ) { // Merwe            
                    filter->SetDesignParametersScaledSet( settings.alpha_sr[0], settings.beta_sr[0], settings.kappa_sr[0] );
                } else {
                    assert( false );
                }
            } else if( filter->GetKalmanType() == KalmanFilters::TKalmanType::KT_Cubature ) {
                filter->SetDesignParametersCubatureBaseSet();
            }
        } else {
            assert( false );                
        }
    }
    //--------------------------------------------------------------------------
    // Начало рабочих циклов
    //--------------------------------------------------------------------------

    int cycle_max = settings.MCruns * N * settings.Filters.size();
    int cycle = 0;
    int prev_percent = -1;
    int Nfixed = N;// 0; // Укороченное время реализации!

    for( auto &p_ : settings.Probabilities ) {
        std::string name_p = "_P_" + to_string_with_precision( p_, 1 ) + "_dt_" + std::to_string( settings.DeltaT );
        
        SPML::Timing::CTimeKeeper TotalWorkTime;
        TotalWorkTime.StartTimer();
        
        //----------------------------------------------------------------------
        for( uint32_t seed_ = settings.MCseed; seed_ < ( settings.MCseed + settings.MCruns ); seed_++ ) {
            if( ( seed_ == 0 ) && ( settings.MCruns == 1 ) ) {
                auto rand_seed = rd();
                generator.seed( rand_seed ); // Выставить зерно ГСЧ случайным, если seed = 0
                std::cout << "random seed = " << rand_seed << std::endl;
            } else {
                generator.seed( seed_ ); // Выставить зерно ГСЧ
            }
            std::string name_seed_p = "_seed_" + std::to_string( seed_ ) + "_" + name_p;
            //------------------------------------------------------------------
            // Начальное положение X, Y, Delta, P
            true_X.col(0) = arma::vec( settings.x_start );
            true_Y.col(0) = observationModel( true_X.col(0) );
            
            double U = random_0_1( generator ); // Случайное число от 0 до 1
            if( ( U > Pg ) || ( Pg == 0.0 ) ) {
                measured_Y.col(0) = true_Y.col(0) + arma::vec{
                    noise_Rm_Y_R( generator ),
                    noise_Rm_Y_Az( generator ),
                    noise_Rm_Y_Vr( generator )            
                };
            } else {
                measured_Y.col(0) = true_Y.col(0) + arma::vec{
                    noise_Rg_Y_R( generator ),
                    noise_Rg_Y_Az( generator ),
                    noise_Rg_Y_Vr( generator )            
                };
            }

            ///
            double r = ( measured_Y.col(0) )(0);
            double az = ( measured_Y.col(0) )(1);
            double vf = ( measured_Y.col(0) )(2);
            double v0 = std::abs( vf ); // Vr
            double NV = 2.0;// Более или равен 1 (эмпирически)
            double Vmax = 300.0; // м/с
            double startV = ( Vmax / NV ) + ( ( ( NV - 1.0 ) / ( NV * Vmax ) ) * v0 * v0 );
            int sign = 1;
            if( vf < 0 ) {
                sign = -1;
            }
            double adding = std::acos( vf / startV ) * SPML::Convert::RdToDg;
            double startK = ( az * sign ) + adding;
            //
            // bool sideVariant = false;
            bool sideVariant = true;
            if( !sideVariant ) {
                startV = v0; // --> доплеровская
                startK = az; // --> курс равен азимуту
            }
            //
            arma::vec startX {
                r * std::sin( az * SPML::Convert::DgToRd ),
                r * std::cos( az * SPML::Convert::DgToRd ),
                startV,//
                startK,//                
                1.0e-5
            };
            arma::vec startY = observationModel( startX );
            arma::vec deltaY = measured_Y.col(0) - startY;
            deltaY = checkDeltaMeasurement( deltaY );

            double dispR = resElementR * resElementR / 12.0;
            double sigmaV = ( Vmax - v0 );
            double sigmaK = std::acos( v0 / Vmax ) * SPML::Convert::RdToDg;
            double dispV = sigmaV * sigmaV / 12.0;
            double dispK = sigmaK * sigmaK / 12.0;

            arma::mat H = observationJacobianH( startX );

            arma::vec startP {
                dispR * 1.0,
                dispR * 1.0,
                dispV * 0.25,//1.0,//
                dispK * 0.25,//1.0,//
                ( RMS_X_Ka * RMS_X_Ka ) * 1.0e4
            };
            arma::vec startP_UKF = arma::vec( startP );
            
            // Инициализация фильтров            
            // #pragma omp parallel for           
            // for( auto &filterName : settings.Filters  ) {        
            //     auto filter = filters.at( filterName );                
            for( int f = 0; f < settings.Filters.size(); f++ ) {            
                auto filter = filters.at( settings.Filters[f] );

                filter->SetEstimatedVectorX( startX );
                filter->SetEstimatedVectorY( startY );
                filter->SetMeasuredVectorY( measured_Y.col(0) );
                filter->SetDeltaY( deltaY );

                arma::mat Pdense = arma::mat( SizeX, SizeX );
                Pdense.fill( 1.0e-9 );
                Pdense.diag() = startP;
                if( filter->GetCovType() == KalmanFilters::TCovType::CT_SquareRootMatrices ) {                    
                    Pdense =  arma::chol( Pdense, "lower" );
                }
                filter->SetEstimateCovarianceMatrixP( Pdense );

                if( settings.Debug ) {
                    Pdense.print("Pdense" + filter->GetFilterName() );
                }

                estimated_X.at( filter->GetFilterName() ).col(0) = startX;
                estimated_Y.at( filter->GetFilterName() ).col(0) = startY;
                RMSE_X.at( filter->GetFilterName() ).col(0) += 
                    arma::square( checkDeltaState( estimated_X.at( filter->GetFilterName() ).col(0) - true_X.col(0) ) );
                if( filter->GetCovType() == KalmanFilters::TCovType::CT_FullMatrices ) {
                    estimated_Pdiag.at( filter->GetFilterName() ).col(0) = arma::sqrt( Pdense.diag() );
                    mean_Pdiag.at( filter->GetFilterName() ).col(0) += arma::sqrt( Pdense.diag() );
                    SDCM.at( filter->GetFilterName() )[0] = std::sqrt( arma::trace( Pdense ) );
                } else if( filter->GetCovType() == KalmanFilters::TCovType::CT_SquareRootMatrices ) {
                    estimated_Pdiag.at( filter->GetFilterName() ).col(0) = Pdense.diag();
                    mean_Pdiag.at( filter->GetFilterName() ).col(0) += Pdense.diag();
                    SDCM.at( filter->GetFilterName() )[0] = arma::trace( Pdense );
                }
            }
            //------------------------------------------------------------------
            // Симуляция
            for( int i = 1; i < N; i++ ) { // Цикл по тактам времени
                time[i] = i * settings.DeltaT;

                // Прогноз
                {                
                // #pragma omp parallel for               
                // for( auto &filterName : settings.Filters  ) {        
                //    auto filter = filters.at( filterName );
                for( int f = 0; f < settings.Filters.size(); f++ ) {
                    auto filter = filters.at( settings.Filters[f] );

                    timerPrediction.at( filter->GetFilterName() ).StartTimer();
                    filter->Prediction( settings.DeltaT );      
                    timerPrediction.at( filter->GetFilterName() ).EndTimer();
                    timeMeasuredPrediction.at( filter->GetFilterName() )[i] = 
                        timerPrediction.at( filter->GetFilterName() ).TimeCur() * 1.0e6; // мкс
                } 

                }               
                //------------------------------------------------------------------------------------------------------
                // Зададим маневры
                //------------------------------------------------------------------------------------------------------
                double Radius = 20 * 1e3; // [м] - радиус 20 км типичен дла АВАКС-а
                double AngleOfTurn = 0; // Угол поворота
                double StraightPath = 0; // Прямой участок, [м]

                if( settings.eksperiment == 0 ) { // Прямой полет
                    // Do nothing
                } else if( settings.eksperiment == 1 ) { // Прямой полет, правый поворот на 360
                    AngleOfTurn = SPML::Convert::DgToRd * 360;
                    StraightPath = Radius * 4;
                } else if( settings.eksperiment == 2 ) { // Прямой полет, правый поворот на 180, прямой полет, снова правый поворот на 180
                    AngleOfTurn = SPML::Convert::DgToRd * 180;
                    StraightPath = Radius * 4;
                } else if( settings.eksperiment == 3 ) { // Прямой полет, правый поворот на 270, прямой полет, левый поворот на 270, прямой полет
                    AngleOfTurn = SPML::Convert::DgToRd * 270;
                    StraightPath = 40e3;//Radius * 2;
                } else if( settings.eksperiment == 4 ) {
                    // do nothing
                } else if( settings.eksperiment == 5 ) {
                    AngleOfTurn = SPML::Convert::DgToRd * 360;
                    StraightPath = Radius * 4;
                } else if( settings.eksperiment == 6 ) {
                    AngleOfTurn = 1.e-6;
                    StraightPath = Radius * 2;
                } else {
                    assert( false );
                }
                double AccelerationTime = 100.0; // сек
                double StraightPathTime = StraightPath / settings.x_start[IndX::V];
                double TimeOfTurn = AngleOfTurn * Radius / settings.x_start[IndX::V];
                double dKdT = ( AngleOfTurn * SPML::Convert::RdToDg ) / TimeOfTurn; // градус / сек

                int Time1 = StraightPathTime;
                int Time1a = StraightPathTime / 3.0;
                int Time1b = 2.0 * StraightPathTime / 3.0;
                int Time2 = Time1 + TimeOfTurn;
                int Time3 = Time2 + StraightPathTime;
                int Time4 = Time3 + TimeOfTurn;
//                int Time5 = Time4 + StraightPathTime;

                double X_Ka = 0; // Задаём dK/dt
                double X_a = 0;

                if( settings.eksperiment == 0 ) { // Прямой полет
                    // do nothing!
                } else if( settings.eksperiment == 1 ) { // Прямой полет, правый поворот на 360
                    if( time[i] <= Time1 ) {
                        X_Ka = 0.0;
                    } else if( time[i] > Time1 && time[i] <= Time2 ) {
                        X_Ka = dKdT;
                    } else {
                        Nfixed = i;
                        break; // Прерывание цикла по тактам времени!
                    }
                } else if( settings.eksperiment == 2 ) { // Прямой полет, правый поворот на 180, прямой полет, снова правый поворот на 180
                    if( time[i] <= Time1 ) {
                        X_Ka = 0.0;
                    } else if( time[i] > Time1 && time[i] <= Time2 ) {
                        X_Ka = dKdT;
                    } else if( time[i] > Time2 && time[i] <= Time3 ) {
                        X_Ka = 0.0;
                    } else if( time[i] > Time3 && time[i] <= Time4 ) {
                        X_Ka = dKdT;
                    } else {
                        Nfixed = i;
                        break; // Прерывание цикла по тактам времени!
                    }
                } else if( settings.eksperiment == 3 ) { // Прямой полет, правый поворот на 270, прямой полет, левый поворот на 270, прямой полет
                    if( time[i] <= Time1 ) {
                        X_Ka = 0.0;
                    } else if( time[i] > Time1 && time[i] <= Time2 ) {
                        X_Ka = dKdT;
                    } else if( time[i] > Time2 && time[i] <= Time3 ) {
                        X_Ka = 0.0;
                    } else if( time[i] > Time3 && time[i] <= Time4 ) {
                        X_Ka = -dKdT;
                    } else {
                        Nfixed = i;
                        break; // Прерывание цикла по тактам времени!
                    }
                // } else if( settings.eksperiment == 4 ) {
                //     if( time[i] <= AccelerationTime ) {
                //         X_a = settings.x_start[IndX::a]; // Какое задано в начальном векторе состояния
                //     } else if( time[i] > AccelerationTime ) {
                //         X_a = 0;
                //     } else {
                //         Nfixed = i;
                //         break; // Прерывание цикла по тактам времени!
                //     }
                } else if( settings.eksperiment == 5 ) {
                    if( time[i] <= Time1a ) {
                        X_Ka = 0;
                        X_a = 0;
                    } else if( time[i] > Time1a && time[i] <= Time1b ) {
                        X_a = 1;//2;//
                        X_Ka = 0;
                    } else if( time[i] > Time1b && time[i] <= Time1 ) {
                        X_a = 0;
                        X_Ka = 0;
                    } else if( time[i] > Time1 && time[i] <= Time2 ) {
                        X_a = 0;
                        X_Ka = dKdT;
                    } else if( time[i] > Time2 && time[i] <= Time2 * 1.15 ) {
                        X_a = 0;
                        X_Ka = 0;
                    } else if( time[i] > Time2 * 1.15 && time[i] <= Time2 * 1.35 ) {
                        X_a = 1;//2;//
                        X_Ka = 0;
                    } else if( time[i] > Time2 * 1.35 && time[i] <= Time2 * 1.5 ) {
                        X_a = 0;
                        X_Ka = 0;
                    } else if( time[i] > Time2 * 1.5 && time[i] <= Time2 * 4 ) {
                        X_a = 0;
                        X_Ka = dKdT;
                    } else if( time[i] > Time2 * 4 && time[i] <= Time2 * 5 ) {
                        X_a = 0;
                        X_Ka = 0;
                    } else if( time[i] > Time2 * 5 && time[i] <= Time2 * 6 ) {
                        X_a = 1;//3;//
                        X_Ka = 0;
                    } else {
                        Nfixed = i;
                        break; // Прерывание цикла по тактам времени!
                    }
                } else if( settings.eksperiment == 6 ) {
                    if( time[i] <= Time1 ) {
                        X_Ka = 0;
                        X_a = 0;
                    } else if( time[i] > Time1 && time[i] <= ( Time1 + 60 ) ) {
                        X_a = 2;
                        X_Ka = 0;
                    } else if( time[i] > ( Time1 + 60 ) && time[i] <= ( Time1 + 120 ) ) {
                        X_a = 0;
                        X_Ka = 0;
                    } else if( time[i] > ( Time1 + 120 ) && time[i] <= ( Time1 + 180 ) ) {
                        X_a = -2;
                        X_Ka = 0;
                    } else if( time[i] > ( Time1 + 180 ) && time[i] <= ( Time1 + 240 ) ) {
                        X_a = 0;
                        X_Ka = 0;
                    } else if( time[i] > ( Time1 + 240 ) && time[i] <= ( Time1 + 300 ) ) {
                        X_a = 5;
                        X_Ka = 0;
                    } else if( time[i] > ( Time1 + 300 ) && time[i] <= ( Time1 + 360 ) ) {
                        X_a = 0;
                        X_Ka = 0;
                    } else if( time[i] > ( Time1 + 360 ) && time[i] <= ( Time1 + 420 ) ) {
                        X_a = -5;
                        X_Ka = 0;
                    } else if( time[i] > ( Time1 + 420 ) ) {
                        X_a = 0;
                        X_Ka = 0;
                    } else {
                        Nfixed = i;
                        break; // Прерывание цикла по тактам времени!
                    }
                } else {
                    assert( false );
                }

                arma::vec true_X_prev = true_X.col(i - 1);
                true_X_prev = arma::vec{
                    true_X_prev[IndX::X],
                    true_X_prev[IndX::Y],
                    true_X_prev[IndX::V],
                    true_X_prev[IndX::K],
                    true_X_prev[IndX::Ka]
                };
                arma::vec true_X_next = stateTransitionModel( true_X_prev, settings.DeltaT );
                true_X.col(i) = arma::vec{
                    true_X_next[IndX::X],
                    true_X_next[IndX::Y],
                    true_X_next[IndX::V],
                    true_X_next[IndX::K],
                    true_X_next[IndX::Ka]
                };

                // Measurement Y
                true_Y.col(i) = observationModel( true_X.col(i) );
                double U1 = random_0_1( generator ); // Случайное число от 0 до 1
                if( ( U1 > Pg ) || ( Pg == 0.0 ) ) {
                    measured_Y.col(i) = true_Y.col(i) + arma::vec{
                        noise_Rm_Y_R( generator ),
                        noise_Rm_Y_Az( generator ),
                        noise_Rm_Y_Vr( generator )            
                    };
                } else {
                    measured_Y.col(i) = true_Y.col(i) + arma::vec{
                        noise_Rg_Y_R( generator ),
                        noise_Rg_Y_Az( generator ),
                        noise_Rg_Y_Vr( generator )            
                    };
                }

                // Коррекция
                double randomNum = random_0_1( generator ); // Случайное число от 0 до 1
                if( ( p_ > randomNum ) || p_ == 1.0 ) { // || ( i < 2 )
                    
                    // #pragma omp parallel for                               
                    // for( auto &filterName : settings.Filters  ) {        
                    //    auto filter = filters.at( filterName );
                    for( int f = 0; f < settings.Filters.size(); f++ ) {            
                        auto filter = filters.at( settings.Filters[f] );

                        estimated_Y.at( filter->GetFilterName() ).col(i) = 
                            filter->Y(); // Прогнозный Y
                        arma::vec deltaY = measured_Y.col(i) - estimated_Y.at( filter->GetFilterName() ).col(i);
                        deltaY = checkDeltaMeasurement( deltaY );
                        filter->SetDeltaY( deltaY );

                        // Коррекция
                        timerCorrection.at( filter->GetFilterName() ).StartTimer();
                        filter->Correction( measured_Y.col(i) );
                        timerCorrection.at( filter->GetFilterName() ).EndTimer();

                        timeMeasuredCorrection.at( filter->GetFilterName() )[i] = 
                            timerCorrection.at( filter->GetFilterName() ).TimeCur() * 1.0e6; // мкс

                        arma::mat S = filter->GetInnovationCovarianceMatrixS();
                        if( filter->GetCovType() == KalmanFilters::TCovType::CT_SquareRootMatrices ) {
                            S = S * arma::trans( S );
                        }
                        estimated_Sdiag.at( filter->GetFilterName() ).col(i) = arma::sqrt( arma::abs( S.diag() ) );
                        arma::mat Sinv = arma::inv( S );
                        arma::mat lam = arma::trans( deltaY ) * Sinv * deltaY;
                        double md = std::sqrt( lam[0] );
                        mahalanobis.at( filter->GetFilterName() )[i] = md;

                        cycle++;
                        print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                    }
                } // end if( ( p_ > randomNum ) ) - Correction
                
                // #pragma omp parallel for           
                // for( auto &filterName : settings.Filters  ) {        
                //    auto filter = filters.at( filterName );
                for( int f = 0; f < settings.Filters.size(); f++ ) {            
                    auto filter = filters.at( settings.Filters[f] );

                    timeMeasuredSumm.at( filter->GetFilterName() )[i] = 
                        timeMeasuredPrediction.at( filter->GetFilterName() )[i] + 
                        timeMeasuredCorrection.at( filter->GetFilterName() )[i];
                    estimated_X.at( filter->GetFilterName() ).col(i) = filter->X();
                    estimated_Y.at( filter->GetFilterName() ).col(i) = filter->Y();

                    arma::mat estimatedP = filter->GetEstimatedCovarianceMatrixP();
                    if( filter->GetCovType() == KalmanFilters::TCovType::CT_FullMatrices ) {
                        estimatedP = arma::chol( estimatedP, "lower" );
                    }

                    estimated_Pdiag.at( filter->GetFilterName() ).col(i) = estimatedP.diag();
                    mean_Pdiag.at( filter->GetFilterName() ).col(i) += estimatedP.diag();
                    SDCM.at( filter->GetFilterName() )[i] = arma::trace( estimatedP );

                    delta_Y.at( filter->GetFilterName() ).col(i) = filter->GetDeltaY();
                    RMSE_X.at( filter->GetFilterName() ).col(i) += 
                        arma::square( 
                            checkDeltaState( estimated_X.at( filter->GetFilterName() ).col(i) - true_X.col(i) ) 
                        );
                }
            } // end Цикл по тактам времени
            //------------------------------------------------------------------
            // Суммирование времени для усреднения
            for( auto &filter : settings.Filters ) {
                timerPredictionAverage.at( filter ) += timerPrediction.at( filter ).TimePerOp() * 1.0e6; // Прогноз
                timerCorrectionAverage.at( filter ) += timerCorrection.at( filter ).TimePerOp() * 1.0e6; // Коррекция
            }
            //------------------------------------------------------------------
            // Построение графиков
            //------------------------------------------------------------------
            if( settings.MCruns == 1 ) {
                std::cout << "\nПостроение графиков..." << std::endl;
            }

            //------------------------------------------------------------------
            // Укорочение времени
            //------------------------------------------------------------------
            true_X.reshape( SizeX, Nfixed );
            true_Y.reshape( SizeY, Nfixed );
            measured_Y.reshape( SizeY, Nfixed );
            time.resize( Nfixed );
            for( auto &filter : settings.Filters ) {
                ( estimated_X.at( filter ) ).reshape( SizeX, Nfixed );
                ( estimated_Y.at( filter ) ).reshape( SizeY, Nfixed );
                ( estimated_Pdiag.at( filter ) ).reshape( SizeX, Nfixed );
                ( estimated_Sdiag.at( filter ) ).reshape( SizeY, Nfixed );
                ( delta_Y.at( filter ) ).reshape( SizeY, Nfixed );
                ( mahalanobis.at( filter ) ).resize( Nfixed );
                ( SDCM.at( filter ) ).resize( Nfixed );
                ( RMSE_X.at( filter ) ).reshape( SizeX, Nfixed );
                ( mean_Pdiag.at( filter ) ).reshape( SizeX, Nfixed );
                ( RMSE_X.at( filter ) ).reshape( SizeX, Nfixed );                
                ( timeMeasuredPrediction.at( filter ) ).resize( Nfixed );
                ( timeMeasuredCorrection.at( filter ) ).resize( Nfixed );
                ( timeMeasuredSumm.at( filter ) ).resize( Nfixed );
            }

            if( settings.Graphs_0_RMSE_1 == 0 ) {
            if( settings.Debug ) {
                std::cout << "seed = " << seed_ << "/" << settings.MCruns << " P = " << p_ << std::endl;
            }
            //------------------------------------------------------------------
            // X-Y, R-Az            
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::figure_size( pixels_width, pixels_height );
                plt::subplot( 1, 2, 1 );
            }
            plt::title( "а) Координаты X-Y" );
            plt::xlabel( "X, км" );
            plt::ylabel( "Y, км" );
            plt::plot(
                std::vector<double>{ ( true_X.row( IndX::X ) )[0] },
                std::vector<double>{ ( true_X.row( IndX::Y ) )[0] },
                true_start );
            plt::plot(
                arma::conv_to< std::vector<double> >::from( true_X.row( IndX::X ) ),
                arma::conv_to< std::vector<double> >::from( true_X.row( IndX::Y ) ),
                true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot(
                    arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row( IndX::X ) ),
                    arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row( IndX::Y ) ),
                    estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
            if( ylim_yes ) {
                plt::xlim( 98.0, 103.0 );
                plt::ylim( 197.0, 203.0 );
            }        
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
            plt::title( "б) Координаты R-Az" );
            plt::xlabel( "Az, градусы" );
            plt::ylabel( "R, км" );
            plt::plot(
                arma::conv_to< std::vector<double> >::from( measured_Y.row( IndY::Az ) ),
                arma::conv_to< std::vector<double> >::from( measured_Y.row( IndY::R ) ),
                marks_keywords );
            plt::plot(
                arma::conv_to< std::vector<double> >::from( true_Y.row( IndY::Az ) ),
                arma::conv_to< std::vector<double> >::from( true_Y.row( IndY::R ) ),
                true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot(
                    arma::conv_to< std::vector<double> >::from( ( estimated_Y.at( filter ) ).row( IndY::Az ) ),
                    arma::conv_to< std::vector<double> >::from( ( estimated_Y.at( filter ) ).row( IndY::R ) ),
                        estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
            if( ylim_yes ) {
                plt::xlim( 26.4, 27.0 );
                plt::ylim( 220.0, 227.5 );
            }        
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
            //------------------------------------------------------------------
            // mahalanobis 2
            plt::figure_size( pixels_width, pixels_height );

            plt::title( "Расстояние Махаланобиса" );
            plt::xlabel( "Время, с" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, mahalanobis.at( filter ), { { "color", graphColors.at( filter ) }, { "linestyle", "-" },
                    { "linewidth", "2" }, { "label", filter } } );
            }

            double argChiSquared = 4.0; //0.997
            plt::axhline( argChiSquared, 0, settings.SimulationTime - settings.DeltaT, { { "color", "black" }, { "linestyle", "--" },
                { "linewidth", "2" }, { "label", "Порог по уровню вероятности ~0.997" } } );
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
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
            //------------------------------------------------------------------
            // State 3
            //----------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::figure_size( pixels_width, pixels_height );
                plt::subplot( 3, 2, 1 );
            }
            plt::title( "а) Оценка координаты X" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "X, км");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row( IndX::X ) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row( IndX::X ) ), estimated_keywords.at( filter ) );
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
        //        plt::ylim( 98.0, 103.0 );
            if( ylim_yes ) {
                plt::ylim( 99.7, 101.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State_X" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //----------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 2 );
            }
            plt::title( "б) Оценка координаты Y");
            plt::xlabel( "Время, с" );
            plt::ylabel( "Y, км");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row( IndX::Y ) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row( IndX::Y ) ), estimated_keywords.at( filter ) );
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
        //        plt::ylim( 197.0, 203.0 );
            if( ylim_yes ) {
                plt::ylim( 199.5, 201.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State_Y" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //----------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 3 );
            }
            plt::title( "в) Оценка полной скорости");
            plt::xlabel( "Время, с" );
            plt::ylabel( "V, м/с");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row( IndX::V ) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row( IndX::V ) ), estimated_keywords.at( filter ) );
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
        //        plt::ylim( 94.0, 130.0 );
            if( ylim_yes ) {
                plt::ylim( 94.0, 120.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State_V" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //----------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 4 );
            }
            plt::title( "д) Оценка курса");
            plt::xlabel( "Время, с" );
            plt::ylabel( "К, градусы");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row( IndX::K ) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row( IndX::K ) ), estimated_keywords.at( filter ) );
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( -1.0, 60.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_3_State_K" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //----------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 6 );
            }
            plt::title( "е) Оценка скорости изменения курса" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "dK/dt, градусы/с");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( true_X.row( IndX::Ka ) ), true_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_X.at( filter ).row( IndX::Ka ) ), estimated_keywords.at( filter ) );
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true) ;
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( -0.2, 0.3 );
            }

            if( !settings.GraphSeparated ) {        
                plt::legend( loc, settings.LocLegend );
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
            //------------------------------------------------------------------
            // P 4
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::figure_size( pixels_width, pixels_height );
                plt::subplot( 3, 2, 1 );
            }
            plt::title( "СКО координаты X");
            plt::xlabel( "Время, с" );
            plt::ylabel( "СКО X, км");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row( IndX::X ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.01, 0.03 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P_X" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //----------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 2 );
            }
            plt::title( "СКО координаты Y");
            plt::xlabel( "Время, с" );
            plt::ylabel( "СКО Y, км");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row( IndX::Y ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.01, 0.03 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P_Y" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //----------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 3 );
            }
            plt::title( "СКО полной скорости");
            plt::xlabel( "Время, с" );
            plt::ylabel( "СКО V, м/с");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row( IndX::V ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 4.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P_V" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //----------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 4 );
            }
            plt::title( "СКО курса");
            plt::xlabel( "Время, с" );
            plt::ylabel( "СКО K, градусы");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row( IndX::K ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 12.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_4_sqrt_P_K" + name_seed_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //----------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 6 );
            }
            plt::title("СКО скорости изменения курса");
            plt::xlabel( "Время, с" );
            plt::ylabel( "СКО dK/dt, градусы/с");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Pdiag.at( filter ).row( IndX::Ka ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
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
            //------------------------------------------------------------------
            // Measurement 5
            plt::figure_size( pixels_width, pixels_height );

            plt::subplot( 3, 1, 1 );
            plt::title( "Дальность отметки" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "R, км" );
            plt::plot( time, arma::conv_to< std::vector<double> >::from( measured_Y.row( IndY::R ) ), marks_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Y.at( filter ).row( IndY::R ) ), estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );

            plt::subplot( 3, 1, 2 );
            plt::title( "Азимут отметки" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "Az, градусы");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( measured_Y.row( IndY::Az ) ), marks_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Y.at( filter ).row( IndY::Az ) ), estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );

            plt::subplot( 3, 1, 3 );
            plt::title( "Радиальная скорость отметки");
            plt::xlabel( "Время, с" );
            plt::ylabel( "Vr, м/с");
            plt::plot( time, arma::conv_to< std::vector<double> >::from( measured_Y.row( IndY::Vf ) ), marks_keywords );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Y.at( filter ).row( IndY::Vf ) ), estimated_keywords.at( filter ) );
            }
            plt::legend();
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );

            plt::subplots_adjust( keywords );
            name_tmp = "./" + filters_names + "_5_Measurement" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );

//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //------------------------------------------------------------------
            // S matrix
            plt::figure_size( pixels_width, pixels_height );

            plt::subplot( 3, 1, 1 );
            plt::title( "S, корень диаг. элемент дальности" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "dR, мc" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Sdiag.at( filter ).row( IndY::R ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );

            plt::subplot( 3, 1, 2 );
            plt::title( "S, корень диаг. элемент азимута" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "dAz, градусы" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Sdiag.at( filter ).row( IndY::Az ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );

            plt::subplot( 3, 1, 3 );
            plt::title( "S, корень диаг. элемент рад.скорости" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "dVr, м/с");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( estimated_Sdiag.at( filter ).row( IndY::Vf ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );

            plt::subplots_adjust( keywords );
            name_tmp = "./" + filters_names + "_6_sqrtS" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );
//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //------------------------------------------------------------------
            // deltaY
            plt::figure_size( pixels_width, pixels_height );

            plt::subplot( 3, 1, 1 );
            plt::title( "Невязка дальности" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "dR, км" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( delta_Y.at( filter ).row( IndY::R ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );

            plt::subplot( 3, 1, 2 );
            plt::title( "Невязка азимута" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "dAz, градусы" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( delta_Y.at( filter ).row( IndY::Az ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );

            plt::subplot( 3, 1, 3 );
            plt::title( "Невязка радиальной скорости" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "dVr, м/с");
            for( auto &filter : settings.Filters ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( delta_Y.at( filter ).row( IndY::Vf ) ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );

            plt::subplots_adjust( keywords );
            name_tmp = "./" + filters_names + "_5_deltaY" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );
//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //------------------------------------------------------------------
            // SDCM
            plt::figure_size( pixels_width, pixels_height );

            plt::title("SDCM");
            plt::xlabel( "Время, с" );
            plt::ylabel( "SDCM" );
            for( auto &filter : settings.Filters ) {
                plt::plot( time, SDCM.at( filter ), estimated_keywords.at( filter ) );
            }
            plt::legend( legend_loc );
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            name_tmp = "./" + filters_names + "_8_SDCM" + name_seed_p + "." + settings.Format;
            plt::save( name_tmp, dpi );
//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();

            //------------------------------------------------------------------
            // Time
            plt::figure_size( pixels_width, pixels_height );
            plt::title("Время выполнения прогноза/коррекции");
            plt::ylabel( "Время, мкс" );

            std::map<std::string, std::vector<double>> y;
            std::map<std::string, std::vector<double>> x;
            std::vector<double> x_ticks;

            double dy = settings.Filters.size() + 1;
            double fi = 0.0;
            for( auto &filter : settings.Filters ) {
                y.insert( std::make_pair( filter, std::vector<double>{
                    timerPrediction.at( filter ).TimePerOp() * 1.0e6, // Прогноз
                    timerCorrection.at( filter ).TimePerOp() * 1.0e6 // Коррекция
                    // ( timerPrediction.at( filter ).TimePerOp() + timerCorrection.at( filter ).TimePerOp() ) * 1.0e6 // Суммарно
                } ) );
                x.insert( std::make_pair( filter, std::vector<double>{
                    ( 0.0 + fi ),       // Прогноз
                    ( dy + fi )        // Коррекция
                    // ( 2.0 * dy ) + fi   // Суммарно
                } ) );
                plt::bar( x.at( filter ), y.at( filter ), "black", "-", 1.0, estimated_keywords_time.at( filter ) );
                fi += 1.0;
            }
            double first_tick = ( settings.Filters.size() - 1.0 ) / 2.0;
            for( int t = 0; t < 2; t++ ) { // 3 группы - Прогноз, Коррекция, Суммарно
                x_ticks.push_back( first_tick + ( t * dy ) );
            }
            std::vector<std::string> underlabels = { "Экстраполяция", "Коррекция" };
            plt::xticks( x_ticks, underlabels );

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
        
        TotalWorkTime.EndTimer();
        std::cout << "TotalWorkTime: " << TotalWorkTime.TimeSumm() << std::endl;
        
        //----------------------------------------------------------------------
        // RMSE
        if( settings.Graphs_0_RMSE_1 == 1 ) {

            std::cout << "\nПостроение графиков..." << std::endl;

            if( settings.Debug ) {
                std::cout << "P = " << p_ << std::endl;
            }

            for( auto &filter : settings.Filters ) {
                RMSE_X.at( filter ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
                RMSE_X.at( filter ) = arma::sqrt( RMSE_X.at( filter ) );
                mean_Pdiag.at( filter ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
            }
            //--------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::figure_size( pixels_width, pixels_height );
                plt::subplot( 3, 2, 1 );
                plt::subplots_adjust( keywords );
            }
            plt::title( "а) RMSE координаты X" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "RMSE X, км");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::X ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::X ) / RMSE_X.at( "EKF" ).row( IndX::X ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::X ) / ( RMSE_X.at( filter ).row( IndX::X ) )[0]  ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 0.03 ); // 0.05
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //--------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 2 );
            }
            plt::title( "б) RMSE координаты Y");
            plt::xlabel( "Время, с" );
            plt::ylabel( "RMSE Y, км");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::Y ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::Y ) / RMSE_X.at( "EKF" ).row( IndX::Y ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::Y ) / ( RMSE_X.at( filter ).row( IndX::Y )(0) ) ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 0.03 ); // 0.05
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_Y_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //--------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 3 );
            }
            plt::title( "в) RMSE полной скорости");
            plt::xlabel( "Время, с" );
            plt::ylabel( "RMSE V, м/с");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::V ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::V ) / RMSE_X.at( "EKF" ).row( IndX::V ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::V ) / ( RMSE_X.at( filter ).row( IndX::V )(0) ) ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 10.0 ); // 20.0
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_V_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //--------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 4 );
            }
            plt::title( "д) RMSE курса");
            plt::xlabel( "Время, с" );
            plt::ylabel( "RMSE К, градусы");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::K ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::K ) / RMSE_X.at( "EKF" ).row( IndX::K ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::K ) / ( RMSE_X.at( filter ).row( IndX::K )(0) ) ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 25.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_RMSE_X_K_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //--------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 6 );
            }
            plt::title( "е) RMSE скорости изменения курса" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "RMSE dK/dt, градус/с");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::Ka ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::Ka ) / RMSE_X.at( "EKF" ).row( IndX::Ka ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::Ka ) / ( RMSE_X.at( filter ).row( IndX::Ka )(0) ) ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true) ;
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 0.18 ); // 0.7
            }

            if( !settings.GraphSeparated ) {        
                plt::legend( loc, settings.LocLegend );
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
            //----------------------------------------------------------------------------------------------------------
            // Mean P
            //----------------------------------------------------------------------------------------------------------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::figure_size( pixels_width, pixels_height );
                plt::subplot( 3, 2, 1 );
            }
            plt::title( "а) Mean P координаты X" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "Mean P_X, км");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::X ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::X ) / mean_Pdiag.at( "EKF" ).row( IndX::X ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::X ) / ( RMSE_X.at( filter ).row( IndX::X ) )[0]  ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 0.03 ); // 0.05
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_Mean_P_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //--------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 2 );
            }
            plt::title( "б) Mean P координаты Y");
            plt::xlabel( "Время, с" );
            plt::ylabel( "Mean P_Y, км");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::Y ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::Y ) / mean_Pdiag.at( "EKF" ).row( IndX::Y ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::Y ) / ( RMSE_X.at( filter ).row( IndX::Y )(0) ) ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 0.03 ); // 0.05
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_Mean_P_Y_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //--------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 3 );
            }
            plt::title( "в) Mean P полной скорости");
            plt::xlabel( "Время, с" );
            plt::ylabel( "Mean P_V, м/с");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::V ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::V ) / mean_Pdiag.at( "EKF" ).row( IndX::V ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::V ) / ( RMSE_X.at( filter ).row( IndX::V )(0) ) ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 10.0 ); // 20.0
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_Mean_P_V_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //--------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 5 );
            }
            plt::title( "д) Mean P курса");
            plt::xlabel( "Время, с" );
            plt::ylabel( "Mean P_К, градусы");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::K ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::K ) / mean_Pdiag.at( "EKF" ).row( IndX::K ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::K ) / ( RMSE_X.at( filter ).row( IndX::K )(0) ) ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true );
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 25.0 );
            }
            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_Mean_P_K_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }
            //--------
            if( settings.GraphSeparated ) {
                plt::figure_size( pixels_width, pixels_height );
            } else {
                plt::subplot( 3, 2, 6 );
            }
            plt::title( "е) Mean P скорости изменения курса" );
            plt::xlabel( "Время, с" );
            plt::ylabel( "Mean P_dK/dt, градус/с");
            for( auto &filter : settings.Filters ) {
                if( !relatEKF ) {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::Ka ) ), estimated_keywords.at( filter ) );
                } else {
                    plt::plot( time, arma::conv_to< std::vector<double> >::from( mean_Pdiag.at( filter ).row( IndX::Ka ) / mean_Pdiag.at( "EKF" ).row( IndX::Ka ) ), estimated_keywords.at( filter ) );
//                    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X.at( filter ).row( IndX::Ka ) / ( RMSE_X.at( filter ).row( IndX::Ka )(0) ) ), estimated_keywords.at( filter ) );
                }
            }
            if( settings.GraphSeparated ) {
                plt::legend( legend_loc );
            }
            plt::grid( true) ;
//            plt::xlim( 0.0, settings.SimulationTime );
            plt::xlim( 0.0, Nfixed * settings.DeltaT );
            if( ylim_yes ) {
                plt::ylim( 0.0, 0.18 ); // 0.7
            }

            if( !settings.GraphSeparated ) {        
                plt::legend( loc, settings.LocLegend );
            }

            if( settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_Mean_P_Ka_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

            if( !settings.GraphSeparated ) {
                plt::subplots_adjust( keywords );
                name_tmp = "./" + filters_names + "_Mean_P_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_p + "." + settings.Format;
                plt::save( name_tmp, dpi );
            }

//            if( settings.ShowGraphs ) {
//                plt::show();
//            }
//            plt::close();
            //------------------------------------------------------------------
            // Время усреднённое
            for( auto &filter : settings.Filters ) {
                timerPredictionAverage.at( filter ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
                timerCorrectionAverage.at( filter ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
            }
            plt::figure_size( pixels_width, pixels_height );
            std::string my_title = "Время выполнения прогноза/коррекции усреднённое по " + std::to_string( settings.MCruns ) + " реализациям";
//            plt::title( my_title );
            plt::ylabel( "Время, мкс" );

            std::map<std::string, std::vector<double>> y;
            std::map<std::string, std::vector<double>> x;
            std::vector<double> x_ticks;

            double dy = settings.Filters.size() + 1;
            double fi = 0.0;
            for( auto &filter : settings.Filters ) {
                y.insert( std::make_pair( filter, std::vector<double>{
                    timerPredictionAverage.at( filter ), // Прогноз
                    timerCorrectionAverage.at( filter ) // Коррекция
                    // ( timerPredictionAverage.at( filter ) + timerCorrectionAverage.at( filter ) ) // Суммарно
                } ) );
                x.insert( std::make_pair( filter, std::vector<double>{
                    ( 0.0 + fi ),       // Прогноз
                    ( dy + fi )        // Коррекция
                    // ( 2.0 * dy ) + fi   // Суммарно
                } ) );
                plt::bar( x.at( filter ), y.at( filter ), "black", "-", 1.0, estimated_keywords_time.at( filter ) );
                fi += 1.0;
            }
            double first_tick = ( settings.Filters.size() - 1.0 ) / 2.0;
            for( int t = 0; t < 2; t++ ) { // 3 группы - Прогноз, Коррекция, Суммарно
                x_ticks.push_back( first_tick + ( t * dy ) );
            }
            std::vector<std::string> underlabels = { "Экстраполяция", "Коррекция" };
            plt::xticks( x_ticks, underlabels );

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
            //------------------------------------------------------------------
            // Время усреднённое относительно EKF
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "EKF" ) ) {

                plt::figure_size( pixels_width, pixels_height );
                my_title = "Время выполнения прогноза/коррекции (относительно EKF) усреднённое по " + std::to_string( settings.MCruns ) + " реализациям";
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
                    // double summ = ( timerPredictionAverage.at( filter ) + timerCorrectionAverage.at( filter ) ) /
                        // ( timerPredictionAverage.at( "EKF" ) + timerCorrectionAverage.at( "EKF" ) );

                    ofs <<
                        filter << "\t" <<
                        pred << "\t" <<
                        corr << "\t" <<
                        // summ << "\t" <<
                        std::endl;

                    y.insert( std::make_pair( filter, std::vector<double>{
                        pred, // Прогноз
                        corr // Коррекция
                        // summ // Суммарно
                    } ) );
                    x.insert( std::make_pair( filter, std::vector<double>{
                        ( 0.0 + fi ),       // Прогноз
                        ( dy + fi )        // Коррекция
                        // ( 2.0 * dy ) + fi   // Суммарно
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
                for( int t = 0; t < 2; t++ ) { // 3 группы - Прогноз, Коррекция, Суммарно
                    x_ticks.push_back( first_tick + ( t * dy ) );
                }
                std::vector<std::string> underlabels2 = { "Экстраполяция", "Коррекция" };
                plt::xticks( x_ticks, underlabels2 );

    //            plt::legend( legend_loc_time );
                plt::legend( loc_time, { 1, 1 } );
                plt::grid( true );
                plt::subplots_adjust( keywords );
                plt::save( name_tmp, dpi );
            }

            if( settings.ShowGraphs ) {
                plt::show();
            }
            plt::close();
        } // end if( settings.Graphs_0_RMSE_1 == 1 ) {
    } // p - по вероятностям
    //--------------------------------------------------------------------------
    plt::clf();
    plt::cla();
    Py_Finalize();

    // delete filters;

    // TotalWorkTime.EndTimer();
    // std::cout << "TotalWorkTime: " << TotalWorkTime.TimeSumm() << std::endl;

    }
    catch( std::exception &ex )
    {
        std::cout << "Exception occured: " << ex.what() << std::endl;
    }
    std::cout << "Завершено." << std::endl;
}

/// \}

