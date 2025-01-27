//------------------------------------------------------------------------------
///
/// \file       compare_run_rmse.cpp
/// \brief      Сравнение фильтров Калмана
/// \date       26.12.24 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#include <compare_run.h>

void CKalmanFiltersCompare::Run_RMSE_var_params( const CSettings &settings )
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
    std::string filters_names = "COMPARE_EKF";

    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        filters_names += "_SRUKF";
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        filters_names += "_SREUKF";
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
        filters_names += "_SRCKF";
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
        filters_names += "_SRECKF";
    }
    std::map<std::string, std::string>legend_loc { { "loc", "upper right" } };
    // std::string loc = "center right";
    std::string loc = "lower center";

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

    arma::mat RMSE_X_EKF = template_mat_X;
    arma::mat RMSE_X_SRCKF = template_mat_X;
    arma::mat RMSE_X_SRECKF = template_mat_X;
    std::map<std::string, arma::mat> RMSE_X_SRUKF;
    std::map<std::string, arma::mat> RMSE_X_SREUKF;

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

    //--------------------------------------------------------------------------
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

    EKF.SetProcessCovarianceMatrixQdiag( Q_CV % arma::vec( settings.q_koef_ekf ) );    
    EKF.SetObservationCovarianceMatrixRdiag( R );
    //--------------------------------------------------------------------------
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

    SRUKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q_CV % arma::vec( settings.q_koef_ukf ) ) );    
    SRUKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    //--------------------------------------------------------------------------
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

    SREUKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q_CV % arma::vec( settings.q_koef_ukf ) ) );
    SREUKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    //--------------------------------------------------------------------------    
    KalmanFilters::CKalmanSRCKF<SizeX, SizeY> SRCKF;
    SRCKF.SetStateTransitionModel( stateTransitionModel );
    SRCKF.SetObservationModel( observationModel );

    SRCKF.SetCheckBordersStateAfterPrediction( checkBordersState );
    SRCKF.SetCheckBordersStateAfterCorrection( checkBordersState );
    SRCKF.SetCheckBordersMeasurement( checkBordersMeasurement );
    SRCKF.SetCheckDeltaState( checkDeltaState );
    SRCKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

    SRCKF.SetWeightedSumStateSigmas( weightedSumStateSigmas );
    SRCKF.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

    SRCKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q_CV % arma::vec( settings.q_koef_ukf ) ) );    
    SRCKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    SRCKF.SetDesignParametersCubatureBaseSet();
    //--------------------------------------------------------------------------    
    KalmanFilters::CKalmanSRECKF<SizeX, SizeY> SRECKF;
    SRECKF.SetStateTransitionModel( stateTransitionModel );
    SRECKF.SetObservationModel( observationModel );
    SRECKF.SetStateTransitionJacobianF( stateTransitionJacobianF );
    SRECKF.SetObservationJacobianH( observationJacobianH );

    SRECKF.SetCheckBordersStateAfterPrediction( checkBordersState );
    SRECKF.SetCheckBordersStateAfterCorrection( checkBordersState );
    SRECKF.SetCheckBordersMeasurement( checkBordersMeasurement );
    SRECKF.SetCheckDeltaState( checkDeltaState );
    SRECKF.SetCheckDeltaMeasurement( checkDeltaMeasurement );

    SRECKF.SetWeightedSumStateSigmas( weightedSumStateSigmas );
    SRECKF.SetWeightedSumMeasurementSigmas( weightedSumMeasurementSigmas );

    SRECKF.SetProcessCovarianceMatrixQdiag( arma::sqrt( Q_CV % arma::vec( settings.q_koef_ukf ) ) );    
    SRECKF.SetObservationCovarianceMatrixRdiag( arma::sqrt( R ) );
    SRECKF.SetDesignParametersCubatureBaseSet();
    //--------------------------------------------------------------------------            
    if( settings.Set == 0 ) { // Julier
        filters_names += "_Julier";
        for( std::size_t i = 0; i < settings.w0_sr.size(); i++ ) {
            SRUKF.SetDesignParametersMeanSet( settings.w0_sr[i] );            
            SREUKF.SetDesignParametersMeanSet( settings.w0_sr[i] );

            // std::string key = "w0=" + std::to_string( settings.w0_sr[i] );
            std::ostringstream oss;
            oss << "w0=" << std::setw(3) << settings.w0_sr[i];
            std::string key = oss.str();
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                auto srukf = SRUKF;
                filtersSRUKF.insert( std::make_pair( key, srukf ) );
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                auto sreukf = SREUKF;
                filtersSREUKF.insert( std::make_pair( key, sreukf ) );
            }
        }
    } else if( settings.Set == 1 ) { // Merwe
        filters_names += "_Merwe";
        for( std::size_t a = 0; a < settings.alpha_sr.size(); a++ )
        for( std::size_t b = 0; b < settings.beta_sr.size(); b++ )
        for( std::size_t k = 0; k < settings.kappa_sr.size(); k++ ) {
            SRUKF.SetDesignParametersScaledSet( settings.alpha_sr[a], settings.beta_sr[b], settings.kappa_sr[k] );            
            SREUKF.SetDesignParametersScaledSet( settings.alpha_sr[a], settings.beta_sr[b], settings.kappa_sr[k] );

            // std::string key = "a_b_k=" +
            //     std::to_string( settings.alpha_sr[a] ) + "_" +
            //     std::to_string( settings.beta_sr[b] ) + "_" +
            //     std::to_string( settings.kappa_sr[k] );
            std::ostringstream oss;
            oss << "SRUKF a, b, k =" << std::setw(3) << 
                settings.alpha_sr[a] << ", " <<
                settings.beta_sr[b] << ", " <<
                settings.kappa_sr[k];
            std::string key = oss.str();
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                auto srukf = SRUKF;
                filtersSRUKF.insert( std::make_pair( key, srukf ) );
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                auto sreukf = SREUKF;
                filtersSREUKF.insert( std::make_pair( key, sreukf ) );
            }
        }
    } else {
        assert( false );
    }
    //--------------------------------------------------------------------------
    // Начало рабочих циклов
    //--------------------------------------------------------------------------

    std::string name_dt = "_dt_" + std::to_string( settings.DeltaT );
                                         //EKF                                              //SRUKF2
    int cycle_max = settings.MCruns * N * ( 1 + filtersSRUKF.size() + filtersSREUKF.size() + 1 );
    int cycle = 0;
    int prev_percent = -1;

    //--------------------------------------------------------------------------
    for( uint32_t seed_ = settings.MCseed; seed_ < ( settings.MCseed + settings.MCruns ); seed_++ ) {
        if( ( seed_ == 0 ) && ( settings.MCruns == 1 ) ) {
            auto rand_seed = rd();
            generator.seed( rand_seed ); // Выставить зерно ГСЧ случайным, если seed = 0
            std::cout << "random seed = " << rand_seed << std::endl;
        } else {
            generator.seed( seed_ ); // Выставить зерно ГСЧ
        }
        std::string name_seed_p = "_seed_" + std::to_string( seed_ ) + "_" + name_dt;
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
        //----------------------------------------------------------------------
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
        RMSE_X_EKF.col(0) += arma::square( checkDeltaState( startX - true_X.col(0) ) );
        //----------------------------------------------------------------------
        // SRCKF
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
        RMSE_X_SRCKF.col(0) += arma::square( checkDeltaState( startX - true_X.col(0) ) );
        //----------------------------------------------------------------------
        // SRECKF
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
        RMSE_X_SRECKF.col(0) += arma::square( checkDeltaState( startX - true_X.col(0) ) );
        //----------------------------------------------------------------------
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
                PdenseSRUKF.diag() = startP_UKF;
                PdenseSRUKF = arma::chol( PdenseSRUKF, "lower" );
                filter.SetEstimateCovarianceMatrixP( PdenseSRUKF );
                if( settings.Debug ) {
                    PdenseSRUKF.print("PdenseSRUKF:");
                }

                RMSE_X_SRUKF.insert( std::make_pair( key, template_mat_X ) );

                RMSE_X_SRUKF.at( key ).col(0) += arma::square( checkDeltaState( startX - true_X.col(0) ) );
            }
        }
        //----------------------------------------------------------------------
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
                PdenseSREUKF.diag() = startP_UKF;
                PdenseSREUKF = arma::chol( PdenseSREUKF, "lower" );
                filter.SetEstimateCovarianceMatrixP( PdenseSREUKF );
                if( settings.Debug ) {
                    PdenseSREUKF.print("PdenseSREUKF:");
                }

                RMSE_X_SREUKF.insert( std::make_pair( key, template_mat_X ) );

                RMSE_X_SREUKF.at( key ).col(0) += arma::square( checkDeltaState( startX - true_X.col(0) ) );
            }
        }
        //------------------------------------------------------------------
        // Симуляция
        for( int i = 1; i < N; i++ ) { // Цикл по тактам времени
            time[i] = i * settings.DeltaT;

            // Прогноз
            EKF.Prediction( settings.DeltaT );
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
                SRCKF.Prediction( settings.DeltaT );
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {                
                SRECKF.Prediction( settings.DeltaT );
            }

            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                for( auto &item : filtersSRUKF ) {
                    auto &filter = item.second;
                    filter.Prediction( settings.DeltaT );
                }
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                for( auto &item : filtersSREUKF ) {
                    auto &filter = item.second;
                    filter.Prediction( settings.DeltaT );
                }
            }

            // Создание измерений текущего такта:
            true_X.col(i) = stateTransitionModel( true_X.col(i - 1), settings.DeltaT );

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
            arma::vec deltaY;
            deltaY = measured_Y.col(i) - EKF.Y();
            deltaY = checkDeltaMeasurement( deltaY );
            EKF.SetDeltaY( deltaY );
            EKF.Correction( measured_Y.col(i) );            

            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
                deltaY = measured_Y.col(i) - SRCKF.Y();
                deltaY = checkDeltaMeasurement( deltaY );
                SRCKF.SetDeltaY( deltaY );
                SRCKF.Correction( measured_Y.col(i) );            
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {                
                deltaY = measured_Y.col(i) - SRECKF.Y();
                deltaY = checkDeltaMeasurement( deltaY );
                SRECKF.SetDeltaY( deltaY );
                SRECKF.Correction( measured_Y.col(i) );            
            }
            cycle++;
            print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения

            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                for( auto &item : filtersSRUKF ) {
                    auto &filter = item.second;
                    deltaY = measured_Y.col(i) - filter.Y();
                    deltaY = checkDeltaMeasurement( deltaY );
                    filter.SetDeltaY( deltaY );
                    filter.Correction( measured_Y.col(i) );
                    cycle++;
                    print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                }
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                for( auto &item : filtersSREUKF ) {
                    auto &filter = item.second;
                    deltaY = measured_Y.col(i) - filter.Y();
                    deltaY = checkDeltaMeasurement( deltaY );
                    filter.SetDeltaY( deltaY );
                    filter.Correction( measured_Y.col(i) );
                    cycle++;
                    print_percent( cycle, cycle_max, prev_percent ); // Напечатать проценты выполнения
                }
            }

            // RMSE
            RMSE_X_EKF.col(i) += arma::square( checkDeltaState( EKF.X() - true_X.col(i) ) );            
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
                RMSE_X_SRCKF.col(i) += arma::square( checkDeltaState( SRCKF.X() - true_X.col(i) ) );            
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
                RMSE_X_SRECKF.col(i) += arma::square( checkDeltaState( SRECKF.X() - true_X.col(i) ) );            
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
                for( auto &item : filtersSRUKF ) {
                    auto &key = item.first;
                    auto &filter = item.second;
                    RMSE_X_SRUKF.at( key ).col(i) += arma::square( checkDeltaState( filter.X() - true_X.col(i) ) );
                }
            }
            if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
                for( auto &item : filtersSREUKF ) {
                    auto &key = item.first;
                    auto &filter = item.second;
                    RMSE_X_SREUKF.at( key ).col(i) += arma::square( checkDeltaState( filter.X() - true_X.col(i) ) );
                }
            }
        } // end Цикл по тактам времени
    } // seed
    //----------------------------------------------------------------------
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

    RMSE_X_SRCKF *= ( 1.0 / static_cast<double>( settings.MCruns ) );
    RMSE_X_SRCKF = arma::sqrt( RMSE_X_SRCKF );

    RMSE_X_SRECKF *= ( 1.0 / static_cast<double>( settings.MCruns ) );
    RMSE_X_SRECKF = arma::sqrt( RMSE_X_SRECKF );

    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
            RMSE_X_SRUKF.at( key ) *= ( 1.0 / static_cast<double>( settings.MCruns ) );
            RMSE_X_SRUKF.at( key ) = arma::sqrt( RMSE_X_SRUKF.at( key ) );
        }
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
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
    plt::title( "а) RMSE координаты X" );
    plt::xlabel( "Время, с" );
    plt::ylabel( "RMSE X, км");
    int i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            if( filtersSRUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::X ) ),
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::X ) ),
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "" }, { "label", "SRUKF" } } );
            }
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            // std::string color_hex_string = "#" + stream.str() + "0000"; // RED
            std::string color_hex_string = "#FFFF" + stream.str(); // YELLOW

            if( filtersSREUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::X ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::X ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SREUKF" } } );
            }
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row( IndX::X ) ),
        { { "color", "darkblue" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "EKF" } } );
    
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRCKF.row( IndX::X ) ),
            { { "color", "indianred" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRCKF" } } );
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRECKF.row( IndX::X ) ),
            { { "color", "mediumpurple" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRECKF" } } );
    }
    
    if( settings.GraphSeparated ) {
        plt::legend( legend_loc );
    }

    plt::grid( true );
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 0.03 );
    }
    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_X_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }

    if( settings.GraphSeparated ) {
        plt::figure_size( pixels_width, pixels_height );
    } else {
        plt::subplot( 3, 2, 2 );
    }
    plt::title( "б) RMSE координаты Y");
    plt::xlabel( "Время, с" );
    plt::ylabel( "RMSE Y, км");
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            if( filtersSRUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::Y ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::Y ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRUKF" } } );
            }
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            // std::string color_hex_string = "#" + stream.str() + "0000"; // RED
            std::string color_hex_string = "#FFFF" + stream.str(); // YELLOW
            if( filtersSREUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::Y ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::Y ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SREUKF" } } );
            }
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row( IndX::Y ) ), //estimated_keywords.at( filter ) );
        { { "color", "darkblue" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "EKF" } } );

    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRCKF.row( IndX::Y ) ),
            { { "color", "indianred" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRCKF" } } );
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRECKF.row( IndX::Y ) ),
            { { "color", "mediumpurple" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRECKF" } } );
    }        

    if( settings.GraphSeparated ) {
        plt::legend( legend_loc );
    }

    plt::grid( true );
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 0.03 );
    }
    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_Y_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }

    if( settings.GraphSeparated ) {
        plt::figure_size( pixels_width, pixels_height );
    } else {
        plt::subplot( 3, 2, 3 );
    }
    plt::title( "в) RMSE полной скорости");
    plt::xlabel( "Время, с" );
    plt::ylabel( "RMSE V, м/с");
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            if( filtersSRUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::V ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::V ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRUKF" } } );
            }
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            // std::string color_hex_string = "#" + stream.str() + "0000"; // RED
            std::string color_hex_string = "#FFFF" + stream.str(); // YELLOW

            if( filtersSREUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::V ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::V ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "1" }, { "label", "SREUKF" } } );
            }
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row( IndX::V ) ), //estimated_keywords.at( filter ) );
        { { "color", "darkblue" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "EKF" } } );

    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRCKF.row( IndX::V ) ),
            { { "color", "indianred" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRCKF" } } );
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRECKF.row( IndX::V ) ),
            { { "color", "mediumpurple" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRECKF" } } );
    }         

    if( settings.GraphSeparated ) {
        plt::legend( legend_loc );
    }

    plt::grid( true );
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 10.0 ); // 15
    }
    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_V_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }
    ///
    if( settings.GraphSeparated ) {
        plt::figure_size( pixels_width, pixels_height );
    } else {
        plt::subplot( 3, 2, 4 );
    }
    plt::title( "г) RMSE курса");
    plt::xlabel( "Время, с" );
    plt::ylabel( "RMSE К, градусы");
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            if( filtersSRUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::K ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::K ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRUKF" } } );
            }
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            // std::string color_hex_string = "#" + stream.str() + "0000"; // RED
            std::string color_hex_string = "#FFFF" + stream.str(); // YELLOW

            if( filtersSREUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::K ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::K ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SREUKF" } } );
            }
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row( IndX::K ) ), //estimated_keywords.at( filter ) );
        { { "color", "darkblue" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "EKF" } } );

    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRCKF.row( IndX::K ) ),
            { { "color", "indianred" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRCKF" } } );
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRECKF.row( IndX::K ) ),
            { { "color", "mediumpurple" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRECKF" } } );
    }         

    if( settings.GraphSeparated ) {
        plt::legend( legend_loc );
    }    

    plt::grid( true );
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 25.0 );
    }
    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_K_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
        plt::save( name_tmp, dpi );
    }


    if( settings.GraphSeparated ) {
        plt::figure_size( pixels_width, pixels_height );
    } else {
        plt::subplot( 3, 2, 6 );
    }
    plt::title( "д) RMSE скорости изменения курса" );
    plt::xlabel( "Время, с" );
    plt::ylabel( "RMSE dK/dt, градусы/с");
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRUKF" ) ) {
        for( auto &item : filtersSRUKF ) {
            auto &key = item.first;
            int color = color_SRUKF_start + color_SRUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            std::string color_hex_string = "#00" + stream.str() + "00"; // GREEN

            if( filtersSRUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::Ka ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRUKF.at( key ).row( IndX::Ka ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRUKF" } } );
            }
            i++;
        }
    }
    i = 0;
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SREUKF" ) ) {
        for( auto &item : filtersSREUKF ) {
            auto &key = item.first;
            int color = color_SREUKF_start + color_SREUKF_step * i;
            std::stringstream stream;
            stream << std::hex << color;
            // std::string color_hex_string = "#" + stream.str() + "0000"; // RED
            std::string color_hex_string = "#FFFF" + stream.str(); // YELLOW

            if( filtersSREUKF.size() > 1 ) {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::Ka ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", key } } );
            } else {
                plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SREUKF.at( key ).row( IndX::Ka ) ), //estimated_keywords.at( key ) );
                    { { "color", color_hex_string }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SREUKF" } } );
            }
            i++;
        }
    }
    plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_EKF.row( IndX::Ka ) ), //estimated_keywords.at( filter ) );
        { { "color", "darkblue" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "EKF" } } );
    
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRCKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRCKF.row( IndX::Ka ) ),
            { { "color", "indianred" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRCKF" } } );
    }
    if( std::count( settings.Filters.begin(), settings.Filters.end(), "SRECKF" ) ) {
        plt::plot( time, arma::conv_to< std::vector<double> >::from( RMSE_X_SRECKF.row( IndX::Ka ) ),
            { { "color", "mediumpurple" }, { "linestyle", "-" }, { "linewidth", "2" }, { "label", "SRECKF" } } );
    }           

    if( settings.GraphSeparated ) {
        plt::legend( legend_loc );
    }

    plt::grid( true) ;
    plt::xlim( 0.0, settings.SimulationTime );
    if( ylim_yes ) {
        plt::ylim( 0.0, 0.18 );
    }

    if( !settings.GraphSeparated ) {        
        plt::legend( loc, settings.LocLegend );
    }

    if( settings.GraphSeparated ) {
        plt::subplots_adjust( keywords );
        name_tmp = "./" + filters_names + "_Ka_" + "_RMSE_X_MCruns_" + std::to_string( settings.MCruns ) + "_" + name_dt + "." + settings.Format;
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
    //--------------------------------------------------------------------------
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

