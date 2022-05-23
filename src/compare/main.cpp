//----------------------------------------------------------------------------------------------------------------------
///
/// \file       main.cpp
/// \brief      Точка входа
/// \date       06.04.22 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

// System includes:
#include <iostream>
#include <boost/program_options.hpp>
#ifdef NDEBUG // We need assert even in a release build for test code.
# undef NDEBUG
#endif
#include <cassert>

// Project includes:
#include <compare.h>

///
/// \brief Перегрузка оператора << для std::vector<double>
/// \param os - поток, куда выводим
/// \param vec - выводимый вектор
/// \return Поток
///
std::ostream& operator<<( std::ostream &os, const std::vector<double> &vec )
{
    for( auto item : vec ) {
        os << item << " ";
    }
    return os;
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief main - Основная функция
/// \param argc - Количество аргументов командной строки
/// \param argv - Аргументы командной строки
/// \return 0 - штатная работа, 1 - ошибка
///
int main( int argc, char *argv[] )
{
    CKalmanFiltersCompare test; // Класс тестирования фильтров
    CSettings settings; // Настройки программы
//    arma::arma_version ver;
//    std::cout << "Armadillo library version: "<< ver.as_string() << std::endl;
    //------------------------------------------------------------------------------------------------------------------
    // Зададим параметры запуска приложения
    namespace po = boost::program_options;
    po::options_description desc( "Параметры приложения", 220 ); // 220 - задает ширину строки вывода в терминал
    desc.add_options()
        ( "help", "Показать параметры справки и выйти" )
        ( "debug", "Выводить номер такта работы и время" )
//        ( "seed,s", po::value<uint32_t>( &settings.Seed )->default_value( settings.Seed ),
//            "Задать зерно ГСЧ для имитационного сценария (unsigned int > 0, если = 0, то случайное зерно)" )
        ( "par", po::value<std::vector<double>>( &settings.MatPlotParams )->multitoken()->default_value( std::vector<double>{ 0.12, 0.2, 0.9, 0.88, 0.2, 0.2 }, "0.12, 0.2, 0.9, 0.88, 0.2, 0.2" ),
            "Задать параметры matplotlib: left, bottom, right, top, wspace, hspace" )
        ( "prob", po::value<std::vector<double>>( &settings.Probabilities )->multitoken()->default_value( std::vector<double>{1.0}, "1.0" ),
            "Задать вероятности целевой отметки" )
        ( "dt", po::value<double>( &settings.DeltaT )->default_value( settings.DeltaT ),
            "Задать шаг работы по времени в секундах" )
        ( "time", po::value<double>( &settings.SimulationTime )->default_value( settings.SimulationTime ),
            "Задать полное время имитации в секундах" )
        ( "graph", po::value<uint32_t>( &settings.Graphs_0_RMSE_1 )->default_value( settings.Graphs_0_RMSE_1 ),
            "0 - строить графики состояния, 1 - строить график RMSE" )
        ( "size", po::value<std::vector<double>>( &settings.Size )->multitoken()->default_value( std::vector<double>{ 195, 130 }, "195, 130" ),
            "Задать размеры графиков (ширина/высота) в мм при dpi=300" )
        ( "sep", "Раздельные графики" )
        ( "eps", "Сохранять в eps формате (иначе в png, если ключ не задан)" )
        ( "show", "Показать графики, потом записать в файлы (иначе только вывести в файлы, не показывая)" )
        ( "mc", po::value<uint32_t>( &settings.MCruns )->default_value( settings.MCruns ),
            "Число запусков реализаций ГСЧ" )
#ifdef EKF_
        ( "ekf", "Расширенный Фильтр Калмана, РФК (Extended Kalman Filter, EKF)" )
#endif     
#ifdef SREKF_
        ( "srekf", "Квадратно-Корневой Расширенный Фильтр Калмана, КК-РФК (Square Root Extended Kalman Filter, SREKF)" )
#endif
#ifdef UKF_
        ( "ukf", "Сигма-точечный (ансцентный) Фильтр Калмана, СТФК (Unscented Kalman Filter, UKF)" )
#endif
#ifdef SRUKF_
        ( "srukf", "Квадратно-Корневой Сигма-точечный Фильтр Калмана, КК-СТФК (Square Root Unscented Kalman Filter, SRUKF)" )
        ( "srukfb", "Квадратно-Корневой Сигма-точечный Фильтр Калмана (блочная реализация), КК-СТФКБ (Square Root Unscented Kalman Filter Block, SRUKFB)" )
#endif
#ifdef CKF_
        ( "ckf", "Кубатурный Фильтр Калмана, КФК (Cubature Kalman Filter, CKF)" )
#endif
#ifdef SRCKF_
        ( "srckf", "Квадратно-Корневой Кубатурный Фильтр Калмана, КК-КФК (Square Root Cubature Kalman Filter, SRCKF)" )
        ( "srckfb", "Квадратно-Корневой Кубатурный Фильтр Калмана (блочная реализация), КК-КФКБ (Square Root Cubature Kalman Filter Block, SRCKFB)" )
#endif
#ifdef ECKF_
        ( "eckf", "Расширенный Кубатурный Фильтр Калмана, РКФК (Extended Cubature Kalman Filter, ECKF)" )
#endif
#ifdef SRECKF_
        ( "sreckf", "Квадратно-Корневой Расширенный Кубатурный Фильтр Калмана, КК-РКФК (Square Root Extended Cubature Kalman Filter, SRECKF)" )
        ( "sreckfb", "Квадратно-Корневой Расширенный Кубатурный Фильтр Калмана (блочная реализация), КК-РКФКБ (Square Root Extended Cubature Kalman Filter Block, SRECKFB)" )
#endif
#ifdef EUKF_
        ( "eukf", "Расширенный Сигма-точечный (ансцентный) Фильтр Калмана, РСТФК (Extended Unscented Kalman Filter, EUKF)" )
#endif
#ifdef SREUKF_
        ( "sreukf", "Квадратно-Корневой Расширенный Сигма-точечный Фильтр Калмана, КК-РСТФК (Square Root Extended Unscented Kalman Filter, SREUKF)" )
        ( "sreukfb", "Квадратно-Корневой Расширенный Сигма-точечный Фильтр Калмана (блочная реализация), КК-РСТФКБ (Square Root Extended Unscented Kalman Filter Block, SREUKFB)" )
#endif
        ( "all", "Все фильтры" )

        ( "rmse", "построить RMSE для EKF и SRUKF или SREUKF отдельно (надо задать --srukf и/или --sreukf)" )

        ( "set", po::value<uint32_t>( &settings.Set )->default_value( settings.Set ),
            "Задать способ установки сигма-точек: 0-Julier, 1-Merwe" )        

        ( "w0", po::value<std::vector<double>>( &settings.w0 )->multitoken()->default_value( std::vector<double>{ 0.3 }, "0.3" ),
            "Задать параметр(ы) веса w0 для сигма-точечных фильтров" )
        ( "alpha", po::value<std::vector<double>>( &settings.alpha )->multitoken()->default_value( std::vector<double>{ 0.687369 }, "0.687369" ), //
            "Задать параметр(ы) alpha для сигма-точечных фильтров" )
        ( "beta", po::value<std::vector<double>>( &settings.beta )->multitoken()->default_value( std::vector<double>{ 2.0 }, "2" ),
            "Задать параметр(ы) beta для сигма-точечных фильтров" )
        ( "kappa", po::value<std::vector<double>>( &settings.kappa )->multitoken()->default_value( std::vector<double>{ -2.0 }, "-2" ),
            "Задать параметр(ы) kappa для сигма-точечных фильтров" )

        ( "w0_sr", po::value<std::vector<double>>( &settings.w0_sr )->multitoken()->default_value( std::vector<double>{ 0.3 }, "0.3" ),
            "Задать параметр(ы) веса w0 для квадратно-корневых сигма-точечных фильтров" )
        ( "alpha_sr", po::value<std::vector<double>>( &settings.alpha_sr )->multitoken()->default_value( std::vector<double>{ 0.45 }, "0.45" ), //0.687369
            "Задать параметр(ы) alpha для квадратно-корневых сигма-точечных фильтров" )
        ( "beta_sr", po::value<std::vector<double>>( &settings.beta_sr )->multitoken()->default_value( std::vector<double>{ 2.0 }, "2" ),
            "Задать параметр(ы) beta для квадратно-корневых сигма-точечных фильтров" )
        ( "kappa_sr", po::value<std::vector<double>>( &settings.kappa_sr )->multitoken()->default_value( std::vector<double>{ -2.0 }, "-2" ),
            "Задать параметр(ы) kappa для квадратно-корневых сигма-точечных фильтров" )
    ;
    po::options_description cla; // Аргументы командной строки (сommand line arguments)
    cla.add( desc );
    po::variables_map vm;
//    po::store( po::command_line_parser( argc, argv ).options( cla ).run(), vm );
    po::store( po::parse_command_line( argc, argv, cla, po::command_line_style::unix_style ^ po::command_line_style::allow_short ), vm );
    po::notify( vm );
    //------------------------------------------------------------------------------------------------------------------
    // Обработаем аргументы запуска приложения
    if( vm.count( "help" ) ) {
        std::cout << desc << std::endl;
        return EXIT_SUCCESS;
    }
//    if( vm.count( "seed" ) ) {
//        settings.Seed = vm["seed"].as<uint32_t>();
//        if( settings.Seed > 0 ) {
////            imitator.SetupSeed( settings.Seed );
//        }
//        if( settings.Seed == 0 ) {
//            std::random_device rd;
//            std::mt19937 gen( rd() );
//            std::uniform_int_distribution<> distrib( 1, INT32_MAX );
//            uint32_t randomSeed = distrib( gen );
////            imitator.SetupSeed( randomSeed );
//        }
//    }
    if( vm.count( "debug" ) ) {
        settings.Debug = true;
    }
    if( vm.count( "par" ) ) {
        settings.MatPlotParams = vm["par"].as<std::vector<double>>();
        assert( settings.MatPlotParams.size() == 6 );
    }
    if( vm.count( "prob" ) ) {
        settings.Probabilities = vm["prob"].as<std::vector<double>>();
        for( auto &p : settings.Probabilities ) {
            assert( p > 0.0 && p <= 1.0 );
        }
    }
    if( vm.count( "dt" ) ) {
        settings.DeltaT = vm["dt"].as<double>();
        assert( settings.DeltaT > 0.0 );
    }
    if( vm.count( "time" ) ) {
        settings.SimulationTime = vm["time"].as<double>();
        assert( settings.SimulationTime > 0.0 );
    }
    if( vm.count( "graph" ) ) {
        settings.Graphs_0_RMSE_1 = vm["graph"].as<uint32_t>();
        assert( settings.Graphs_0_RMSE_1 == 0 || settings.Graphs_0_RMSE_1 == 1 );
    }
    if( vm.count( "size" ) ) {
        settings.Size = vm["size"].as<std::vector<double>>();
        assert( settings.Size.size() == 2 );
    }
    if( vm.count( "sep" ) ) {
        settings.GraphSeparated = true;
    }
    if( vm.count( "eps" ) ) {
        settings.Format = "eps";
    }
    if( vm.count( "show" ) ) {
        settings.ShowGraphs = true;
    }
    if( vm.count( "mc" ) ) {
        settings.MCruns = vm["mc"].as<uint32_t>();
    }
#ifdef EKF_
    if( vm.count( "ekf" ) ) {
        settings.Filters.push_back("EKF");
    }
#endif
#ifdef SREKF_
    if( vm.count( "srekf" ) ) {
        settings.Filters.push_back("SREKF");
    }
#endif
#ifdef UKF_
    if( vm.count( "ukf" ) ) {
        settings.Filters.push_back("UKF");
    }
#endif
#ifdef SRUKF_
    if( vm.count( "srukf" ) ) {
        settings.Filters.push_back("SRUKF");
    }
    if( vm.count( "srukfb" ) ) {
        settings.Filters.push_back("SRUKFB");
    }
#endif
#ifdef CKF_
    if( vm.count( "ckf" ) ) {
        settings.Filters.push_back("CKF");
    }
#endif
#ifdef SRCKF_
    if( vm.count( "srckf" ) ) {
        settings.Filters.push_back("SRCKF");
    }
    if( vm.count( "srckfb" ) ) {
        settings.Filters.push_back("SRCKFB");
    }
#endif
#ifdef EUKF_
    if( vm.count( "eukf" ) ) {
        settings.Filters.push_back("EUKF");
    }
#endif
#ifdef SREUKF_
    if( vm.count( "sreukf" ) ) {
        settings.Filters.push_back("SREUKF");
    }
    if( vm.count( "sreukfb" ) ) {
        settings.Filters.push_back("SREUKFB");
    }
#endif
#ifdef ECKF_
    if( vm.count( "eckf" ) ) {
        settings.Filters.push_back("ECKF");
    }
#endif
#ifdef SRECKF_
    if( vm.count( "sreckf" ) ) {
        settings.Filters.push_back("SRECKF");
    }
    if( vm.count( "sreckfb" ) ) {
        settings.Filters.push_back("SRECKFB");
    }
#endif
    if( vm.count( "all" ) ) {
        settings.Filters.clear();
#ifdef EKF_
        settings.Filters.push_back("EKF");
#endif
#ifdef SREKF_
        settings.Filters.push_back("SREKF");
#endif
#ifdef UKF_
        settings.Filters.push_back("UKF");
#endif
#ifdef SRUKF_
        settings.Filters.push_back("SRUKF");                
        settings.Filters.push_back("SRUKFB");
#endif
#ifdef CKF_
        settings.Filters.push_back("CKF");
#endif
#ifdef SRCKF_
        settings.Filters.push_back("SRCKF");     
        settings.Filters.push_back("SRCKFB");
#endif
#ifdef EUKF_
        settings.Filters.push_back("EUKF");
#endif
#ifdef SREUKF_
        settings.Filters.push_back("SREUKF");
        settings.Filters.push_back("SREUKFB");
#endif
#ifdef ECKF_
        settings.Filters.push_back("ECKF");
#endif
#ifdef SRECKF_
        settings.Filters.push_back("SRECKF");     
        settings.Filters.push_back("SRECKFB");
#endif
    }

    if( vm.count( "set" ) ) {
        settings.Set = vm["set"].as<uint32_t>();
        assert( settings.Set == 0 || settings.Set == 1 );
    }
    if( settings.Set == 0 ) { // Julier
        if( vm.count( "w0" ) ) {
            settings.w0 = vm["w0"].as<std::vector<double>>();
        }

        if( vm.count( "w0_sr" ) ) {
            settings.w0_sr = vm["w0_sr"].as<std::vector<double>>();
        }
    } else if( settings.Set == 1 ) { // Merwe
        if( vm.count( "alpha" ) ) {
            settings.alpha = vm["alpha"].as<std::vector<double>>();
        }
        if( vm.count( "beta" ) ) {
            settings.beta = vm["beta"].as<std::vector<double>>();
        }
        if( vm.count( "kappa" ) ) {
            settings.kappa = vm["kappa"].as<std::vector<double>>();
        }

        if( vm.count( "alpha_sr" ) ) {
            settings.alpha_sr = vm["alpha_sr"].as<std::vector<double>>();
        }
        if( vm.count( "beta_sr" ) ) {
            settings.beta_sr = vm["beta_sr"].as<std::vector<double>>();
        }
        if( vm.count( "kappa_sr" ) ) {
            settings.kappa_sr = vm["kappa_sr"].as<std::vector<double>>();
        }
    } else {
        assert( false );
    }
    //------------------------------------------------------------------------------------------------------------------
    // Запуск тестов
    if( vm.count( "rmse" ) ) {
        test.Run_RMSE_EKF_SRUKF_SREUKF( settings );
    } else {
        test.RunMain( settings );
    }
    return 0;
}

/// \}
