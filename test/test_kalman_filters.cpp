//----------------------------------------------------------------------------------------------------------------------
///
/// \file       test_kalman_filters.cpp
/// \brief      Тестирование задач о назначениях
/// \date       20.05.22 - создан
/// \author     Соболев А.А.
///

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_kalman_filters

#include <iostream>
#include <armadillo>
#include <boost/test/unit_test.hpp>
#include <vector>
#include <map>
#include <matplotlibcpp.h>

//#include <compare.h>
#include <kalman_filter_linear.h>

const bool print_to_console = false;//true;//

const size_t dimX = 2;
const size_t dimY = 3;
/*
BOOST_AUTO_TEST_CASE( test_LKF_insert_delete_map )
{
    CKalmanFiltersCompare c;
    KalmanFilters::CKalmanLKF<dimX, dimY> lkf1, lkf2, lkf3;
    arma::mat A1( dimX, dimX ), A2( dimX, dimX ), A3( dimX, dimX );
    A1.fill( 1.0 );
    A2.fill( 2.0 );
    A3.fill( 3.0 );
    lkf1.SetEstimateCovarianceMatrixP( A1 );
    lkf2.SetEstimateCovarianceMatrixP( A2 );
    lkf3.SetEstimateCovarianceMatrixP( A3 );

    std::map<int, KalmanFilters::CKalmanLKF<dimX, dimY> > filters;
    filters.insert( std::make_pair( 1, lkf1 ) );
    filters.insert( std::make_pair( 2, lkf2 ) );
    filters.insert( std::make_pair( 3, lkf3 ) );

    if( print_to_console ) {
        for( auto &f : filters ) {
            std::cout << f.first << std::endl;
            (f.second).GetEstimatedCovarianceMatrixP().print();
        }
    }

    for(auto it = filters.begin(); it != filters.end(); ) {
        if( it->first == 2 ) {
            it = filters.erase( it );
        } else {
            ++it;
        }
    }

    if( print_to_console ) {
        for( auto &f : filters ) {
            std::cout << f.first << std::endl;
            (f.second).GetEstimatedCovarianceMatrixP().print();
        }
    }

    BOOST_CHECK_EQUAL( filters.count( 2 ), 0 );
}

BOOST_AUTO_TEST_CASE( test_UKF_insert_delete_map )
{
    CKalmanFiltersCompare c;
    KalmanFilters::CKalmanUKF<dimX, dimY> lkf1, lkf2, lkf3;
    arma::mat A1( dimX, dimX ), A2( dimX, dimX ), A3( dimX, dimX );
    A1.fill( 1.0 );
    A2.fill( 2.0 );
    A3.fill( 3.0 );
    lkf1.SetEstimateCovarianceMatrixP( A1 );
    lkf2.SetEstimateCovarianceMatrixP( A2 );
    lkf3.SetEstimateCovarianceMatrixP( A3 );

    std::map<int, KalmanFilters::CKalmanUKF<dimX, dimY> > filters;
    filters.insert( std::make_pair( 1, lkf1 ) );
    filters.insert( std::make_pair( 2, lkf2 ) );
    filters.insert( std::make_pair( 3, lkf3 ) );

    if( print_to_console ) {
        for( auto &f : filters ) {
            std::cout << f.first << std::endl;
            (f.second).GetEstimatedCovarianceMatrixP().print();
        }
    }

    for(auto it = filters.begin(); it != filters.end(); ) {
        if( it->first == 2 ) {
            it = filters.erase( it );
        } else {
            ++it;
        }
    }

    if( print_to_console ) {
        for( auto &f : filters ) {
            std::cout << f.first << std::endl;
            (f.second).GetEstimatedCovarianceMatrixP().print();
        }
    }

    BOOST_CHECK_EQUAL( filters.count( 2 ), 0 );

}

BOOST_AUTO_TEST_CASE( test_CKalmanSREUKFB_insert_delete_map )
{
    CKalmanFiltersCompare c;
    KalmanFilters::CKalmanSREUKFB<dimX, dimY> lkf1, lkf2, lkf3;
    arma::mat A1( dimX, dimX ), A2( dimX, dimX ), A3( dimX, dimX );
    A1.fill( 1.0 );
    A2.fill( 2.0 );
    A3.fill( 3.0 );
    lkf1.SetEstimateCovarianceMatrixP( A1 );
    lkf2.SetEstimateCovarianceMatrixP( A2 );
    lkf3.SetEstimateCovarianceMatrixP( A3 );

    std::map<int, KalmanFilters::CKalmanUKF<dimX, dimY> > filters;
    filters.insert( std::make_pair( 1, lkf1 ) );
    filters.insert( std::make_pair( 2, lkf2 ) );
    filters.insert( std::make_pair( 3, lkf3 ) );

    if( print_to_console ) {
        for( auto &f : filters ) {
            std::cout << f.first << std::endl;
            (f.second).GetEstimatedCovarianceMatrixP().print();
        }
    }

    for(auto it = filters.begin(); it != filters.end(); ) {
        if( it->first == 2 ) {
            it = filters.erase( it );
        } else {
            ++it;
        }
    }

    if( print_to_console ) {
        for( auto &f : filters ) {
            std::cout << f.first << std::endl;
            (f.second).GetEstimatedCovarianceMatrixP().print();
        }
    }

    BOOST_CHECK_EQUAL( filters.count( 2 ), 0 );

}
*/
BOOST_AUTO_TEST_CASE( test_LKF_insert_delete_vector )
{
//    CKalmanFiltersCompare c;
    KalmanFilters::CKalmanLKF<dimX, dimY> lkf1, lkf2, lkf3;
    arma::mat A1( dimX, dimX ), A2( dimX, dimX ), A3( dimX, dimX );
    A1.fill( 1.0 );
    A2.fill( 2.0 );
    A3.fill( 3.0 );
    lkf1.SetEstimateCovarianceMatrixP( A1 );
    lkf2.SetEstimateCovarianceMatrixP( A2 );
    lkf3.SetEstimateCovarianceMatrixP( A3 );

    std::vector< KalmanFilters::CKalmanLKF<dimX, dimY> > filters;
    filters.push_back( lkf1 );
    filters.push_back( lkf2 );
    filters.push_back( lkf3 );

    if( print_to_console ) {
        for( auto &f : filters ) {
            f.GetEstimatedCovarianceMatrixP().print();
        }
    }

    for( auto it = filters.begin(); it != filters.end(); ) {
        if( ( ( *it ).GetEstimatedCovarianceMatrixP() ).at(0,0) == 2 ) {
            it = filters.erase( it );
        } else {
            ++it;
        }
    }

    if( print_to_console ) {
        for( auto &f : filters ) {
            f.GetEstimatedCovarianceMatrixP().print();
        }
    }

    BOOST_CHECK_EQUAL( filters.size(), 2 );
}
/*
BOOST_AUTO_TEST_CASE( test_UKF_insert_delete_vector )
{
    CKalmanFiltersCompare c;
    KalmanFilters::CKalmanUKF<dimX, dimY> lkf1, lkf2, lkf3;
    arma::mat A1( dimX, dimX ), A2( dimX, dimX ), A3( dimX, dimX );
    A1.fill( 1.0 );
    A2.fill( 2.0 );
    A3.fill( 3.0 );
    lkf1.SetEstimateCovarianceMatrixP( A1 );
    lkf2.SetEstimateCovarianceMatrixP( A2 );
    lkf3.SetEstimateCovarianceMatrixP( A3 );

    std::vector< KalmanFilters::CKalmanUKF<dimX, dimY> > filters;
    filters.push_back( lkf1 );
    filters.push_back( lkf2 );
    filters.push_back( lkf3 );

    if( print_to_console ) {
        for( auto &f : filters ) {
            f.GetEstimatedCovarianceMatrixP().print();
        }
    }

    for( auto it = filters.begin(); it != filters.end(); ) {
        if( ( ( *it ).GetEstimatedCovarianceMatrixP() ).at(0,0) == 2 ) {
            it = filters.erase( it );
        } else {
            ++it;
        }
    }

    if( print_to_console ) {
        for( auto &f : filters ) {
            f.GetEstimatedCovarianceMatrixP().print();
        }
    }

    BOOST_CHECK_EQUAL( filters.size(), 2 );

}

BOOST_AUTO_TEST_CASE( test_SREUKFB_insert_delete_vector )
{
    CKalmanFiltersCompare c;
    KalmanFilters::CKalmanSREUKFB<dimX, dimY> lkf1, lkf2, lkf3;
    arma::mat A1( dimX, dimX ), A2( dimX, dimX ), A3( dimX, dimX );
    A1.fill( 1.0 );
    A2.fill( 2.0 );
    A3.fill( 3.0 );
    lkf1.SetEstimateCovarianceMatrixP( A1 );
    lkf2.SetEstimateCovarianceMatrixP( A2 );
    lkf3.SetEstimateCovarianceMatrixP( A3 );

    std::vector< KalmanFilters::CKalmanLKF<dimX, dimY> > filters;
    filters.push_back( lkf1 );
    filters.push_back( lkf2 );
    filters.push_back( lkf3 );

    if( print_to_console ) {
        for( auto &f : filters ) {
            f.GetEstimatedCovarianceMatrixP().print();
        }
    }

    for( auto it = filters.begin(); it != filters.end(); ) {
        if( ( ( *it ).GetEstimatedCovarianceMatrixP() ).at(0,0) == 2 ) {
            it = filters.erase( it );
        } else {
            ++it;
        }
    }

    if( print_to_console ) {
        for( auto &f : filters ) {
            f.GetEstimatedCovarianceMatrixP().print();
        }
    }

    BOOST_CHECK_EQUAL( filters.size(), 2 );
}





*/
