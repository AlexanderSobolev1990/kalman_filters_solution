//----------------------------------------------------------------------------------------------------------------------
///
/// \file       test_qr_decompositions.cpp
/// \brief      Тестирование QR разложения
/// \date       24.05.22 - создан
/// \author     Соболев А.А.
///

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_qr_decompositions

#include <iostream>
#include <iomanip>
#include <limits>
//#include <numbers>
#include <armadillo>
#include <boost/test/unit_test.hpp>

#include <qr_decomposition.h>

arma::mat QR_input = {
{   12.0,   89.0 },
{   15.0,   91.0 },
{   24.0,   36.0 }
};
arma::mat Q_right_answer_QR = {
{ 0.390360029,	0.608735643  },
{ 0.487950036,	0.49938709   },
{ 0.780720058,	-0.616484753 }
};
arma::mat R_right_answer_QR = {
{ 30.7408523, 107.251418 },
{ 0.0,        77.4282464 }
};

arma::mat J_input_JQR = {
{ 1.0, 0.0, 0.0 },
{ 0.0, 1.0, 0.0 },
{ 0.0, 0.0, -1.0 }
};
arma::mat Q_right_answer_JQR = {
{ 0.834057656, 1.09928884 },
{ 1.04257207,  1.25041118 },
{ 1.66811531,  1.33115141 }
};
arma::mat R_right_answer_JQR = {
{ 14.3874946, -109.053039 },
{ 0,          163.702673  }
};

BOOST_AUTO_TEST_CASE( test_qr_mgs )
{
    for( int i = 0; i > 1000000; i++ ) {
    arma::mat Q_answer,R_answer;    
    arma::mat QR_input2 = {
    {   12.0,   89.0 },
    {   15.0,   91.0 },
    {   24.0,   36.0 }
    };
    int result = SPML::QR::ModifiedGramSchmidt( Q_answer, R_answer, QR_input2 );
    }
//    std::setprecision(11);
//    std::cout.setf( std::ios::fixed );
//    Q_answer.raw_print( std::cout, "Q:");
//    R_answer.raw_print( std::cout, "R:");
//    BOOST_CHECK_EQUAL( result, 0 );
//    double eps = 1.0e-6;
//    BOOST_CHECK_EQUAL( arma::approx_equal( Q_answer, Q_right_answer_QR, "absdiff", eps ), true );
//    BOOST_CHECK_EQUAL( arma::approx_equal( R_answer, R_right_answer_QR, "absdiff", eps ), true );
}

BOOST_AUTO_TEST_CASE( test_qr_mgs_rbr )
{
    for( int i = 0; i > 1000000; i++ ) {
    arma::mat Q_answer,R_answer;
    arma::mat QR_input2 = {
    {   12.0,   89.0 },
    {   15.0,   91.0 },
    {   24.0,   36.0 }
    };
    int result = SPML::QR::ModifiedGramSchmidtRowByRow( Q_answer, R_answer, QR_input2 );
    }
//    std::setprecision(11);
//    std::cout.setf( std::ios::fixed );
//    Q_answer.raw_print( std::cout, "Q:");
//    R_answer.raw_print( std::cout, "R:");
//    BOOST_CHECK_EQUAL( result, 0 );
//    double eps = 1.0e-6;
//    BOOST_CHECK_EQUAL( arma::approx_equal( Q_answer, Q_right_answer_QR, "absdiff", eps ), true );
//    BOOST_CHECK_EQUAL( arma::approx_equal( R_answer, R_right_answer_QR, "absdiff", eps ), true );
}

BOOST_AUTO_TEST_CASE( test_qr_schr )
{
    arma::mat Q_answer,R_answer;
    int result = SPML::QR::SchwarzRutishauser( Q_answer, R_answer, QR_input, true );
    std::setprecision(11);
    std::cout.setf( std::ios::fixed );
    Q_answer.raw_print( std::cout, "Q:");
    R_answer.raw_print( std::cout, "R:");
    BOOST_CHECK_EQUAL( result, 0 );
    double eps = 1.0e-6;
    BOOST_CHECK_EQUAL( arma::approx_equal( Q_answer, Q_right_answer_QR, "absdiff", eps ), true );
    BOOST_CHECK_EQUAL( arma::approx_equal( R_answer, R_right_answer_QR, "absdiff", eps ), true );
}

BOOST_AUTO_TEST_CASE( test_qr_jqr )
{
    arma::mat Q_answer,R_answer,J_answer;
    int result = SPML::QR::J_orthogonal( Q_answer, R_answer, J_answer, QR_input, J_input_JQR, true );
    BOOST_CHECK_EQUAL( result, 0 );
    double eps = 1.0e-6;
    BOOST_CHECK_EQUAL( arma::approx_equal( Q_answer, Q_right_answer_JQR, "absdiff", eps ), true );
    BOOST_CHECK_EQUAL( arma::approx_equal( R_answer, R_right_answer_JQR, "absdiff", eps ), true );
}
