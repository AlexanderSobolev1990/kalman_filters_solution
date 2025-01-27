//------------------------------------------------------------------------------
///
/// \file       qr_decomposition.cpp
/// \brief      QR-разложение матриц
/// \date       01.04.21 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#include <qr_decomposition.h>

namespace SPML /// Специальная библиотека программных модулей (СБ ПМ)
{

arma::vec sgn2( arma::vec v )
{
    arma::vec res( v.n_elem );
    for( auto i = 0; i < v.n_elem; i++ ) {
        if( v[i] >= 0 ) {
            res[i] = 1;
        } else {
            res[i] = -1;
        }
    }
    return res;
}

//------------------------------------------------------------------------------    
namespace QR /// QR-разложение матриц
{

void MGS( arma::mat &Q, arma::mat &R, const arma::mat &A )
{
    int m = A.n_rows;
    int n = A.n_cols;
    Q = arma::mat( m, n, arma::fill::zeros );
    R = arma::mat( n, n, arma::fill::zeros );
    for( int k = 0; k < n; k++ ) {
        Q.col(k) = A.col(k);
        for( int i = 0; i < k; i++ ) {
            R(i, k) = arma::dot( arma::trans( Q.col(i) ), Q.col(k) );
            Q.col(k) = Q.col(k) - R(i, k) * Q.col(i);
        }
        R(k, k) = ( arma::norm( Q.col(k) ) );
        Q.col(k) = Q.col(k) / R(k, k);
    }    
}

//------------------------------------------------------------------------------
void MGS_1( arma::mat &R, const arma::mat &A )
{
    arma::mat Q;
    MGS( Q, R, A );

    // int m = A.n_rows;
    // int n = A.n_cols;
    // arma::mat Q( m, n, arma::fill::zeros );
    // R = arma::mat( n, n, arma::fill::zeros );
    // for( int k = 0; k < n; k++ ) {
    //     Q.col(k) = A.col(k);
    //     for( int i = 0; i < k; i++ ) {
    //         R(i, k) = arma::dot( arma::trans( Q.col(i) ), Q.col(k) );
    //         Q.col(k) = Q.col(k) - R(i, k) * Q.col(i);
    //     }
    //     R(k, k) = ( arma::norm( Q.col(k) ) );
    //     Q.col(k) = Q.col(k) / R(k, k);
    // }    
}

//------------------------------------------------------------------------------
int SchwarzRutishauser( arma::mat &Q, arma::mat &R, const arma::mat &A, 
    bool econ, double tol )
{
    int m = A.n_rows;
    int n = A.n_cols;
    Q = arma::mat( m, n, arma::fill::zeros );
    R = arma::mat( n, n, arma::fill::zeros );
    int orig, flag;
    double s = 0.0;
    double t = 0.0;
    int max_iter = 20;
    int iter = 0;
    for( int k = 0; k < n; k++ ) {
        for( int j = 0; j < m; j++ ) {
            Q(j, k) = A(j, k);
        }
        orig = 1;
        flag = 1;
        while( flag ) {
            iter++;
            if( iter > max_iter ) {
                break;
            }
            t = 0.0;
            for( int i = 0; i < k; i++ ) {
                s = 0.0;
                for( int j = 0; j < m; j++ ) {
                    s += ( Q(j, i) * Q(j, k) );
                }
                if( orig ) {
                    R(i, k) = s;
                } else {
                    R(i, k) = R(i, k) + s;
                }
                t += ( s * s );
                for( int j = 0; j < m; j++ ) {
                    Q(j, k) -= ( s * Q(j, i) );
                }
            } // end for i
            s = 0.0;
            for( int j = 0; j < m; j++ ) {
                s += ( Q(j, k) * Q(j, k) );
            }
            t += s;
            flag = 0;
//            if( ( s < ( t / 100.0 ) ) || ( std::abs( s - ( t / 100.0 ) ) < DBL_EPSILON || ( t * DBL_EPSILON == 0 ) )  ) { // s <= t // WORKED!
            if( ( s < ( t / 100.0 ) ) || ( std::abs( s - ( t / 100.0 ) ) < tol || ( std::abs( t ) < tol ) )  ) { // s <= t
                orig = 0;
                if( std::abs( s ) < tol ) {
//                if( ( s * DBL_EPSILON ) == 0.0 ) {
                    s = 0.0;
                } else {
                    flag = 1;
                }
            } // end if
        }
        s = std::sqrt( s );
        R(k, k) = s;
//        if( std::abs( s ) > DBL_EPSILON ) {
        if( std::abs( s ) > tol ) {
            s = 1.0 / s;
        }
        for( int j = 0; j < m; j++ ) {
            Q(j, k) *= s;
        }
//        std::cout << iter << std::endl;
    }
    if( econ && ( m > n ) ) {
        R = R( arma::span( 0, ( n - 1 ) ), arma::span( 0, ( n - 1 ) ) );
        Q = Q( arma::span(), arma::span( 0, ( n - 1 ) ));
    }
    return 0;
}

}

//------------------------------------------------------------------------------    
namespace LQ /// LQ-разложение матриц
{

void MGS( arma::mat &Q, arma::mat &L, const arma::mat &A )
{
    int m = A.n_rows;
    int n = A.n_cols;    
    Q = arma::mat( m, n, arma::fill::zeros );
    L = arma::mat( m, m, arma::fill::zeros );
    for( int k = 0; k < m; k++ ) {    
        Q.row(k) = A.row(k);
        for( int i = 0; i < k; i++ ) {        
            L(k, i) = arma::dot( Q.row(i), arma::trans( Q.row(k) ) );
            Q.row(k) = Q.row(k) - L(k, i) * Q.row(i);
        }
        L(k, k) = ( arma::norm( Q.row(k) ) );
        Q.row(k) = Q.row(k) / L(k, k);
    }    
}

//------------------------------------------------------------------------------
void MGS_1( arma::mat &L, const arma::mat &A )
{
    arma::mat Q;
    MGS( Q, L, A );
    // int m = A.n_rows;
    // int n = A.n_cols;    
    // arma::mat Q = arma::mat( m, n, arma::fill::zeros );
    // L = arma::mat( m, m, arma::fill::zeros );
    // for( int k = 0; k < m; k++ ) {    
    //     Q.row(k) = A.row(k);
    //     for( int i = 0; i < k; i++ ) {        
    //         L(k, i) = arma::dot( Q.row(i), arma::trans( Q.row(k) ) );
    //         Q.row(k) = Q.row(k) - L(k, i) * Q.row(i);
    //     }
    //     L(k, k) = ( arma::norm( Q.row(k) ) );
    //     Q.row(k) = Q.row(k) / L(k, k);
    // }    
}

}

//------------------------------------------------------------------------------
namespace JQR /// JQR-разложение матриц
{

int JQR( arma::mat &Q, arma::mat &R, arma::vec &Jp, const arma::mat &A, 
    const arma::vec &J, bool econ, double tol )
{
    int m = A.n_rows;
    int n = A.n_cols;

    // arma::vec J_ = J.diag();
    arma::vec J_ = J;

    bool od_method = true; // Use orthogonal-diagonal method for real matrices

    // Initialize matrices
    R = A; // upper-triangular matrix
    Q = arma::mat( m, m, arma::fill::eye ); // J-orthogonal matrix
    
    arma::ivec perm( m ); // row permutations
    for( int i = 0; i < m; i++ ) {
        perm( i ) = i;
    }

    int r = 0; // Column index
    int c = 0; // Row index

    double t = 0.0;
    double z = 0.0;

    int r1 = 0;
    int r2 = 0;
    int retry = 0;

    double a1 = 0.0;
    double a2 = 0.0;
    double a3 = 0.0;
    double da = 0.0;

    std::deque<int> r2_saved;

    arma::rowvec x;
    arma::rowvec y;
    arma::rowvec s;
    arma::rowvec q;

    int min_n_m = std::min( n, m );

    while( c < min_n_m ) {
//        c++;
//        r++;
        r1 = r;
        r2 = r1;
        retry = 0;

        r2_saved.clear();

        while( ( r2 < ( m - 1 ) ) || ( retry > 0 ) ) { // Cancel entries in the column
            if( r2 < ( m - 1 ) ) { // Proceed to cancel entries
                r2++;
            } else { // Retry to cancel saved entries
//                r2 = r2_saved[retry];
                r2 = r2_saved[retry-1];
                // retry--;
                if( retry > 0 ) {
                    retry--;
                    // std::cout << "retry > 2 JQR!" << std::endl;
                    // assert( false );
                    // return 1;
                }
                // retry--;
                std::cout << "Retry JQR!" << std::endl;
//                assert( false );
//                return 1;
            }
            a1 = R( r1, c );
            a2 = R( r2, c );
            if( std::abs( a2 ) <= ( tol * std::abs( a1 ) ) ) { // ORIG
//            if( ( std::abs( a2 ) < ( tol * std::abs( a1 ) ) ) ||
//                ( std::abs( std::abs( a2 ) - ( tol * std::abs( a1 ) ) ) < tol ) ) { // MY FIX 2
//            if( std::abs( a2 ) < tol ) { // MY FIX
                t = sgn2( a1, tol );
                R.row( r1 ) = t * R.row( r1 );
                R.row( r2 ) = t * R.row( r2 );
                R( r2, c ) = 0.0;
                Q.row( r1 ) = t * Q.row( r1 );
                Q.row( r2 ) = t * Q.row( r2 );
            } else { // Non-zero a2
                if( J_( r1 ) == J_( r2 ) ) { // Same sign2, apply orthogonal rotation

                    // Apply Givens rotation using mixed-downdating (from [3])
                    if( std::abs( a1 ) <= ( tol * std::abs( a2 ) ) ) { // ORIG
//                    if( ( std::abs( a1 ) < ( tol * std::abs( a2 ) ) ) ||
//                        ( std::abs( std::abs( a1 ) - ( tol * std::abs( a2 ) ) ) < tol ) ) { // MY FIX 2
//                    if( std::abs( a1 ) < tol ) { // MY FIX
                        t = sgn2( a2, tol );
                        s = R.row( r1 );
                        R.row( r1 ) = t * R.row( r2 );
                        R.row( r2 ) = ( -t ) * s;
                        R( r2, c ) = 0.0;
                        q = Q.row( r1 );
                        Q.row( r1 ) = t * Q.row( r2 );
                        Q.row( r2 ) = ( -t ) * q;
                    } else if( std::abs( a1 ) > std::abs( a2 ) ) {
                        t = a2 / a1;
                        z = sgn2( a1, tol ) * std::sqrt( 1.0 + ( t * t ) );
                        R.row( r1 ) = ( R.row( r1 ) + t * ( R.row( r2 ) ) ) / z;
                        R.row( r2 ) = ( ( -t ) * R.row( r1 ) ) + ( z * R.row( r2 ) );
                        R( r2, c ) = 0.0;
                        Q.row( r1 ) = ( Q.row( r1 ) + ( t * Q.row( r2 ) ) ) / z;
                        Q.row( r2 ) = ( ( -t ) * Q.row( r1 ) ) + ( z * Q.row( r2 ) );
                    } else {
                        t = a1 / a2;
                        z = sgn2( a2, tol ) * std::sqrt( 1.0 + ( t * t ) );
                        R.row( r1 ) = ( ( t * R.row( r1 ) ) + R.row( r2 ) ) / z;
                        R.row( r2 ) = ( -R.row( r1 ) + ( z * R.row( r2 ) ) ) / t;
                        R( r2, c ) = 0.0;
                        Q.row( r1 ) = ( ( t * Q.row( r1 ) ) + Q.row( r2 ) ) / z;
                        Q.row( r2 ) = ( -Q.row( r1 ) + ( z * Q.row( r2 ) ) ) / t;
                    }
                } else { // Different sign2, apply hyperbolic rotation..
                    da = std::abs( a1 ) - std::abs( a2 );
                    if( std::abs( da ) > ( tol * ( std::abs( a1 ) + std::abs( a2 ) ) ) ) {
                        if( da < 0.0 ) { // Permute rows if necessary
                            Q.swap_rows( r2, r1 );
                            Q.swap_cols( r2, r1 );
                            R.swap_rows( r2, r1 );
                            J_.swap_rows( r2, r1 );
                            perm.swap_rows( r2, r1 );
                            a3 = a1;
                            a1 = a2;
                            a2 = a3;
                        }
                        if( std::abs( a2 ) <= ( tol * std::abs( a1 ) ) ) { // ORIG
//                        if( ( std::abs( a2 ) < ( tol * std::abs( a1 ) ) ) ||
//                            ( std::abs( std::abs( a2 ) - ( tol * std::abs( a1 ) ) ) < tol ) ) { // MY FIX 2
//                        if( std::abs( a2 ) < tol ) { // MY FIX
                            t = sgn2( a1, tol );
                            R.row( r1 ) = t * R.row( r1 );
                            R.row( r2 ) = t * R.row( r2 );
                            R( r2, c ) = 0.0;
                            Q.row( r1 ) = t * Q.row( r1 );
                            Q.row( r2 ) = t * Q.row( r2 );
                        } else { // Non-zero a2
                            if( od_method ) {
                                // Apply hyperbolic rotation using orthogonal-diagonal
                                // method (forward numerical stable) (from [2])
                                double a = 0.5 * sgn2( a1, tol );
                                t = a * std::sqrt( ( a1 + a2 ) / ( a1 - a2 ) );
                                z = a * std::sqrt( ( a1 - a2 ) / ( a1 + a2 ) );
                                x = R.row( r1 ) - R.row( r2 );
                                y = R.row( r1 ) + R.row( r2 );
                                x = t * x;
                                y = z * y;
                                R.row( r1 ) = ( x + y );
                                R.row( r2 ) = ( y - x );
                                R( r2, c ) = 0.0;
                                x = Q.row( r1 ) - Q.row( r2 );
                                y = Q.row( r1 ) + Q.row( r2 );
                                x = t * x;
                                y = z * y;
                                Q.row( r1 ) = ( x + y );
                                Q.row( r2 ) = ( y - x );
                            } else {
                                // Apply hyperbolic rotation using mixed-downdating
                                t = a2 / a1;
                                z = sgn2( a1, tol ) * std::sqrt( 1.0 - ( t * t ) );
                                R.row( r1 ) = ( R.row( r1 ) - ( t * R.row( r2 ) ) ) / z;
                                R.row( r2 ) = ( ( -t ) * R.row( r1 ) ) + ( z * R.row( r2 ) );
                                R( r2, c ) = 0.0;
                                Q.row( r1 ) = ( Q.row( r1 ) - ( t * Q.row( r2 ) ) ) / z;
                                Q.row( r2 ) = ( ( -t ) * Q.row( r1 ) ) + ( z * Q.row( r2 ) );
                            }
                        }
                    } else { // No hyperbolic rotation
                        if( ( r2 < ( m - 1 ) ) || ( retry > 0 ) ) { // Still entries remaining
                            // Cancelling all remaining entries first
                            retry++;
                            r2_saved.push_front( r2 );
                        } else { // No hyperbolic rotation
                            std::cout << "Error: JQR failed to apply hyperbolic rotation. Check tolerance or method." << std::endl;
                            assert( false );
                            return 1;
                        }
                    }
                } // end if-else same sign2
            }
        } // while Cancel entries in the column

        c++;
        r++;
    } // while

    // Make last diagonal entry positive
    if( R( r-1, c-1 ) < 0.0 ) {
        Q.row( r-1 ) = -Q.row( r-1 );
        R.row( r-1 ) = -R.row( r-1 );
    }

    // Post-processing    
    Jp = J_;    
    arma::mat tmp = arma::diagmat( J_ ) * Q.t() * arma::diagmat( J_ );
    for( int i = 0; i < m; i++ ) {
        Q.row( perm( i ) ) = tmp.row( i );
    }
    if( econ && ( m > n ) ) {
        R = R( arma::span( 0, ( n - 1 ) ), arma::span( 0, ( n - 1 ) ) );
        Q = Q( arma::span(), arma::span( 0, ( n - 1 ) ) );        
        Jp = Jp( arma::span( 0, ( n - 1 ) ) );
    }
    return 0;
}

//------------------------------------------------------------------------------
int JQR_1( arma::mat &R, const arma::mat &A, const arma::vec &J )
{    
    arma::mat Q;
    arma::vec Jp;
    int res = JQR( Q, R, Jp, A, J, true );
    return res;
}    

}
//------------------------------------------------------------------------------
namespace JLQ /// JLQ-разложение матриц
{

int JLQ( arma::mat &Q, arma::mat &L, arma::vec &Jp, const arma::mat &A, 
    const arma::vec &J, bool econ, double tol )
{
    int res = JQR::JQR( Q, L, Jp, A.t(), J, econ, tol );
    L = L.t();
    return res;
}

//------------------------------------------------------------------------------
int JLQ_1( arma::mat &L, const arma::mat &A, const arma::vec &J )
{   
    int res = JQR::JQR_1( L, A.t(), J );
    L = L.t();
    return res;
}  

}

}
/// \}
