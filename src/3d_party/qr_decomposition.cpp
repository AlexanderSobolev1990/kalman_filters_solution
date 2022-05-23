//----------------------------------------------------------------------------------------------------------------------
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
namespace QR /// QR-разложение матриц
{
//----------------------------------------------------------------------------------------------------------------------
int J_orthogonal_OLD( arma::mat &Q, arma::mat &R, arma::mat &J1, const arma::mat &A, const arma::mat &J0, bool econ, double tol )
{
    int m = A.n_rows;
    int n = A.n_cols;

    arma::vec J = J0.diag();

    bool od_method = true; // Use orthogonal-diagonal method for real matrices

    // Initialize matrices
    R = arma::mat( A ); // upper-triangular matrix
    Q = arma::mat( m, m, arma::fill::eye ); // J-orthogonal matrix

    arma::ivec perm( m ); // row permutations
    for( int i = 0; i < m; i++ ) {
        perm( i ) = i + 1;
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
        c++;
        r++;
        r1 = r;
        r2 = r1;
        retry = 0;

        r2_saved.clear();

        while( ( r2 < m ) || ( retry > 0 ) ) { // Cancel entries in the column
            if( r2 < m ) { // Proceed to cancel entries
                r2++;
            } else { // Retry to cancel saved entries
//                r2 = r2_saved[retry];
                r2 = r2_saved[retry-1];

                if( retry > 0 ) {
                    std::cout << "retry > 2 JQR!" << std::endl;
                    return 1;
                }
//                retry--;
                std::cout << "Retry JQR!" << std::endl;
//                assert( false ); // ???
//                return 1;
            }
            a1 = R( r1-1, c-1 );
            a2 = R( r2-1, c-1 );
            if( std::abs( a2 ) <= ( tol * std::abs( a1 ) ) ) { // ORIG
//            if( ( std::abs( a2 ) < ( tol * std::abs( a1 ) ) ) ||
//                ( std::abs( std::abs( a2 ) - ( tol * std::abs( a1 ) ) ) < tol ) ) { // MY FIX 2
//            if( std::abs( a2 ) < tol ) { // MY FIX
                t = sgn2( a1, tol );
                R.row( r1-1 ) = t * R.row( r1-1 );
                R.row( r2-1 ) = t * R.row( r2-1 );
                R( r2-1, c-1 ) = 0.0;
                Q.row( r1-1 ) = t * Q.row( r1-1 );
                Q.row( r2-1 ) = t * Q.row( r2-1 );
            } else { // Non-zero a2
                if( J( r1-1 ) == J( r2-1 ) ) { // Same sign2, apply orthogonal rotation
                    // Apply Givens rotation using mixed-downdating (from [3])
                    if( std::abs( a1 ) <= ( tol * std::abs( a2 ) ) ) { // ORIG
//                    if( ( std::abs( a1 ) < ( tol * std::abs( a2 ) ) ) ||
//                        ( std::abs( std::abs( a1 ) - ( tol * std::abs( a2 ) ) ) < tol ) ) { // MY FIX 2
//                    if( std::abs( a1 ) < tol ) { // MY FIX
                        t = sgn2( a2, tol );
                        s = R.row( r1-1 );
                        R.row( r1-1 ) = t * R.row( r2-1 );
                        R.row( r2-1 ) = ( -t ) * s;
                        R( r2-1, c-1 ) = 0.0;
                        q = Q.row( r1-1 );
                        Q.row( r1-1 ) = t * Q.row( r2-1 );
                        Q.row( r2-1 ) = ( -t ) * q;
                    } else if( std::abs( a1 ) > std::abs( a2 ) ) {
                        t = a2 / a1;
                        z = sgn2( a1, tol ) * std::sqrt( 1.0 + ( t * t ) );
                        R.row( r1-1 ) = ( R.row( r1-1 ) + t * ( R.row( r2-1 ) ) ) / z;
                        R.row( r2-1 ) = ( ( -t ) * R.row( r1-1 ) ) + ( z * R.row( r2-1 ) );
                        R( r2-1, c-1 ) = 0.0;
                        Q.row( r1-1 ) = ( Q.row( r1-1 ) + ( t * Q.row( r2-1 ) ) ) / z;
                        Q.row( r2-1 ) = ( ( -t ) * Q.row( r1-1 ) ) + ( z * Q.row( r2-1 ) );
                    } else {
                        t = a1 / a2;
                        z = sgn2( a2, tol ) * std::sqrt( 1.0 + ( t * t ) );
                        R.row( r1-1 ) = ( ( t * R.row( r1-1 ) ) + R.row( r2-1 ) ) / z;
                        R.row( r2-1 ) = ( -R.row( r1-1 ) + ( z * R.row( r2-1 ) ) ) / t;
                        R( r2-1, c-1 ) = 0.0;
                        Q.row( r1-1 ) = ( ( t * Q.row( r1-1 ) ) + Q.row( r2-1 ) ) / z;
                        Q.row( r2-1 ) = ( -Q.row( r1-1 ) + ( z * Q.row( r2-1 ) ) ) / t;
                    }
                } else { // Different sign2, apply hyperbolic rotation..
                    da = std::abs( a1 ) - std::abs( a2 );
                    if( std::abs( da ) > ( tol * ( std::abs( a1 ) + std::abs( a2 ) ) ) ) {
                        if( da < 0.0 ) { // Permute rows if necessary
                            Q.swap_rows( r2-1, r1-1 );
                            Q.swap_cols( r2-1, r1-1 );
                            R.swap_rows( r2-1, r1-1 );
                            J.swap_rows( r2-1, r1-1 );
                            perm.swap_rows( r2-1, r1-1 );
                            a3 = a1;
                            a1 = a2;
                            a2 = a3;
                        }
                        if( std::abs( a2 ) <= ( tol * std::abs( a1 ) ) ) { // ORIG
//                        if( ( std::abs( a2 ) < ( tol * std::abs( a1 ) ) ) ||
//                            ( std::abs( std::abs( a2 ) - ( tol * std::abs( a1 ) ) ) < tol ) ) { // MY FIX 2
//                        if( std::abs( a2 ) < tol ) { // MY FIX
                            t = sgn2( a1, tol );
                            R.row( r1-1 ) = t * R.row( r1-1 );
                            R.row( r2-1 ) = t * R.row( r2-1 );
                            R( r2-1, c-1 ) = 0.0;
                            Q.row( r1-1 ) = t * Q.row( r1-1 );
                            Q.row( r2-1 ) = t * Q.row( r2-1 );
                        } else { // Non-zero a2
                            if( od_method ) {
                                // Apply hyperbolic rotation using orthogonal-diagonal
                                // method (forward numerical stable) (from [2])
                                double a = 0.5 * sgn2( a1, tol );
                                t = a * std::sqrt( ( a1 + a2 ) / ( a1 - a2 ) );
                                z = a * std::sqrt( ( a1 - a2 ) / ( a1 + a2 ) );
                                x = R.row( r1-1 ) - R.row( r2-1 );
                                y = R.row( r1-1 ) + R.row( r2-1 );
                                x = t * x;
                                y = z * y;
                                R.row( r1-1 ) = ( x + y );
                                R.row( r2-1 ) = ( y - x );
                                R( r2-1, c-1 ) = 0.0;
                                x = Q.row( r1-1 ) - Q.row( r2-1 );
                                y = Q.row( r1-1 ) + Q.row( r2-1 );
                                x = t * x;
                                y = z * y;
                                Q.row( r1-1 ) = ( x + y );
                                Q.row( r2-1 ) = ( y - x );
                            } else {
                                // Apply hyperbolic rotation using mixed-downdating
                                t = a2 / a1;
                                z = sgn2( a1, tol ) * std::sqrt( 1.0 - ( t * t ) );
                                R.row( r1-1 ) = ( R.row( r1-1 ) - ( t * R.row( r2-1 ) ) ) / z;
                                R.row( r2-1 ) = ( ( -t ) * R.row( r1-1 ) ) + ( z * R.row( r2-1 ) );
                                R( r2-1, c-1 ) = 0.0;
                                Q.row( r1-1 ) = ( Q.row( r1-1 ) - ( t * Q.row( r2-1 ) ) ) / z;
                                Q.row( r2-1 ) = ( ( -t ) * Q.row( r1-1 ) ) + ( z * Q.row( r2-1 ) );
                            }
                        }
                    } else { // No hyperbolic rotation
                        if( ( r2 < m ) || ( retry > 0 ) ) { // Still entries remaining
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
    } // while

    // Make last diagonal entry positive
    if( R( r-1, c-1 ) < 0.0 ) {
        Q.row( r-1 ) = -Q.row( r-1 );
        R.row( r-1 ) = -R.row( r-1 );
    }

    // Post-processing
    J1 = arma::diagmat( J );
    arma::mat tmp = arma::diagmat( J ) * Q.t() * arma::diagmat( J );
    for( int i = 0; i < m; i++ ) {
        Q.row( perm( i )-1 ) = tmp.row( i );
    }
    if( econ && ( m > n ) ) {
        R = R( arma::span( 0, ( n - 1 ) ), arma::span( 0, ( n - 1 ) ) );
        Q = Q( arma::span(), arma::span( 0, ( n - 1 ) ));
        J1 = J1( arma::span( 0, ( n - 1 ) ), arma::span( 0, ( n - 1 ) ) );
    }
    return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int J_orthogonal( arma::mat &Q, arma::mat &R, arma::mat &J1, const arma::mat &A, const arma::mat &J0, bool econ, double tol )
{
    int m = A.n_rows;
    int n = A.n_cols;

    arma::vec J = J0.diag();

    bool od_method = true; // Use orthogonal-diagonal method for real matrices

    // Initialize matrices
    R = arma::mat( A ); // upper-triangular matrix
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

                if( retry > 0 ) {
                    std::cout << "retry > 2 JQR!" << std::endl;
                    return 1;
                }
//                retry--;
                std::cout << "Retry JQR!" << std::endl;
//                assert( false ); // ???
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
                if( J( r1 ) == J( r2 ) ) { // Same sign2, apply orthogonal rotation

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
                            J.swap_rows( r2, r1 );
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
    J1 = arma::diagmat( J );
    arma::mat tmp = arma::diagmat( J ) * Q.t() * arma::diagmat( J );
    for( int i = 0; i < m; i++ ) {
        Q.row( perm( i ) ) = tmp.row( i );
    }
    if( econ && ( m > n ) ) {
        R = R( arma::span( 0, ( n - 1 ) ), arma::span( 0, ( n - 1 ) ) );
        Q = Q( arma::span(), arma::span( 0, ( n - 1 ) ));
        J1 = J1( arma::span( 0, ( n - 1 ) ), arma::span( 0, ( n - 1 ) ) );
    }
    return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int ModifiedGramSchmidt( arma::mat &Q, arma::mat &R, const arma::mat &X )
{
    int n = X.n_rows;
    int p = X.n_cols;
    Q = arma::mat( n, p, arma::fill::zeros );
    R = arma::mat( p, p, arma::fill::zeros );
    for( int k = 0; k < p; k++ ) {
        Q.col(k) = X.col(k);
        for( int i = 0; i < k; i++ ) {
            R(i, k) = arma::dot( arma::trans( Q.col(i) ), Q.col(k) );
            Q.col(k) = Q.col(k) - R(i, k) * Q.col(i);
        }
        R(k, k) = ( arma::norm( Q.col(k) ) );
        Q.col(k) = Q.col(k) / R(k, k);
    }
    return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int SchwarzRutishauser( arma::mat &Q, arma::mat &R, const arma::mat &A, double tol )
{
    int rowsA = A.n_rows;
    int colsA = A.n_cols;
    Q = arma::mat( rowsA, colsA, arma::fill::zeros );
    R = arma::mat( colsA, colsA, arma::fill::zeros );
    int orig, flag;
    double s = 0.0;
    double t = 0.0;
    int max_iter = 20;
    int iter = 0;
    for( int k = 0; k < colsA; k++ ) {
        for( int j = 0; j < rowsA; j++ ) {
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
                for( int j = 0; j < rowsA; j++ ) {
                    s += ( Q(j, i) * Q(j, k) );
                }
                if( orig ) {
                    R(i, k) = s;
                } else {
                    R(i, k) = R(i, k) + s;
                }
                t += ( s * s );
                for( int j = 0; j < rowsA; j++ ) {
                    Q(j, k) -= ( s * Q(j, i) );
                }
            } // end for i
            s = 0.0;
            for( int j = 0; j < rowsA; j++ ) {
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
        for( int j = 0; j < rowsA; j++ ) {
            Q(j, k) *= s;
        }
//        std::cout << iter << std::endl;
    }
    return 0;
}

}
}
/// \}
