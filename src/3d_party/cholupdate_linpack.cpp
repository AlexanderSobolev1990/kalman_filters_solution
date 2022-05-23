//----------------------------------------------------------------------------------------------------------------------
///
/// \file       cholupdate_linpack.cpp
/// \brief      Процедура cholupdate (осуществляет вызов linpack процедур: dchud для update и dchdd для downdate)
/// \date       17.02.22 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#include <cholupdate_linpack.h>

int cholupdate_linpack( arma::mat &A, const arma::vec &X, double v )
{
    assert( A.n_cols == A.n_rows );
    assert( A.n_cols == X.n_elem );

    int p = A.n_cols;
    int ldr = p;
    int nz = 1;
    int ldz = p;
    arma::mat z = arma::mat( ldz, nz, arma::fill::zeros );
    arma::vec y = arma::vec( nz, arma::fill::zeros );
    arma::vec rho = arma::vec( nz, arma::fill::zeros );
    arma::vec c = arma::vec( p, arma::fill::zeros );
    arma::vec s = arma::vec( p, arma::fill::zeros );

    arma::vec vx = X * std::sqrt( std::abs( v ) );

    A = arma::trans( A );

    int res = 0;
    if( v > 0.0 ) {  // Update
        dchud( A.memptr(), ldr, p, vx.memptr(), z.memptr(), ldz, nz, y.memptr(), rho.memptr(), c.memptr(), s.memptr() );
    } else {  // Downdate
        res = dchdd( A.memptr(), ldr, p, vx.memptr(), z.memptr(), ldz, nz, y.memptr(), rho.memptr(), c.memptr(), s.memptr() );
    }
    ///////////////
    // debug печать
//    std::cout << "A:" << std::endl;
//    A.print();

//    std::cout << "Switch A to lower:" << std::endl;
    A = arma::trans( A );
//    A.print();
    ///////////////

    return std::abs( res );
}
/// \}
