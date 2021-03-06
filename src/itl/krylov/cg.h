//                                                                           
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu> 
//          Lie-Quan Lee     <llee@osl.iu.edu>
//
// This file is part of the Iterative Template Library
//
// You should have received a copy of the License Agreement for the
// Iterative Template Library along with the software;  see the
// file LICENSE.  
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
// Modified by U.Mennel and O.Koenig to use BLAS functions Dec 2004
//=========================================================================

#ifndef ITL_CG_H
#define ITL_CG_H

#include "itl/itl.h"

#include "BlasHeaders.h"

namespace itl {

  //: Conjugate Gradient(CG)
  //
  //Solves the symmetric positive definite linear system A x = b.
  //
  //<p>A return value of 0 indicates convergence within the
  //maximum number of iterations (determined by the iter object).
  //A return value of 1 indicates a failure to converge.
  //
  //<p>See M. R. Hestenes nd E. Stiefel, Methods of conjugate gradients for 
  //solving linear system, Journal of Research of the National Bureau of 
  //Standards, 49(1952), pp. 409-436
  //
  //
  //!category: itl,algorithms
  //!component: function
  //!definition: cg.h
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods
  //!tparam: VectorX - Vector 
  //!tparam: VectorB - Vector
  //!tparam: Preconditioner - Current choices are Incomplete Cholesky or identity_preconditioner.
  //!tparam: Iteration - Controls the stopping criteria
  //
  /* required operations: mult,copy,dot_conj,add,scaled */

  
  template < class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration >
    int cg(const Matrix& A, VectorX& x, const VectorB& b, const Preconditioner& M, Iteration& iter)
  {
    
    typedef VectorX TmpVec;
    typename itl_traits<VectorX>::value_type rho(0), rho_1(0), alpha(0), beta(0);
    TmpVec p(size(x)), q(size(x)), r(size(x)), z(size(x));

#ifdef HAVE_BLAS 
    int n = size(x);
#endif
   
    itl::mult(A, itl::scaled(x, -1.0), b, r);     
    
    while (! iter.finished(r)) {
    
      itl::solve(M, r, z);

#ifdef HAVE_BLAS
#ifdef USE_FLOATS
      rho = cblas_sdot(n, &(*r.begin()), 1, &(*z.begin()), 1);
#else
      rho = cblas_ddot(n, &(*r.begin()), 1, &(*z.begin()), 1);
#endif
#else
      rho = itl::dot_conj(r, z); //dense scalar
#endif
      
    if (iter.first())
      itl::copy(z, p);            
    else {
      beta = rho / rho_1;
      itl::add(z, itl::scaled(p, beta), p);
    }
    
    itl::mult(A, p, q); //sparse matrix vector product    

#ifdef HAVE_BLAS
#ifdef USE_FLOATS
    alpha = rho / cblas_sdot(n, &(*p.begin()), 1, &(*q.begin()), 1);
    cblas_saxpy( n, alpha, &(*p.begin()), 1, &(*x.begin()), 1 );
    cblas_saxpy( n, -alpha, &(*q.begin()), 1, &(*r.begin()), 1 );
#else
    alpha = rho / cblas_ddot(n, &(*p.begin()), 1, &(*q.begin()), 1);
    cblas_daxpy( n, alpha, &(*p.begin()), 1, &(*x.begin()), 1 );
    cblas_daxpy( n, -alpha, &(*q.begin()), 1, &(*r.begin()), 1 );
#endif
#else
    alpha = rho / itl::dot_conj(p, q); //dense scalar product
    itl::add(x, itl::scaled(p, alpha), x);  //dense daxpy x=p*alpha + x
    itl::add(r, itl::scaled(q, -alpha), r); //dense daxpy r=-q*alpha + r 
#endif
        
    rho_1 = rho;
    
    ++iter;
    }
    
    return iter.error_code();
  }
  
} /* namespace itl */

#endif

