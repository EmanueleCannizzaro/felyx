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
//
//

#ifndef ITL_LDL_CHOLESKY_H
#define ITL_LDL_CHOLESKY_H

#include <vector>
#include <assert.h>
#include <algorithm>

#include "itl/itl.h"
#include "itl/preconditioner/detail/preconditioner.h"
#include "mtl/meta_if.h"

//#include <boost/timer.hpp>


namespace itl {

  /* VC++ work around for ambiguous function */
  /*
   inline void check_symm(mtl::symmetric_tag) { }

    template < class Shape >
    inline void check_symm(Shape)
    {
      std::cout << "Matrix is not symmetric. Abort." << std::endl;
      assert(0);
    }
   */

  //: Incomplete Cholesky Preconditioner, LDL' variant with zero fill.
  //  For use with symmetric matrices.
  //  This variant is much more stable than cholesky, since it doesn't use any square roots.
  //  
  //
  //<codeblock>
  // Usage:
  //    SymMatrix A;
  //    ldl_cholesky< SymMartix > precond(A);
  //    cg(A, x, b, precond(), iter);
  //</codeblock>
  //
  //Notes: The idea under a concrete Preconditioner such 
  //as Incomplete Cholesky is to create a Preconditioner
  //object to use in iterative methods. 
  //
  //
  //!definition: ldl_cholesky.h
  //!example: cholesky.cc
  //!category: itl,functors
  //!component: type
  //!tparam: Matrix - A symmetric Matrix 
  //
  template < class Matrix >
  class ldl_cholesky {
    typedef typename Matrix::value_type T;
    typedef typename Matrix::orientation Orien;
    typedef Matrix SymMatrix;
    enum { Orien_id = Orien::id };

    typedef typename mtl::matrix< T, mtl::rectangle<>, 
                         mtl::compressed<int, mtl::external>, mtl::column_major >::type Matrix1; 

    typedef typename mtl::matrix< T, mtl::rectangle<>, 
                         mtl::compressed<int, mtl::external>, mtl::row_major >::type Matrix2; 

    typedef typename mtl::IF< EQUAL < Orien_id, mtl::ROW_MAJOR >::RET, 
                         Matrix2,
            typename mtl::IF< EQUAL < Orien_id, mtl::COL_MAJOR >::RET, 
                         Matrix1,
			 mtl::generators_error
    >::RET
    >::RET TriMatrix;


  public:

    typedef preconditioner3 < Matrix1, Matrix2, mtl::lower, mtl::upper> Precond;
    typedef preconditioner1< Matrix1, Matrix2, mtl::lower, mtl::upper> Left;
    typedef preconditioner2< Matrix1, Matrix2, mtl::lower, mtl::upper> Right;
  
    ldl_cholesky(const SymMatrix& A) 
      : Tri_val(A.nnz()-A.nrows()), Tri_ind(A.nnz()-A.nrows()), Tri_ptr(A.nrows()+1), Adiag(A.nrows()), diag(A.nrows(),0), Y(A.nrows(),0)
    {
      //typedef typename mtl::matrix_traits<SymMatrix>::shape Shape;
      //check_symm(Shape()); // JGS change to compile-time test
      do_cholesky(A, Orien()); 
    }
  private:
 
    void do_cholesky(const SymMatrix& A, mtl::row_tag) {
      using std::copy;
      if ( A.is_upper() ) {
	int Tri_loc= 0;
	Tri_ptr[0] = 0;
	
	typename SymMatrix::const_iterator A_i = A.begin();
	
        while (mtl::not_at(A_i,A.end())) {
          typename SymMatrix::OneD A_row = *A_i;
          
          int i=A_i.index(); 
          Adiag[i] = *(A_row.begin());
          std::copy(A_row.begin()+1, A_row.end(), Tri_val.begin()+Tri_loc);
          std::copy(A_row.nz_struct().begin()+1, A_row.nz_struct().end(), 
		    Tri_ind.begin()+Tri_loc);
          
          Tri_loc += A_row.nnz()-1;
          Tri_ptr[i+1] = Tri_loc;
          
          ++A_i;
        }
        
	int d, g, h, i, j, k, n = A.nrows();
        
        for (k = 0; k < n; k++) {
          d = Tri_ptr[k];
          diag[k] = 1/(Adiag[k] - diag[k]);
	  
          for (i = d; i < Tri_ptr[k+1]; i++) {
            Y[Tri_ind[i]] = Tri_val[i];
	    Tri_val[i] = Tri_val[i]*diag[k];
          }
          
          for (i = d; i < Tri_ptr[k+1]; i++) {
            h = Tri_ind[i];
            g = i;
            diag[Tri_ind[i]] += Tri_val[i]*Y[Tri_ind[i]];
            
            for (j = Tri_ptr[h]; j < Tri_ptr[h+1]; j++)
              for ( ; g < Tri_ptr[k+1] && Tri_ind[g] <= Tri_ind[j]; g++)
                if (Tri_ind[g] == Tri_ind[j])
                  Tri_val[j] -= Y[Tri_ind[i]] * Tri_val[g];
          }
        }
        
        std::cout << "--< " << my_timer.elapsed() << " >-- " << ": Preconditioning done. ";
	std::cout << "Non-Zeros: " << Tri_loc << endl;

        Tri = TriMatrix(A.nrows(), A.ncols(), 
                        Tri_loc, &Tri_val[0], &Tri_ptr[0], 
                        &Tri_ind[0]);        
      } else { 
        std::cout << "Warning: It is not so efficient as symmetric row-wise upper Matrix" << std::endl;
        assert(0);
      }    
    }

    void do_cholesky(const SymMatrix& A, mtl::column_tag) {
      using std::copy;
      if ( A.is_upper() ) {
	std::cout << "Warning: It is not so efficient as symmetric column-wise lower Matrix" << std::endl;
	assert(0);
      } else {
	int Tri_loc= 0;
	Tri_ptr[0] = 0;
	
	typename SymMatrix::const_iterator A_i = A.begin();
	
        while (mtl::not_at(A_i,A.end())) {
          typename SymMatrix::OneD A_row = *A_i;
          
          int i=A_i.index(); 
          Adiag[i] = *(A_row.begin());
          std::copy(A_row.begin()+1, A_row.end(), Tri_val.begin()+Tri_loc);
          std::copy(A_row.nz_struct().begin()+1, A_row.nz_struct().end(), 
		    Tri_ind.begin()+Tri_loc);
          
          Tri_loc += A_row.nnz()-1;
          Tri_ptr[i+1] = Tri_loc;
          
          ++A_i;
        }
        
	int d, g, h, i, j, k, n = A.nrows();
        
        for (k = 0; k < n; k++) {
          d = Tri_ptr[k];
          diag[k] = 1/(Adiag[k] - diag[k]);
	  
          for (i = d; i < Tri_ptr[k+1]; i++) {
            Y[Tri_ind[i]] = Tri_val[i];
	    Tri_val[i] = Tri_val[i]*diag[k];
          }
          
          for (i = d; i < Tri_ptr[k+1]; i++) {
            h = Tri_ind[i];
            g = i;
            diag[Tri_ind[i]] += Tri_val[i]*Y[Tri_ind[i]];
            
            for (j = Tri_ptr[h]; j < Tri_ptr[h+1]; j++)
              for ( ; g < Tri_ptr[k+1] && Tri_ind[g] <= Tri_ind[j]; g++)
                if (Tri_ind[g] == Tri_ind[j])
                  Tri_val[j] -= Y[Tri_ind[i]] * Tri_val[g];
          }
        }
        
        std::cout << "--< " << my_timer.elapsed() << " >-- " << ": Preconditioning done. ";
	std::cout << "Non-Zeros: " << Tri_loc << endl;

        Tri = TriMatrix(A.nrows(), A.ncols(), 
                        Tri_loc, &Tri_val[0], &Tri_ptr[0], 
                        &Tri_ind[0]);
      }
    }

    inline Precond pre3(mtl::row_tag) {
      return  Precond(mtl::trans(Tri), Tri, diag);
    }
    
    inline Precond pre3(mtl::column_tag) {
      return  Precond(Tri, mtl::trans(Tri), diag);
    }

    inline Left pre1(mtl::row_tag) {
      return  Left(mtl::trans(Tri), Tri);
    }
    inline Left pre1(mtl::column_tag) {
      return  Left(Tri, mtl::trans(Tri));
    }
    inline Right pre2(mtl::row_tag) {
      return  Right(mtl::trans(Tri), Tri);
    }
    inline Right pre2(mtl::column_tag) {
      return  Right(Tri, mtl::trans(Tri));
    }


  public:
    //:return a left or right Preconditioner object
    inline Precond operator()() { 
      return pre3(Orien());
    }
    //: return a left part of Split Preconditioner object
    inline Left left() { return pre1(Orien()); }
    //: return a right part of Split Preconditioner object
    inline Right right() { return pre2(Orien()); }

    void print() {
    }

  private:
    TriMatrix Tri;
    std::vector<T> Tri_val;
    std::vector<int> Tri_ind;
    std::vector<int> Tri_ptr;
    std::vector<T> Adiag;
    mtl::dense1D<T> diag;  //storage for the diagonal elements
    std::vector<T> Y;      //storage for multipliers used during calculation
    boost::timer my_timer;
  };


}

#endif
