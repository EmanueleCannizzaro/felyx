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

#ifndef ITL_SAINV_H
#define ITL_SAINV_H

#include <vector>
#include <assert.h>
#include <algorithm>

#include "itl/itl.h"
#include "itl/preconditioner/detail/preconditioner.h"
#include "mtl/meta_if.h"

//#include <boost/timer.hpp>

#define MAX(A,B) ( (A) > (B) ? (A):(B))
#define MIN(A,B) ( (A) < (B) ? (A):(B))
#define ABS(a) ((a)>=0 ? (a) : -(a))
 
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


  //: Stabilized ,outer-product AINV.
  //  For use with symmetric matrices.
  //  This factorization never brakes down, and is well defined for any
  //  kind of symmtric positive definite Matrix.
  //  
  //
  //<codeblock>
  // Usage:
  //    SymMatrix A;
  //    sainv< SymMartix > precond(A);
  //    cg(A, x, b, precond(), iter);
  //</codeblock>
  //
  //Notes: The idea under a concrete Preconditioner such 
  //as Incomplete Cholesky is to create a Preconditioner
  //object to use in iterative methods. 
  //
  //
  //!definition: sainv.h
  //!example: cholesky.cc
  //!category: itl,functors
  //!component: type
  //!tparam: Matrix - A symmetric Matrix 
  //
  template < class Matrix >
  class sainv {
    typedef typename Matrix::value_type T;
    typedef typename Matrix::orientation Orien;
    typedef Matrix SymMatrix;
    enum { Orien_id = Orien::id };

    typedef typename mtl::matrix< T, mtl::rectangle<>, 
                         mtl::array< compressed<> >, mtl::column_major >::type Matrix2; 

    typedef typename mtl::matrix< T, mtl::rectangle<>, 
                         mtl::array< compressed<> >, mtl::row_major >::type Matrix1; 

    typedef typename mtl::IF< EQUAL < Orien_id, mtl::ROW_MAJOR >::RET, 
                         Matrix1,
            typename mtl::IF< EQUAL < Orien_id, mtl::COL_MAJOR >::RET, 
                         Matrix2,
			 mtl::generators_error
    >::RET
    >::RET TriMatrix;


  public:

    typedef preconditioner4 < Matrix1, Matrix2, mtl::lower, mtl::upper> Precond;
    typedef preconditioner1< Matrix1, Matrix2, mtl::lower, mtl::upper> Left;
    typedef preconditioner2< Matrix1, Matrix2, mtl::lower, mtl::upper> Right;
  
    sainv(const SymMatrix& A, T tol) 
      : Tri(A.nrows(), A.ncols()), diag(A.nrows(),0), norm(A.nrows(),0)
    {
      //typedef typename mtl::matrix_traits<SymMatrix>::shape Shape;
      //check_symm(Shape()); // JGS change to compile-time test
      do_sainv(A, tol, Orien()); 
    }
  private:
 
    void 
    do_sainv(const SymMatrix& A, double tol, mtl::row_tag)
    {
      if ( A.is_upper() ) {
      
      int j,i;
	int n = A.ncols();
	T z;
	
	typename SymMatrix::const_iterator A_i = A.begin();
	
	while (mtl::not_at(A_i,A.end())) {
          typename SymMatrix::OneD A_row = *A_i;
          
          int i=A_i.index(); 
          A_row.filter(0); //must invoke this to prevent too many drops
          
          //calculate the norm vector
          norm[i] = mtl::one_norm(A_row);
	  norm[i] *= tol;
	  norm[i] /= A_row.nnz();
          
          ++A_i;
        }

	typename SymMatrix::OneD::iterator A_it;
	typename TriMatrix::OneD::iterator Tri_it, Tri_it_j;
	typename compressed1D<T>::iterator Tri_col_it, temp_it;
	
	//Load identity
	typename TriMatrix::iterator Tri_i = Tri.begin();
	while (Tri_i != Tri.end()) {
	  typename TriMatrix::OneD Tri_col = *Tri_i;
	  Tri_col.push_back(Tri_i.index(),1);
	  ++Tri_i;
	}

	for (i=0;i < n; ++i) {
	
	  //extract col ii from Tri
	  typename TriMatrix::iterator Tri_i = Tri.begin()+i;
	  typename TriMatrix::OneD Tri_col = *Tri_i;
          
          /*
	  compressed1D<T> Tri_col_i;
          Tri_col_it=Tri_col.begin();
          Tri_col_i.push_back(Tri_col_it.index(), *Tri_col_it);
          ++Tri_col_it;
	  for (; Tri_col_it != Tri_col.end(); ++Tri_col_it) {
	    //if (ABS(*Tri_col_it) > norm[i]) {
	      Tri_col_i.push_back(Tri_col_it.index(), *Tri_col_it);
	    //}
	  }
	  */
          
	  compressed1D<T> Tri_col_i = compressed1D<T>(Tri_col.nz_struct().begin(), Tri_col.nz_struct().end(), Tri_col.nnz());
	  copy(Tri_col.begin(), Tri_col.end(), Tri_col_i.begin());
	  
	  compressed1D<T>u;
	
	  //typename TriMatrix::OneD::const_iterator u_it = u.begin();
	  //calculate u=Az_i

	  //for (Tri_col_it=Tri_col_i.begin(); Tri_col_it != Tri_col_i.end(); ++Tri_col_it) {
	  //ii=Tri_col_it.index();
	  mtl::mult(A, Tri_col_i, u);
	  
          
	  z = mtl::dot(u, Tri_col_i);
          if (z==0) std::cout << "WARNING: Pivot brakedown!" << endl;
          //st
	  diag[i]=1/z;

	  
	  mtl::scale(u, diag[i]);
	  
	  dense1D<T> c(n,0);
	  //calculating coeffs
	  for (j=i+1; j<n; ++j) {

	    typename TriMatrix::iterator Tri_j = Tri.begin()+j;
	    typename TriMatrix::OneD Tri_col = *Tri_j;
	    //z=0;
	    
            
	    c[j] = -mtl::dot(u, Tri_col);
	  }
	  
	  for (j=i+1;j<n;++j) {
	  
	    typename TriMatrix::iterator Tri_j = Tri.begin()+j;
            typename TriMatrix::OneD Tri_col = *Tri_j;
	    
	    compressed1D<T> temp;
            if (c[j]) {
	      //z = ABS(tol/c[j]);
	      for (Tri_col_it=Tri_col_i.begin(); Tri_col_it != Tri_col_i.end(); ++Tri_col_it) {
	        if (ABS(*Tri_col_it) > z) {
		  temp.push_back(Tri_col_it.index(), *Tri_col_it*c[j]);
	        }
	      }
            }
            
            //Postfiltering to reduce memory needs
            if (j == i+1) {
              mtl::add(Tri_col, temp, temp);
              Tri_col.clear();
              temp_it=temp.begin();
              //don't filter diagonal element
              Tri_col.push_back(temp_it.index(), *temp_it);
              ++temp_it;
              for (; temp_it != temp.end(); ++temp_it) {
	        if (ABS(*temp_it) > tol) {
  		  Tri_col.push_back(temp_it.index(), *temp_it);
	        }
	      }
            } else {
              mtl::add(Tri_col, temp, Tri_col);
            }
	  }
        }
        std::cout << "--< " << my_timer.elapsed() << " >-- " << ": Preconditioning done. ";
	std::cout << "Non-Zeros: " << Tri.nnz() << endl;
	//std::cout << "Items dropped: " << drop << endl;
	//std::cout << "Total items calculated: " << calcs << endl;
      } 
      else { 
	std::cout << "Warning: It is not so efficient as symmetric row-wise upper Matrix" << std::endl;
	assert(0);
      }
    }

    void 
    do_sainv(const SymMatrix& A, double tol, mtl::column_tag)
    {
      using std::copy;
      if ( A.is_upper() ) {
        std::cout << "Warning: It is not so efficient as symmetric column-wise lower Matrix" << std::endl;
        assert(0);
      } 
      else {
        
	int j,i;
	int n = A.ncols();
	T z;
	
	typename SymMatrix::const_iterator A_i = A.begin();
	
	while (mtl::not_at(A_i,A.end())) {
          typename SymMatrix::OneD A_row = *A_i;
          
          int i=A_i.index(); 
          A_row.filter(0); //must invoke this to prevent too many drops
          
          //calculate the norm vector
          norm[i] = mtl::one_norm(A_row);
	  norm[i] *= tol;
	  norm[i] /= A_row.nnz();
          
          ++A_i;
        }

	typename SymMatrix::OneD::iterator A_it;
	typename TriMatrix::OneD::iterator Tri_it, Tri_it_j;
	typename compressed1D<T>::iterator Tri_col_it, temp_it;
	
	//Load identity
	typename TriMatrix::iterator Tri_i = Tri.begin();
	while (Tri_i != Tri.end()) {
	  typename TriMatrix::OneD Tri_col = *Tri_i;
	  Tri_col.push_back(Tri_i.index(),1);
	  ++Tri_i;
	}

	for (i=0;i < n; ++i) {
	
	  //extract col ii from Tri
	  typename TriMatrix::iterator Tri_i = Tri.begin()+i;
	  typename TriMatrix::OneD Tri_col = *Tri_i;
          
          /*
	  compressed1D<T> Tri_col_i;
          Tri_col_it=Tri_col.begin();
          Tri_col_i.push_back(Tri_col_it.index(), *Tri_col_it);
          ++Tri_col_it;
	  for (; Tri_col_it != Tri_col.end(); ++Tri_col_it) {
	    //if (ABS(*Tri_col_it) > norm[i]) {
	      Tri_col_i.push_back(Tri_col_it.index(), *Tri_col_it);
	    //}
	  }
	  */
          
	  compressed1D<T> Tri_col_i = compressed1D<T>(Tri_col.nz_struct().begin(), Tri_col.nz_struct().end(), Tri_col.nnz());
	  copy(Tri_col.begin(), Tri_col.end(), Tri_col_i.begin());
	  
	  compressed1D<T>u;
	
	  //typename TriMatrix::OneD::const_iterator u_it = u.begin();
	  //calculate u=Az_i

	  //for (Tri_col_it=Tri_col_i.begin(); Tri_col_it != Tri_col_i.end(); ++Tri_col_it) {
	  //ii=Tri_col_it.index();
	  mtl::mult(A, Tri_col_i, u);
	  
          
	  z = mtl::dot(u, Tri_col_i);
          if (z==0) std::cout << "WARNING: Pivot brakedown!" << endl;

	  diag[i]=1/z;

	  
	  mtl::scale(u, diag[i]);
	  
	  dense1D<T> c(n,0);
	  //calculating coeffs
	  for (j=i+1; j<n; ++j) {

	    typename TriMatrix::iterator Tri_j = Tri.begin()+j;
	    typename TriMatrix::OneD Tri_col = *Tri_j;
	    //z=0;
	    
            
	    c[j] = -mtl::dot(u, Tri_col);
	  }
	  
	  for (j=i+1;j<n;++j) {
	  
	    typename TriMatrix::iterator Tri_j = Tri.begin()+j;
            typename TriMatrix::OneD Tri_col = *Tri_j;
	    
	    compressed1D<T> temp;
            if (c[j]) {
	      //z = ABS(tol/c[j]);
	      for (Tri_col_it=Tri_col_i.begin(); Tri_col_it != Tri_col_i.end(); ++Tri_col_it) {
	        if (ABS(*Tri_col_it) > z) {
		  temp.push_back(Tri_col_it.index(), *Tri_col_it*c[j]);
	        }
	      }
            }
            
            //Postfiltering to reduce memory needs
            if (j == i+1) {
              mtl::add(Tri_col, temp, temp);
              Tri_col.clear();
              temp_it=temp.begin();
              //don't filter diagonal element
              Tri_col.push_back(temp_it.index(), *temp_it);
              ++temp_it;
              for (; temp_it != temp.end(); ++temp_it) {
	        if (ABS(*temp_it) > tol) {
  		  Tri_col.push_back(temp_it.index(), *temp_it);
	        }
	      }
            } else {
              mtl::add(Tri_col, temp, Tri_col);
            }
	  }
	}
	std::cout << "--< " << my_timer.elapsed() << " >-- " << ": Preconditioning done. ";
	std::cout << "Non-Zeros: " << Tri.nnz() << endl;
	
      }
    }
    

    inline Precond pre4(mtl::column_tag) {
      return  Precond(mtl::trans(Tri), Tri, diag);
    }
    
    inline Precond pre4(mtl::row_tag) {
      return  Precond(Tri, mtl::trans(Tri), diag);
    }

    /*
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
    */

  public:
    //:return a left or right Preconditioner object
    inline Precond operator()() { 
      return pre4(Orien());
    }
    //: return a left part of Split Preconditioner object
    //inline Left left() { return pre1(Orien()); }
    //: return a right part of Split Preconditioner object
    //inline Right right() { return pre2(Orien()); }

    void print() {
    }

  private:
    TriMatrix Tri;
    mtl::dense1D<T> diag;  //storage for the diagonal elements
    std::vector<T> norm;
    boost::timer my_timer;
  };


}

#endif
