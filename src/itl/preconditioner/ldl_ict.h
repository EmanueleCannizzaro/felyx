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

#ifndef ITL_LDL_ICT_H
#define ITL_LDL_ICT_H

#include <vector>
#include <assert.h>
#include <algorithm>

#include "itl/itl.h"
#include "itl/preconditioner/detail/preconditioner.h"
#include "mtl/meta_if.h"


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

  //: Incomplete Cholesky Preconditioner, LDL' variant with dual threshold.
  //  For use with symmetric matrices.
  //  This preconditioner is much more accurate than the zero fill variant.
  //  
  //
  //<codeblock>
  // Usage:
  //    SymMatrix A;
  //    ldl_cholesky< SymMartix > precond(A, lfil, tol);
  //    cg(A, x, b, precond(), iter);
  //</codeblock>
  //
  //Notes: The idea under a concrete Preconditioner such 
  //as Incomplete Cholesky is to create a Preconditioner
  //object to use in iterative methods. 
  //
  //
  //!definition: ldl_cholesky_t.h
  //!example: cholesky.cc
  //!category: itl,functors
  //!component: type
  //!tparam: Matrix - A symmetric Matrix 
  //!tparam: lfil - Level of fill
  //!tparam: force - continue factorization, if negative pivot occured
  //!tparam: tol - Drop tolerance
  template < class Matrix >
  class ldl_ict {
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
  
    ldl_ict(const SymMatrix& A, int lfil=10, T tol=1e-6) 
      : Tri_ptr(A.nrows()+1), 
      diag(A.nrows()),
      w(A.nrows()),
      jw(A.nrows()),
      jr(A.nrows())
    {
      int alloc = MIN(A.nnz()+A.nrows()*lfil, Tri_val.max_size());
      
      if (alloc < 0) {  //bad lfil value??
        alloc = Tri_val.max_size();
      }
      Tri_val.resize(alloc);
      Tri_ind.resize(alloc);
      //typedef typename mtl::matrix_traits<SymMatrix>::shape Shape;
      //check_symm(Shape()); // JGS change to compile-time test
      do_cholesky(A, Orien(), lfil, tol); 
    }
  private:
 
    void 
    do_cholesky(const SymMatrix& A, mtl::row_tag, int lfil, T tol)
    {
      using std::copy;
      if ( A.is_upper() ) {
      int Tri_loc= 0;
        Tri_ptr[0] = 0;
        
        typename SymMatrix::const_iterator A_i = A.begin();
                
        while (mtl::not_at(A_i,A.end())) {
          typename SymMatrix::OneD A_row = *A_i;

          A_row.filter(0); //must invoke this to prevent too many drops
          
          ++A_i;
        }
        
                
        int d, l, jj, n = A.nrows();

        std::vector<int> jj_index(n,0);
        std::vector<int> jj_pos(n,0);
        typename std::vector<int>::iterator jj_index_it;
        typename std::vector<int>::iterator jj_index_end;
        typename std::vector<int>::iterator jj_pos_it;
        typename mtl::dense1D<T>::iterator diag_it = diag.begin();
        typename SymMatrix::const_iterator j = A.begin();
        typename SymMatrix::OneD::iterator i;
        T z, t, c,c_norm;
        
        for (j = A.begin(); mtl::not_at(j,A.end()); j++) {
          typename SymMatrix::OneD j_row = *j;
          jj = j.index();
          
          wnum=0;
          i = j_row.begin();
          *diag_it = *i - *diag_it;
          
          i++;
          //fill work structure with non-zeros
          std::copy(i, j_row.end(), w.begin());
          std::copy(j_row.nz_struct().begin()+1, j_row.nz_struct().end(), jw.begin());
          
          for (; i != j_row.end(); ++i) {
            ++wnum;
            jr[i.index()] = wnum;
          }
          
        
          int nonzeros = wnum;
        
          c_norm = mtl::one_norm(*j)/ (*j).nnz()*tol;
          typename std::vector<int>::iterator tri_ptr_it = Tri_ptr.begin();
          typename mtl::dense1D<T>::iterator diag_it2 = diag.begin();
          typename std::vector<int>::iterator jr_it;
        
          
          jj_index_it = jj_index.begin();
          jj_index_end = jj_index_it+jj;
          jj_pos_it = jj_pos.begin();
          for (jj_index_it = jj_index.begin(); jj_index_it != jj_index_end; ++jj_index_it)
           {
            ++tri_ptr_it;
            if (*jj_index_it == jj) { // must be this element or not in row
              z = Tri_val[*jj_pos_it]* *diag_it2;
              c = ABS(c_norm/z);
              ++*jj_pos_it;
              *jj_index_it = Tri_ind[*jj_pos_it];
           
              for (l = *jj_pos_it; l < *tri_ptr_it; ++l) {
                jr_it = jr.begin()+Tri_ind[l];
                if (*jr_it) {
                  //already in work vector
                  w[(*jr_it)-1] -= Tri_val[l]*z;
                }
                else {
                  //insert this element, if large enough
                  if (ABS(Tri_val[l]) > c) {
                    w[wnum] = -Tri_val[l]*z;
                    jw[wnum] = Tri_ind[l];
                    ++wnum;
                    *jr_it = wnum;
                  }
                }
              }
            }
            ++diag_it2;
            ++jj_pos_it;
          }
          
                  
          d = MIN(j_row.nnz()+lfil-1, wnum);
          if (d < wnum  && d > nonzeros) {
            qsplit(nonzeros, wnum, d);
          }
          
          
          typename std::vector<int>::iterator jw_it = jw.begin();
          typename std::vector<int>::iterator jw_end = jw_it+d;
          std::sort(jw_it, jw_end);
          
          Tri_loc = Tri_ptr[jj];
          typename std::vector<T>::iterator tri_val_it = Tri_val.begin()+Tri_loc;
          typename std::vector<int>::iterator tri_ind_it = Tri_ind.begin()+Tri_loc;
          
          
          
          *jj_pos_it = Tri_loc;
          *jj_index_it = *jw.begin();
          
          std::copy(jw_it, jw_end, tri_ind_it); 
          for (;jw_it != jw_end; ++jw_it) {
          
            t = w[jr[*jw_it] -1];
              z = t/ *diag_it;
              *tri_val_it = z;
              ++tri_val_it;
              ++Tri_loc;
              diag[*jw_it] += t *z;
          }
          Tri_ptr[jj+1] = Tri_loc;
          
          
          std::fill(jr.begin(), jr.end(),0);
          ++diag_it;
        }
        
        for (diag_it = diag.begin(); diag_it != diag.end(); ++diag_it ) {
          *diag_it = 1/ *diag_it;
        }
        
        std::cout << "--< " << my_timer.elapsed() << " >-- " << ": Preconditioning done. ";
        std::cout << "Non-Zeros: " << Tri_loc << endl;

        Tri_ind.erase(Tri_ind.begin()+Tri_loc, Tri_ind.end());
        Tri_val.erase(Tri_val.begin()+Tri_loc, Tri_val.end());

        Tri = TriMatrix(A.nrows(), A.ncols(), 
                        Tri_loc, &Tri_val[0], &Tri_ptr[0], 
                        &Tri_ind[0]);
    
      } else { 
        std::cout << "Warning: It is not so efficient as symmetric row-wise upper Matrix" << std::endl;
        assert(0);
      }
    }
    void 
    do_cholesky(const SymMatrix& A, mtl::column_tag, int lfil, T tol)
    {
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
          
          A_row.filter(0); //must invoke this to prevent too many drops
          
          ++A_i;
        }
        
                
        int d, l, jj, n = A.nrows();

        std::vector<int> jj_index(n,0);
        std::vector<int> jj_pos(n,0);
        typename std::vector<int>::iterator jj_index_it;
        typename std::vector<int>::iterator jj_index_end;
        typename std::vector<int>::iterator jj_pos_it;
        typename mtl::dense1D<T>::iterator diag_it = diag.begin();
        typename SymMatrix::const_iterator j = A.begin();
        typename SymMatrix::OneD::iterator i;
        T z, t, c,c_norm;
        
        for (j = A.begin(); mtl::not_at(j,A.end()); j++) {
          typename SymMatrix::OneD j_row = *j;
          jj = j.index();
          
          wnum=0;
          i = j_row.begin();
          *diag_it = *i - *diag_it;
          i++;
          //fill work structure with non-zeros
          std::copy(i, j_row.end(), w.begin());
          std::copy(j_row.nz_struct().begin()+1, j_row.nz_struct().end(), jw.begin());
          
          for (; i != j_row.end(); ++i) {
            ++wnum;
            jr[i.index()] = wnum;
          }
          
          int nonzeros = wnum;
        
          c_norm = mtl::one_norm(*j)/ (*j).nnz()*tol;
          typename std::vector<int>::iterator tri_ptr_it = Tri_ptr.begin();
          typename mtl::dense1D<T>::iterator diag_it2 = diag.begin();
          typename std::vector<int>::iterator jr_it;
        
          jj_index_it = jj_index.begin();
          jj_index_end = jj_index_it+jj;
          jj_pos_it = jj_pos.begin();
          for (jj_index_it = jj_index.begin(); jj_index_it != jj_index_end; ++jj_index_it)
           {
            ++tri_ptr_it;
            if (*jj_index_it == jj) { // must be this element or not in row
              z = Tri_val[*jj_pos_it]* *diag_it2;
              c = ABS(c_norm/z);
              ++*jj_pos_it;
              *jj_index_it = Tri_ind[*jj_pos_it];
           
              for (l = *jj_pos_it; l < *tri_ptr_it; ++l) {
                jr_it = jr.begin()+Tri_ind[l];
                if (*jr_it) {
                  //already in work vector
                  w[(*jr_it)-1] -= Tri_val[l]*z;
                }
                else {
                  //insert this element, if large enough
                  if (ABS(Tri_val[l]) > c) {
                    w[wnum] = -Tri_val[l]*z;
                    jw[wnum] = Tri_ind[l];
                    ++wnum;
                    *jr_it = wnum;
                  }
                }
              } 
            }
            ++diag_it2;
            ++jj_pos_it;
          }
          
                  
          d = MIN(j_row.nnz()+lfil-1, wnum);
          if (d < wnum  && d > nonzeros) {
            qsplit(nonzeros, wnum, d);
          }
        
          
          typename std::vector<int>::iterator jw_it = jw.begin();
          typename std::vector<int>::iterator jw_end = jw_it+d;
          std::sort(jw_it, jw_end);
        
          Tri_loc = Tri_ptr[jj];
          typename std::vector<T>::iterator tri_val_it = Tri_val.begin()+Tri_loc;
          typename std::vector<int>::iterator tri_ind_it = Tri_ind.begin()+Tri_loc;
          
          
        
          *jj_pos_it = Tri_loc;
          *jj_index_it = *jw.begin();
        
          std::copy(jw_it, jw_end, tri_ind_it);
          for (;jw_it != jw_end; ++jw_it) {
            t = w[jr[*jw_it] -1];
            z = t/ *diag_it;
            *tri_val_it = z;
            ++tri_val_it;
            ++Tri_loc;
            diag[*jw_it] += t *z;
            
          }

          Tri_ptr[jj+1] = Tri_loc;
          
          std::fill(jr.begin(), jr.end(),0);

          ++diag_it;
        }
        
        for (diag_it = diag.begin(); diag_it != diag.end(); ++diag_it ) {
          *diag_it = 1/ *diag_it;
        }
                
        std::cout << "--< " << my_timer.elapsed() << " >-- " << ": Preconditioning done. ";
        std::cout << "Non-Zeros: " << Tri_loc << endl;
        
        Tri_ind.erase(Tri_ind.begin()+Tri_loc, Tri_ind.end());
        Tri_val.erase(Tri_val.begin()+Tri_loc, Tri_val.end());

        Tri = TriMatrix(A.nrows(), A.ncols(), 
                        Tri_loc, &Tri_val[0], &Tri_ptr[0], 
                        &Tri_ind[0]);
      }
    }
    


    
    inline void qsplit(int m, int n, int ncut) {
      int first, last, mid;
      int j,t;
      T temp, pkey;

      ncut--;
      first = m;
      last = n-1;
      if (ncut < first || ncut > last) 
        return;
      
      /* outer loop while mid != ncut */
      while (1)
        {
          mid = first;
          pkey = ABS(w[mid]);
          for (j=first+1; j<=last; j++)
            {
              if (ABS(w[j]) > pkey)
                {
                  mid = mid+1;
                 
                  /* interchange */
                  temp = w[mid];         /* swap mid and j */
                  w[mid] = w[j];
                  w[j] = temp;
                  
                  t = jr[jw[mid]];
                  jr[jw[mid]] = jr[jw[j]];
                  jr[jw[j]]= t;
                  
                  t = jw[mid];
                  jw[mid] = jw[j];
                  jw[j] = t;
                }
            }
          
          /* interchange */
          temp = w[mid];         /* swap mid and first */
          w[mid] = w[first];
          w[first] = temp;
          
          t = jr[jw[mid]];
          jr[jw[mid]] = jr[jw[first]];
          jr[jw[first]]= t;
          
          t = jw[mid];
          jw[mid] = jw[first];
          jw[first] = t;
          
          /* test for while loop */
          if (mid == ncut) 
            return;
          
          if (mid > ncut)
            last = mid-1;
          else
            first = mid+1;
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
    mtl::dense1D<T> diag;
    std::vector<T> w;
    std::vector<int> jw;
    std::vector<int> jr;
    int wnum;
    boost::timer my_timer; 
  };


}

#endif
