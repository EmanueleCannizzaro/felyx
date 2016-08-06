//-----------------------------------------------------------------------------
// SkylineSolver.h
//
// begin     : Mai 10 2002
// copyright : (c) 2001 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {okoenig, wintermantel, nzehnder}@imes.mavt.ethz.ch
// www       : www.structures.ethz.ch
/* 
   This file is part of FELyX (Finite Element Library eXperiment).
   
   FELyX is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   FELyX is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with FELyX; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
//-----------------------------------------------------------------------------

#ifndef SKYLINE_SOLVER_H
#define SKYLINE_SOLVER_H SKYLINE_SOLVER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "BlasHeaders.h"

// Declaration of functions
namespace mtl {

  template<class SymSkyMatrix, class rowVec>
  int sky_decomposition(SymSkyMatrix&, rowVec&);
  
  template<class SymSkyMatrix, class myVec, class rowVec>
  void sky_backsubstitution(SymSkyMatrix &, myVec &, rowVec &);
  
  template<class SymSkyMatrix, class myVec, class rowVec>
  int skyline_solve(SymSkyMatrix &, myVec &, rowVec &);
  
}

////
// Implementation of skyline solver, using atlas dot product
////
#ifdef HAVE_BLAS

namespace mtl {
  
  template<class SymSkyMatrix, class rowVec>
  inline int sky_decomposition(SymSkyMatrix &A, rowVec &ch_)
  {
    
    typename SymSkyMatrix::iterator matrixit1, matrixit2;
    typename SymSkyMatrix::OneD::iterator oneDit1, oneDit5;
    typename rowVec::iterator rowit1, rowit2, rowit3;
    typedef typename SymSkyMatrix::value_type valueT;
    typedef typename SymSkyMatrix::size_type sizeT;
    
    rowVec DI(ch_.size());  
    sizeT i, j, dij;
    valueT s, pivot;
    int err_msg = 0;
    
    ////////////////////////////////////////////////////////
    //Initializing distance vector
    //      
      
    i = 1;
    rowit2 = ch_.begin();
    for(rowit1 = DI.begin(); rowit1 != DI.end(); ++rowit1, ++rowit2, ++i)
      *rowit1 = i - *rowit2;

      
    ////////////////////////////////////////////////////////////
    //Factorization
    //
      
    //loop 1
    rowit2 = DI.begin();
      
    for( matrixit1 = A.begin(); matrixit1 != A.end(); ++matrixit1, ++rowit2)
      { 
        //loop 2
        oneDit1 = (*matrixit1).begin();
        j = oneDit1.index();
        rowit3 = DI.begin() + j;
        oneDit5 = ( (*matrixit1).end() - 1 );
        for(; oneDit1 != oneDit5; ++oneDit1, ++rowit3, ++j )
          {
            dij = std::max(*rowit2, *rowit3);

            // Calling the cblas dot product, either the routin for floats or for doubles
#ifdef USE_FLOATS
            (*oneDit1) -= cblas_sdot( j - dij,
                                      &(*( (*matrixit1).begin() + dij - *rowit2 ) ),
                                      1,
                                      &(* (  (*(A.begin() + j)).begin() + dij - *rowit3 ) ),
                                      1 );
#else
            (*oneDit1) -= cblas_ddot( j - dij,
                                      &( *( (*matrixit1).begin() + dij - *rowit2 ) ),
                                      1,
                                      &( *(  (*(A.begin() + j)).begin() + dij - *rowit3 ) ),
                                      1 );
#endif      
          }//end loop 2
        
        pivot = *( (*matrixit1).end() - 1); 
        
        oneDit1 = (*matrixit1).begin(); 
        matrixit2 = A.begin() + oneDit1.index();
        rowit1 = ch_.begin() + oneDit1.index();
        
        //loop 3 
        oneDit5 = (*matrixit1).end() - 1;    
        for(; oneDit1 != oneDit5; ++oneDit1, ++matrixit2, ++rowit1 ) 
          {
            
            s =  (*oneDit1) * ( *( (*matrixit2).begin() + (*rowit1) - 1 ) );
            pivot -= s * (*oneDit1);
            *oneDit1 = s;
            
          }//end loop 3
        
        if(abs(pivot) < 1.0e-15)
          {
            err_msg = 1;
            break;
          } 
        
        *( (*matrixit1).end() - 1) = 1.0 / pivot;
        
      }//end loop 1
     
      return err_msg;
  }

  template<class SymSkyMatrix, class myVec, class rowVec>
void sky_backsubstitution(SymSkyMatrix &A, myVec &b, rowVec &ch_)
        {
    typename SymSkyMatrix::iterator matrixit1;
    typename SymSkyMatrix::reverse_iterator matrixrit1;
    typename SymSkyMatrix::OneD::iterator oneDit1, oneDit5;
    typename rowVec::iterator rowit1, rowit2, rowit3;
    typename myVec::iterator myit1, myit2;
    typedef typename SymSkyMatrix::value_type valueT;
    typedef typename SymSkyMatrix::size_type sizeT;
    
    rowVec DI(ch_.size());  
    sizeT i;
    
    ////////////////////////////////////////////////////////
    //Initializing distance vector
    //      
      
    i = 1;
    rowit2 = ch_.begin();
    for(rowit1 = DI.begin(); rowit1 != DI.end(); ++rowit1, ++rowit2, ++i)
      *rowit1 = i - *rowit2;      

          myit2 = b.begin();
          rowit2 = DI.begin();
          for(matrixit1 = A.begin(); matrixit1 != A.end(); ++matrixit1, ++myit2, ++rowit2 )
            {
#ifdef USE_FLOATS
              *myit2 -= cblas_sdot( ch_[ matrixit1.index() ] -1,
                              &(*( (*matrixit1).begin() )) , 
                              1, 
                              &(*( b.begin() + (*rowit2) ) ),           
                                    1);
#else
              *myit2 -= cblas_ddot( ch_[ matrixit1.index() ] -1,
                              &(*( (*matrixit1).begin() )) , 
                              1, 
                              &(*( b.begin() + (*rowit2) ) ),           
                                    1);
#endif      
            }
          
          matrixit1 = A.begin();
          for(myit1 = b.begin(); myit1 != b.end(); ++myit1, ++matrixit1)
            {
              (*myit1) *= *( (*matrixit1).end() - 1);
            }
          
          myit1 = b.end() - 1;
          rowit1 = ch_.end() - 1;
          rowit2 = DI.end() - 1;
          for(matrixrit1 = A.rbegin(); matrixrit1 != A.rend(); ++matrixrit1, --myit1, --rowit1, --rowit2)
            {
              myit2 = b.begin();
              oneDit5 = (*matrixrit1).end()-1;
              for( oneDit1 = (*matrixrit1).begin(); oneDit1 != oneDit5; ++oneDit1, ++myit2 )
                *( myit2 + *rowit2 )  -= (*myit1) * (*oneDit1);
            }
          
        }
  
  template<class SymSkyMatrix, class myVec, class rowVec>
  int skyline_solve(SymSkyMatrix &A, myVec &b, rowVec &ch_)
    {
      int err_msg = 0;
      err_msg = sky_decomposition(A, ch_);
      sky_backsubstitution(A, b, ch_);
      return err_msg;
    }

} // of namespace

#else

////
// Implementation of Skyline Solver without BLAIS functionality
////

namespace mtl {
    
  template<class SymSkyMatrix, class rowVec>
  int sky_decomposition(SymSkyMatrix &A, rowVec &ch_)
  {
    
    typename SymSkyMatrix::iterator matrixit1, matrixit2;
    typename SymSkyMatrix::OneD::iterator oneDit1, oneDit2, oneDit3, oneDit4, oneDit5;
    typename rowVec::iterator rowit1, rowit2, rowit3;
    typedef typename SymSkyMatrix::value_type valueT;
    typedef typename SymSkyMatrix::size_type sizeT;
    
    rowVec DI(ch_.size());  
    sizeT i, j, dij;
    valueT s, pivot;
    int err_msg = 0;
    
    /////////////////////////////////////////////////////////////
      //Initializing distance vector
      //
      
      i = 1;
      rowit2 = ch_.begin();
      for(rowit1 = DI.begin(); rowit1 != DI.end(); ++rowit1)
        {
          *rowit1 = i - *rowit2;
          
          ++ rowit2;
      ++i;
        }
      
      
      
      ////////////////////////////////////////////////////////////
      //Factorization
      //
      
      //loop 1
      rowit2 = DI.begin();
      
      for(matrixit1 = A.begin(); matrixit1 != A.end(); ++matrixit1)
        { 
          //loop 2
          oneDit1 = (*matrixit1).begin();
          j = oneDit1.index();
          rowit3 = DI.begin() + j;
          oneDit5 = ( (*matrixit1).end() - 1 );
          for(; oneDit1 != oneDit5; ++oneDit1)
            {
              dij = std::max(*rowit2, *rowit3);
              
              oneDit3 = (*(A.begin() + j)).begin() + dij - *rowit3;
              oneDit4 = ((*matrixit1).begin() + j - (*rowit2));

              for( oneDit2 = (*matrixit1).begin() + dij - *rowit2; oneDit2 != oneDit4; ++oneDit2 ,++oneDit3 ) 
                {
                  (*oneDit1) -= (*oneDit2) * (*oneDit3);
                  
                }//end subloop
              
              ++rowit3;
              ++j;
            }//end loop 2
          
          pivot = *( (*matrixit1).end() - 1); 
      
          oneDit1 = (*matrixit1).begin(); 
          matrixit2 = A.begin() + oneDit1.index();
          rowit1 = ch_.begin() + oneDit1.index();
          
          //loop 3 
          oneDit5 = (*matrixit1).end() - 1;    
          for(; oneDit1 != oneDit5; ++oneDit1 ) 
            {
              s =  (*oneDit1) * ( *( (*matrixit2).begin() + (*rowit1) - 1 ) );
              pivot -= s * (*oneDit1);
              *oneDit1 = s;
              
              ++matrixit2;
              ++rowit1;
            }//end loop 3
          
          if(abs(pivot) < 1.0e-15)
            {
              err_msg = 1;
              break;
            } 
          *( (*matrixit1).end() - 1) = 1.0 / pivot;
          ++rowit2;       
        }//end loop 1
           
      return err_msg;
  }
  
template<class SymSkyMatrix, class myVec, class rowVec>
void sky_backsubstitution(SymSkyMatrix &A, myVec &b, rowVec &ch_)
{
    typename SymSkyMatrix::iterator matrixit1;
    typename SymSkyMatrix::reverse_iterator matrixrit1;
    typename SymSkyMatrix::OneD::iterator oneDit1, oneDit2, oneDit3, oneDit4, oneDit5;
    typename rowVec::iterator rowit1, rowit2, rowit3;
    typename myVec::iterator myit1, myit2;
    typedef typename SymSkyMatrix::value_type valueT;
    typedef typename SymSkyMatrix::size_type sizeT;
    
    rowVec DI(ch_.size());  
    sizeT i;
    valueT s;
        
    /////////////////////////////////////////////////////////////
      //Initializing distance vector
      //
      
      i = 1;
      rowit2 = ch_.begin();
      for(rowit1 = DI.begin(); rowit1 != DI.end(); ++rowit1)
        {
          *rowit1 = i - *rowit2;
          
          ++ rowit2;
      ++i;
        }
      
      
      
      ////////////////////////////////////////////////////////////

          myit2 = b.begin();
          rowit2 = DI.begin();
          for(matrixit1 = A.begin(); matrixit1 != A.end(); ++matrixit1)
            {
              s = 0;
              
              myit1 = b.begin() + (*rowit2);
              oneDit5 = (*matrixit1).end() - 1;
              for(oneDit1 = (*matrixit1).begin(); oneDit1 != oneDit5; ++oneDit1)
                {
                  s += (*oneDit1) * (*myit1);
                  
                  ++myit1;
                }
              
              *myit2 -= s;
              
              ++myit2;
              ++rowit2;
            }
          
          matrixit1 = A.begin();
          for(myit1 = b.begin(); myit1 != b.end(); ++myit1)
            {
              (*myit1) *= *( (*matrixit1).end() - 1);
              
              ++matrixit1;
            }
          
          myit1 = b.end() - 1;
          rowit1 = ch_.end() - 1;
          rowit2 = DI.end() - 1;
          for(matrixrit1 = A.rbegin(); matrixrit1 != A.rend(); ++matrixrit1)
            {
              myit2 = b.begin();
              oneDit5 = (*matrixrit1).end()-1;
              for( oneDit1 = (*matrixrit1).begin(); oneDit1 != oneDit5; ++oneDit1 )
                {
                  *( myit2 + *rowit2 )  -= (*myit1) * (*oneDit1);
                  
                  ++myit2;
                }
              
              --myit1;
              --rowit1;
              --rowit2;
            }
          
        }

  template<class SymSkyMatrix, class myVec, class rowVec>
  int skyline_solve(SymSkyMatrix &A, myVec &b, rowVec &ch_)
    {
      int err_msg = 0;
      err_msg = sky_decomposition(A, ch_);
      sky_backsubstitution(A, b, ch_);
      return err_msg;
    }

} // of namespace

#endif // of #ifdef HAVE_BLAS

#endif
