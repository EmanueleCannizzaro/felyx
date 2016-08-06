/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// PreProcessing.h
//
// begin     : May 21 2003
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

#ifndef CG_h
#define CG_h CG_h


namespace mtl{		// Put classes into namespace fe_base

// declaration of functions
template<class SymSkyMatrix, class myVec, class rowVec>
int myCG(const SymSkyMatrix &, myVec &, rowVec &, double , unsigned );
template<class SymSkyMatrix, class myVec, class rowVec>
int myCG(const SymSkyMatrix &, myVec &, rowVec &, double );


typedef dense1D<float_type> Vector;


//! Conjugate Gradient Algorithm
/*!
    \param A Skyline Matrix
    \param b Starting Vector
    \param x Solution Vector
    \param tol Maximum 2-norm of Residual Vector
    \param max_it Maximum Number of Iterations
*/
template<class SymSkyMatrix, class myVec, class rowVec>
int myCG(const SymSkyMatrix &A, myVec &b, rowVec &x, double tol, unsigned max_it){

  unsigned n = b.size();
  double e=0, q=0, res0=0, res=1e12;

  Vector r(n,0);
  Vector r_old(n,0);
  Vector p(n,0);
  Vector z(n,0);

  mult(A,x,r);
  add(scaled(b,-1.0),r);
  res0 = two_norm(r);
  res=res0;
  add(scaled(r,-1.0),p);

  for (unsigned it=0; it<max_it && res>tol; ++it) {

    if (it>0) {
      e = dot(r,r)/dot(r_old,r_old);
      scale(p,e);
      add(scaled(r,-1.0),p);
    }
    mult(A,p,z);
    q = dot(r,r)/dot(p,z);
    add(scaled(p,q),x);
    mtl::copy(r,r_old);
    add(scaled(z,q),r);
    res = two_norm(r);
    //cout << "it = " << it << ", res = " << res << endl;
  }
  return 0;
}


template<class SymSkyMatrix, class myVec, class rowVec>
int myCG(const SymSkyMatrix &A, myVec &b, rowVec &x, double tol){

  return myCG(A,b,x,tol,30);
}

}
#endif
