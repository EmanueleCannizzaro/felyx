// -*- c++ -*-
//-----------------------------------------------------------------------------
// mtl_felyx_utils.h
//
// begin     : June 2002
// copyright : (c) 2001 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {okoenig, wintermantel, nzehnder}@imes.mavt.ethz.ch
// www       : www.imes.ethz.ch/st
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
   

#ifndef mtl_felyx_utils_h
#define mtl_felyx_utils_h mtl_felyx_utils_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//mtl-includes
#include <mtl/matrix.h>
#include <mtl/mtl.h>
#include <mtl/dense1D.h>
#include <mtl/utils.h>
#include <mtl/lu.h>
#include "mtl-specializations/mtl_specializations.h"

using namespace std;

namespace mtl{ //put the following in to namespace mtl
  
  // Typedef for floating point values
#ifdef USE_FLOATS
  typedef float float_type;
#else
  typedef double float_type;
#endif
  
  
  //Typedef for mtl-matrices
  typedef matrix<float_type, rectangle<>, dense<>, row_major>::type Dense_Matrix;
  typedef matrix<float_type, symmetric<lower>, array<dense<> >, row_major>::type Symmetric_Dense_Matrix;
  typedef dense1D<float_type> Dense_Vector; 
  
  
  //Some functions....
  
  /*!This function takes two matrices as arguments. The first one is the matrix one 
    wants to invert. It will not be changed by the function. The second argument is the 
    inverted matrix.
  */
  template<class MatrixType>
  int lu_inversion(const MatrixType& _Matrix, MatrixType& _invMatrix){
    int size = _Matrix.nrows();
    MatrixType LU(size, size);
    dense1D<int> pvector(size);
    copy(_Matrix, LU);
    lu_factor(LU, pvector);
    lu_inverse(LU, pvector, _invMatrix);
    
    return 0;
  }
  
  template<class MatrixType>
  int inversion(const MatrixType& _Matrix, MatrixType& _invMatrix){
    int size = _Matrix.nrows();
    if (size==2)
      {
	double det = _Matrix(0,0)*_Matrix(1,1) - _Matrix(0,1)*_Matrix(1,0);

	_invMatrix(0,0) =  _Matrix(1,1);
	_invMatrix(0,1) =  -_Matrix(0,1);
	_invMatrix(1,0) =  -_Matrix(1,0);
	_invMatrix(1,1) =  _Matrix(0,0);

	scale(_invMatrix, 1.0/det);
      }
    else if (size==3)
      {
	double det = _Matrix(0,0)*_Matrix(1,1)*_Matrix(2,2) - _Matrix(0,0)*_Matrix(1,2)*_Matrix(2,1) - _Matrix(0,1)*_Matrix(1,0)*_Matrix(2,2) + _Matrix(0,2)*_Matrix(1,0)*_Matrix(2,1) + _Matrix(0,1)*_Matrix(1,2)*_Matrix(2,0) - _Matrix(0,2)*_Matrix(1,1)*_Matrix(2,0);
	
	_invMatrix(0,0) =  _Matrix(1,1)*_Matrix(2,2) - _Matrix(1,2)*_Matrix(2,1);
	_invMatrix(0,1) =  _Matrix(0,2)*_Matrix(2,1) - _Matrix(0,1)*_Matrix(2,2);
	_invMatrix(0,2) =  _Matrix(0,1)*_Matrix(1,2) - _Matrix(0,2)*_Matrix(1,1);
	_invMatrix(1,0) =  _Matrix(1,2)*_Matrix(2,0) - _Matrix(1,0)*_Matrix(2,2);
	_invMatrix(1,1) =  _Matrix(0,0)*_Matrix(2,2) - _Matrix(0,2)*_Matrix(2,0);
	_invMatrix(1,2) =  _Matrix(0,2)*_Matrix(1,0) - _Matrix(0,0)*_Matrix(1,2);
	_invMatrix(2,0) =  _Matrix(1,0)*_Matrix(2,1) - _Matrix(1,1)*_Matrix(2,0);
	_invMatrix(2,1) =  _Matrix(0,1)*_Matrix(2,0) - _Matrix(0,0)*_Matrix(2,1);
	_invMatrix(2,2) =  _Matrix(0,0)*_Matrix(1,1) - _Matrix(0,1)*_Matrix(1,0);
	
	scale(_invMatrix, 1.0/det);
      }
    else
      {
	cerr << endl << "#####################################################";
	cerr << endl << "ERROR: in inversion. Matrix dimension must agree. ";
	cerr << endl << "#####################################################";
      }
    return 0;
  }

  /*!This function calculates Res = A^T *B * A
    Assuming the Matrix A is not square, but the Matrix B is.
    The function takes 3 arguments like a normal multiplication.
    The Matrices A, B, and the resulting matrix Res.
  */
  template<class MatrixType>
  int mult_AT_B_A( const MatrixType A, const MatrixType B, MatrixType Res ){
    
    MatrixType Res1( B.nrows(), A.ncols() );
    //MatrixType AT( A.ncols(), A.nrows() );

    //mtl::set_value( Res1, float_type(0) );
    //mtl::set_value( AT, float_type(0) );
    mtl::set_value( Res, float_type(0) );

    mtl::mult( B, A, Res1 );
    //transpose(A, AT );
    mtl::mult( trans(A), Res1, Res );

    return 0;
  }

  template<class MatrixType>
  int mult_AT_B_A_add( const MatrixType A, const MatrixType B, MatrixType Res ){
    
    MatrixType Res1( B.nrows(), A.ncols() );

    mtl::mult( B, A, Res1 );
    mtl::mult(trans(A), Res1, Res );

    return 0;
  }
  
  /*!This function calculates the determinant of a matrix
   */
  template<class MatrixType>
    double lu_det(const MatrixType& _Matrix){
    int size = _Matrix.nrows();
    MatrixType LU(size,size);
    dense1D<int> pvector(size);
    copy(_Matrix, LU);
    lu_factor(LU, pvector);
    
    //For every permutation of rows the product of 
    //diagonal elements has to be multyplied by -1
    int sign = 1;
    for (int i = 0 ; i < size ; ++i){
      if ( pvector[i] != i+1 ) {
	sign *= -1; 
      }
    }
    
    typename MatrixType::iterator iter;
    iter = LU.begin();
    
    double det = 1.0;
    
    for ( iter = LU.begin(); iter != LU.end(); ++iter ){
      unsigned ind = iter.index();     
      det *= LU(ind, ind);
    }

    det *= sign;
    return det;
  }

  template<class MatrixType>
  double determinant(const MatrixType& _Matrix){
    int size = _Matrix.nrows();
    double det=0;
    if (size==2)
      {
	det = _Matrix(0,0)*_Matrix(1,1) - _Matrix(0,1)*_Matrix(1,0);
      }
    else if (size==3)
      {
	det = _Matrix(0,0)*_Matrix(1,1)*_Matrix(2,2) - _Matrix(0,0)*_Matrix(1,2)*_Matrix(2,1) - _Matrix(0,1)*_Matrix(1,0)*_Matrix(2,2) + _Matrix(0,2)*_Matrix(1,0)*_Matrix(2,1) + _Matrix(0,1)*_Matrix(1,2)*_Matrix(2,0) - _Matrix(0,2)*_Matrix(1,1)*_Matrix(2,0);
      }
    else
      {
	cerr << endl << "#####################################################";
	cerr << endl << "ERROR: in determinant. Matrix dimension must agree. ";
	cerr << endl << "#####################################################";
      }
    return det;
  }

  //!The cross product for 3 dimensional vectors
   template<class VecX, class VecY, class VecZ>
   //template<class VecX>
  int cross_prod_3d(const VecX X, const VecY Y, VecZ Z){
     //int cross_prod_3d( VecX X,  VecX Y, VecX Z){

    MTL_ASSERT( X.size() == 3 && Y.size() == 3 && Z.size() == 3 , "mtl::cross_prod_3d()");
    
    Z[0] = X[1]*Y[2] - X[2]*Y[1];
    Z[1] = X[2]*Y[0] - X[0]*Y[2];
    Z[2] = X[0]*Y[1] - X[1]*Y[0];

    return 0;
    
  }

  //!The scalar product for 3 dimensional vectors
  template<class VectorType>
  double scalar_prod_3d(const VectorType& _a, const VectorType& _b){
    return (_a[0]*_b[0] + _a[1]*_b[1] + _a[2]*_b[2]);
  }


  //! This function checks, if a matrix is symmetric or not.
  /*! Symmetry is granted, if two elements K(a,b) and K(b,a) have a
    relative error of less than epsError and are both smaller
    than epsValue. The noiselevel toggles the message of symmetry on
    and off
  */
  template<class MatrixType>
  bool symmetryCheck(MatrixType& _myMatr, int _noiselevel = 1, double _epsError = 1e-10, double _epsValue = 1e-14){
    if ( _myMatr.nrows() != _myMatr.ncols() ){
      cout << endl << "Error in mtl::symmetryCheck: Matrix is not square";
      exit(0);
    }

    double relErr = 0.0;
    int symcheck = 0;
    for (int m = 0 ; m < _myMatr.nrows() ; ++m ){
      for (int n = 0 ; n < m ; ++n ){
	if (_myMatr(m,n) != _myMatr(n,m)){
	  relErr = (_myMatr(m,n) - _myMatr(n,m)) / _myMatr(n,m);
	  if ( abs(relErr) > _epsError && _myMatr(m,n) > _epsValue && _myMatr(n,m) > _epsValue){
	    symcheck = 1;
	    cout << endl << "***** The matrix is NOT symmetric within the following range!!!" << endl
		 <<         "***** The relative Error of element (" << n << "," << m << ") is: " << relErr << endl
		 <<         "***** The values are " << _myMatr(m,n) << "\t" << _myMatr(n,m) << endl
		 <<         "***** The threshold for the relative Error is: " << _epsError << endl
		 <<         "***** The threshold for the matrix elements is: " << _epsValue << endl;
	  }
	}  
      }
    }
    if (_noiselevel != 0 && symcheck == 0)
    cout << endl
	 << "***** The matrix is symmetric within the following range!!! " << endl
	 << "***** The threshold for the relative error is: " << _epsError << endl
	 << "***** The threshold for the elements values is: " << _epsValue << endl;

    return !symcheck;
  }
  //!Checking the orthogonality of a matrix
  template<class MatrixType>
  bool orthogonalityCheck( MatrixType& _myMatr, int _noiselevel = 1 ){

    int rows = _myMatr.nrows(), cols = _myMatr.ncols();

    bool ret = 0;

    if ( rows == cols ){
      //some matrices
      Dense_Matrix OrthoTest(rows,cols), resultMatrix(rows,cols), Identity(rows,cols);
      mtl::set_value(Identity,0.0);
      mtl::set_diagonal(Identity,-1.0);
      copy(_myMatr, OrthoTest);

      transpose(OrthoTest);
      mult(_myMatr,OrthoTest,resultMatrix);
      add(Identity,resultMatrix);
      
      ret = zeroCheck( resultMatrix , 0 );

      if ( ret && _noiselevel != 0 )
	cout << endl << " The matrix is orthogonal " << endl;
      if ( !ret && _noiselevel != 0 )
	cout << endl << " The matrix is not orthogonal " << endl;
    }
    else {
      cerr << endl << "ERROR in orthogonalityCheck() : The matrix is not square" << endl;
      ret = 0;
    }
    return ret;
  }

  //!Checking wheter a matrix is equal to zero within a certain threshold
  template<class MatrixType>
  bool zeroCheck( MatrixType& _myMatr, int _noiselevel = 1 , double _limit = 1e-14){
    int count = 0;
    unsigned i,k;
    for ( i = 0 ; i < _myMatr.nrows() ; ++i ){
      for ( k = 0 ; k < _myMatr.ncols() ; ++k ){
	if ( _myMatr(i,k) > _limit )
	  ++count;
      }
    }
    if ( count == 0 && _noiselevel != 0 )
      cout << endl << " The matrix is zero within the limit " << _limit << endl; 
    if ( count != 0 && _noiselevel != 0 )
      cout << endl << " The matrix is not zero "
	   << endl << " there are " << count << " elements greater than " << _limit << endl;

    return !count;
  }

  //!Checking wheter two matrices are equal within a certain limit
  template<class MatrixType>
  bool MatrixEquality( MatrixType& _a, MatrixType& _b, int _noiselevel = 1, double _limit = 1e-4){
    if (_a.nrows() != _b.nrows() || _a.ncols() != _b.ncols())
      cerr << endl << "ERROR in MatrixEquality: The matrices have different sizes" << endl;
    double err = 0.0, dummy = 0.0;
    unsigned i,k;
    bool res = true;
    for ( i = 0 ; i < _a.nrows() ; ++i ){
      for ( k = 0 ; k < _a.ncols() ; ++k ){
	dummy = _a(i,k) - _b(i,k);
	if ( abs(dummy) > err && abs(dummy) > _limit)
	  err = dummy;
      }
    }

    if ( err < _limit && _noiselevel == 1 )
      cout << endl << " The matrices are equal within the limit " << _limit << endl; 

    if ( err > _limit && _noiselevel == 1 )
      cout << endl << " The matrices are NOT equal within the limit " << _limit << endl; 

    if (err > _limit && _noiselevel == 2 ){
      cout << endl << " The matrices are NOT equal within the limit " << _limit 
	   << endl << " The matrices are: " << endl;
      print_all_matrix(_a);
      print_all_matrix(_b);
    }

    if ( err > _limit ) res = false;

    return res;

  }

template<class VecType>
void printv(VecType a) {
  cout << "[" << endl;
  cout.precision(6);
  for (unsigned i=0;i<a.size();i++)
    {
  cout << "[" << a[i] << "]" << endl;
    }
  cout << "]" << endl;
}

template<class MatrixType>
void print_envelope2dense(MatrixType& _M ){
  typename MatrixType::iterator 		matrixit;
  typename MatrixType::OneD::iterator 	oneDit;
  cout.setf(ios::right, ios::adjustfield );
  cout.precision(6);
  Dense_Matrix D(_M.nrows(),_M.ncols());
  unsigned i(0),j(0);
  for(matrixit = _M.begin(); matrixit != _M.end(); ++matrixit){ 
    oneDit = (*matrixit).begin();
    i=0;
    while ( i < oneDit.index() ){
      ++i;
    }
    // print vals of matrix
    for (; oneDit != (*matrixit).end(); ++oneDit){
      D(j,i)=*oneDit; i++;
    }
    j++;
  }
  print_all_matrix(D);
}

} //end of namespace

#endif //mtl_felyx_utils_h

