//-----------------------------------------------------------------------------
// BaseSolid.cc
//
// begin     : Nov 20 2001
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
   

#include "BaseSolid.h"


using namespace fe_base;

Dense_Matrix BaseSolid::CalcEM() {

  //Initialization of all variables
  double detJ, fdummy;
  
  //The size of the square ESM is given as follows
  unsigned  nSize = GetEMSize();
  unsigned dim = GetDofSet().ElementDimension();

  //Initialization of some Matrices 
  Dense_Matrix 
    J(dim,dim),                     /* Jacobimatrix */
    B((dim-1)*3,nSize),                     /* B-Matrix */
    K(nSize,nSize),                 /* Stifnessmatrix */
    dummyMatr1(6,nSize),            /* Dummymatrix */
    dummyMatr2(nSize,nSize),
    shapefunc(GetNodeCount(), GetDofSet().count()), /* derived - shapefuncs */
    IntPoints = GetIntPoints(),
    D, D_copy;  

  //Get the material data
  //checking wether the elements are planes ore solids
  //Switching between plane-strain and plain-stress has
  //to be done here.
  if ( dim == 3 ) {
    D.resize(6,6); D_copy.resize(6,6);
    mtl::set_value(D,0.0);
    copy(MaterialPtr -> Get(Material::ThreeDSolid), D);
  }
  else if ( dim == 2 ) {
    D.resize(3,3); D_copy.resize(3,3);
    mtl::set_value(D,0.0);
    copy(MaterialPtr -> Get(Material::PlaneStress12), D);
  }  
  else {
    cerr << "\nERROR in Element::CalcEM()" << endl;
    cerr << "DofSet::ElementDimension returned " << endl;
    cerr << "something else than 2 or 3 " << endl;
    exit(0);
    }

  //Calculate the function at every Gausspoint.
  for (unsigned i = 0 ; i < IntPoints.nrows(); ++i) {
    //Derived shape functions N with respect to  xi, eta and zeta coordinates
    //computed at the point represented by the row i of the matrix
    //IntPoints. They are stored in shapefunc.
    EvalDerivedShapeFunc(shapefunc, IntPoints, i);

    
    //Compute the Jacobian at the ith integration point
    EvalJacobi(J, shapefunc, IntPoints, i);

    
    //Transforme the shapefunctions to global x, y, z
    //coordinates
    stCoord2globalCoord(J, shapefunc);

   
    //----------------------------------------------------------
    //Choose one of the following functions for the integration
    //on the element (ansys uses full integration on solid186)
    //To use the bBar method with this element, the method has to 
    //be adjusted (different bBar therm calculation...)
    //----------------------------------------------------------
    SetBMatrix(shapefunc,B,FullIntegration);

    //----------------------------------------------------
    //Computation of the integrand at the ith integration
    //point and addition to the stifnessmatrix
    //----------------------------------------------------
    detJ = determinant(J);
     
    //IntPoints(i,3) is the weight of the ith point.
    //It is included for to be able to change the 
    //degree of integration easely


    fdummy = detJ*IntPoints(i,IntPoints.ncols()-1); 

    mtl::copy(D, D_copy);

    scale(D_copy,fdummy);
    mult_AT_B_A_add(B,D_copy,K);
  }

  //  cout << endl << "Checking the ESM for Symmetry";
  //  symmetryCheck(K,1,1e-3,1e-9);

  return K; //The stifnessmatrix of the element....


} // of CalcEM

Dense_Matrix BaseSolid::CalcEmm() {

  //Initialization of all variables
  double detJ, fdummy;

  //The size of the square ESM is given as follows
  unsigned  nSize = GetEMSize();
  unsigned dim = GetDofSet().ElementDimension();

  //Initialization of some Matrices 
  Dense_Matrix 
    J(dim,dim),                     /* Jacobimatrix */
    N(dim,nSize),                     /* B-Matrix */
    M(nSize,nSize),                 /* Stifnessmatrix */
    shapefunc(GetNodeCount(), 1), /* derived - shapefuncs */
    derivedshapefunc(GetNodeCount(), GetDofSet().count()),
    IntPoints = GetIntPoints();

  //Calculate the function at every Gausspoint.
  for (unsigned i = 0 ; i < IntPoints.nrows(); ++i) {

    EvalShapeFunc(shapefunc, IntPoints, i);
    EvalDerivedShapeFunc(derivedshapefunc, IntPoints, i);
     
    //Compute the Jacobian at the ith integration point
    EvalJacobi(J, derivedshapefunc, IntPoints, i);
    
    SetNMatrix(shapefunc,N);

    detJ = determinant(J);

    fdummy = detJ*IntPoints(i,IntPoints.ncols()-1)*MaterialPtr->Get("rho"); 
    
    mult(trans(N),scaled(N,fdummy),M);

  }

  return M; 
} // of CalcEMM


//compute the Jacobian at a specific point given as the row IP of a matrix(..,3)
//depending on the displacement dof in z direction different algorithms for 2D or 3D are used
//
void BaseSolid::EvalJacobi(Dense_Matrix J_, Dense_Matrix shapefunc_, Dense_Matrix EvalPoints_, unsigned IP_){

  //initialize the matrix to zero   
  mtl::set_value(J_,0.0);

  if (GetDofSet().IsElement3DSolid()){

    for (unsigned k = 0 ; k < GetNodeCount() ; ++k){
      J_(0,0) += shapefunc_(k,0)*NodeVec[k]->Cx;
      J_(0,1) += shapefunc_(k,0)*NodeVec[k]->Cy;
      J_(0,2) += shapefunc_(k,0)*NodeVec[k]->Cz;
      J_(1,0) += shapefunc_(k,1)*NodeVec[k]->Cx;
      J_(1,1) += shapefunc_(k,1)*NodeVec[k]->Cy;
      J_(1,2) += shapefunc_(k,1)*NodeVec[k]->Cz;
      J_(2,0) += shapefunc_(k,2)*NodeVec[k]->Cx;
      J_(2,1) += shapefunc_(k,2)*NodeVec[k]->Cy;
      J_(2,2) += shapefunc_(k,2)*NodeVec[k]->Cz;
    }
  }
  else {
    for (unsigned k = 0 ; k < GetNodeCount() ; ++k)
      {
	J_(0,0) += shapefunc_(k,0)*NodeVec[k]->Cx;
	J_(0,1) += shapefunc_(k,0)*NodeVec[k]->Cy;
	J_(1,0) += shapefunc_(k,1)*NodeVec[k]->Cx;
	J_(1,1) += shapefunc_(k,1)*NodeVec[k]->Cy;
      }
  }
}



void BaseSolid::stCoord2globalCoord(const Dense_Matrix Jacobi_,  Dense_Matrix shapefunc_){

  //inversion of the Jacobian
  unsigned size = Jacobi_.nrows();
  Dense_Matrix invJacobi(size, size);
  inversion(Jacobi_, invJacobi);

  if (GetDofSet().IsElement3DSolid()){
    Dense_Vector dummyvec1(size), dummyvec2(size);
    for (unsigned k = 0 ; k < GetNodeCount() ; k++)
      {
	dummyvec1[0] = shapefunc_(k,0);
	dummyvec1[1] = shapefunc_(k,1);
	dummyvec1[2] = shapefunc_(k,2);
	

	//dummyvec2 = invJacobi*dummyvec1;
	mtl::set_value(dummyvec2,0.0);
	mult(invJacobi, dummyvec1, dummyvec2);
	
	shapefunc_(k,0) = dummyvec2[0];
	shapefunc_(k,1) = dummyvec2[1];
	shapefunc_(k,2) = dummyvec2[2];
      }
  }
  else{
    Dense_Vector dummyvec1(size), dummyvec2(size);
    for (unsigned k = 0 ; k < GetNodeCount() ; k++)
      {
	dummyvec1[0] = shapefunc_(k,0);
	dummyvec1[1] = shapefunc_(k,1);
	
	//dummyvec2 = invJacobi*dummyvec1;
	mtl::set_value(dummyvec2,0.0);
	mult(invJacobi, dummyvec1, dummyvec2);

	
	shapefunc_(k,0) = dummyvec2[0];
	shapefunc_(k,1) = dummyvec2[1];
      }
  }
}


void BaseSolid::SetBMatrix( const Dense_Matrix D, Dense_Matrix B, BMatrixTypes type ){

  mtl::set_value(B,0.0);

  switch (type)
    {
    case FullIntegration :
      {
	if (GetDofSet().IsElement3DSolid())
	  {
	    for (unsigned k = 0 ; k < GetNodeCount() ; k++)
	      {
		B(0,(k+1)*3-3) = D(k,0);
		B(1,(k+1)*3-2) = D(k,1);
		B(2,(k+1)*3-1) = D(k,2);
		B(3,(k+1)*3-3) = D(k,1);
		B(3,(k+1)*3-2) = D(k,0);
		B(4,(k+1)*3-2) = D(k,2);
		B(4,(k+1)*3-1) = D(k,1);
		B(5,(k+1)*3-3) = D(k,2);
		B(5,(k+1)*3-1) = D(k,0);
	      }
	  }
	else
	  {
	    for (unsigned k = 0 ; k < GetNodeCount() ; k++)
	      {
		B(0,(k+1)*2-2) = D(k,0);
		B(1,(k+1)*2-1) = D(k,1);
		B(2,(k+1)*2-2) = D(k,1);
		B(2,(k+1)*2-1) = D(k,0);
	      }
	  }
	break;	
      }
    case bBarMethod :
      {
	if (GetDofSet().IsElement3DSolid())
	  {
	    //Set the dilatational Part of the differential equation.
	    //presently it is computet in for the origin (constant over one element)
	    Dense_Matrix IntegrPoint(1,3), J(3,3), bBar(GetNodeCount(),GetDofSet().count());
	    
	    EvalDerivedShapeFunc(bBar, IntegrPoint, 0);
	    EvalJacobi(J,bBar, IntegrPoint, 0);
	    stCoord2globalCoord(J, bBar);
	    
	    double a = 2.0/3.0, b = 1.0/3.0;
	    
	    for (unsigned k = 0 ; k < GetNodeCount() ; k++)
	      {
		B(0,(k+1)*3-3) = a*D(k,0)+b*bBar(k,0);
		B(0,(k+1)*3-2) = b*bBar(k,1)-b*D(k,1);
		B(0,(k+1)*3-1) = b*bBar(k,2)-b*D(k,2);
		B(1,(k+1)*3-3) = b*bBar(k,0)-b*D(k,0);
		B(1,(k+1)*3-2) = a*D(k,1)+b*bBar(k,1);
		B(1,(k+1)*3-1) = b*bBar(k,2)-b*D(k,2);
		B(2,(k+1)*3-3) = b*bBar(k,0)-b*D(k,0);
		B(2,(k+1)*3-2) = b*bBar(k,1)-b*D(k,1);
		B(2,(k+1)*3-1) = a*D(k,2)+b*bBar(k,2);
		B(3,(k+1)*3-3) = D(k,1);
		B(3,(k+1)*3-2) = D(k,0);
		B(4,(k+1)*3-2) = D(k,2);
		B(4,(k+1)*3-1) = D(k,1);
		B(5,(k+1)*3-3) = D(k,2);
		B(5,(k+1)*3-1) = D(k,0);
	      }
	  }
	else
	  {
	    //Set the dilatational Part of the differential equation.
	    //presently it is computet for the origin (constant over one element)
	    Dense_Matrix IntegrPoint(1,2), J(2,2), bBar(GetNodeCount(),GetDofSet().count());
	    
	    EvalDerivedShapeFunc(bBar, IntegrPoint, 0);
	    EvalJacobi(J, bBar, IntegrPoint, 0);
	    stCoord2globalCoord(J, bBar);
	    
	    double a = 2.0/3.0, b = 1.0/3.0;
	    
	    for (unsigned k = 0 ; k < GetNodeCount() ; k++)
	      {
		B(0,(k+1)*2-2) = a*D(k,0)+b*bBar(k,0);
		B(0,(k+1)*2-1) = b*bBar(k,1)-b*D(k,1);
		B(1,(k+1)*2-2) = b*bBar(k,0)-b*D(k,0);
		B(1,(k+1)*2-1) = a*D(k,1)+b*bBar (k,1);
		B(2,(k+1)*2-2) = D(k,1);
		B(2,(k+1)*2-1) = D(k,0);
	      }
	  }
	break;
      }
    case ShellIntegration :
      {
	cerr << endl << "####         Error in BaseSolid::SetBMatrix          ####"
	     << endl << "#### The Element calling the function is not a Shell! ####" 
	     << endl << "#### It is a " << GetName() <<  " number " << GetId() << " ##### " << endl;
	
	break;
      }
    }
}


////
// Eval mass of an element, based on volume calculation
////
double BaseSolid::EvalMass() const{
  FELYX_RUNTIME_ASSERT( MaterialPtr != NULL,"BaseSolid::EvalMass()");
  return EvalVolume() * MaterialPtr->Get("rho");
}


void BaseSolid::SetNMatrix( const Dense_Matrix S, Dense_Matrix N ){

  mtl::set_value(N,0.0);
  if (GetDofSet().IsElement3DSolid())
    {
      for (unsigned k = 0 ; k < GetNodeCount() ; k++)
	{
	  N(0,(k+1)*3-3) = S(k,0);
	  N(1,(k+1)*3-2) = S(k,0);
	  N(2,(k+1)*3-1) = S(k,0);
	}
    }
  else
    {
      for (unsigned k = 0 ; k < GetNodeCount() ; k++)
	{
	  N(0,(k+1)*2-2) = S(k,0);
	  N(1,(k+1)*2-1) = S(k,0);
	}
    }
}


////
// Eval nodal stresses and store'em in vector< dense_vector > of Element
////
void BaseSolid::EvalStresses(string StressType) {

  unsigned gnc = GetNodeCount();
  unsigned gds = GetDofSet().count();
  unsigned dof = gnc*gds;
  unsigned nSize = GetEMSize();
  unsigned dim = GetDofSet().ElementDimension();
  Dense_Vector u(dof);
  Dense_Vector epsilon((dim-1)*3);
  Dense_Matrix B((dim-1)*3,nSize);
  Dense_Matrix shapefunc(GetNodeCount(), GetDofSet().count());
  Dense_Matrix IntPoints = GetIntPoints();
  Dense_Matrix D,     J(dim,dim);
  Dense_Matrix sigma_ip(IntPoints.nrows(),(dim-1)*3);
  Dense_Matrix sigma_nodes(gnc,(dim-1)*3);
  unsigned nip = IntPoints.nrows();
  Dense_Vector p(nip);
  Dense_Matrix dummyM1(nip,nip), inv_dummyM1(nip,nip), dummyM2(gnc,nip);

  IntPointStress.resize(nip,(dim-1)*3);
  //for (unsigned k=0; k<gnc; ++k)
    //   if (NodePtr[k]->counter==0) NodePtr[k]->Stresses.resize((dim-1)*3);


  if ( dim == 3 ) {
    D.resize(6,6);
    copy(MaterialPtr -> Get(Material::ThreeDSolid), D);
  }
  else if ( dim == 2 ) {
    D.resize(3,3);
    copy(MaterialPtr -> Get(Material::PlaneStress12), D);
  }  
  else {
    cerr << "\nERROR in Element::CalcEM()" << endl;
    cerr << "DofSet::ElementDimension returned " << endl;
    cerr << "something else than 2 or 3 " << endl;
    exit(0);
    }

  for (unsigned k = 0 ; k < gnc ; ++k){
    for (unsigned l = 0 ; l < gds ; ++l){
      u[k*gds+l] =  NodeVec[k]->GetDeformations()[l];
    }
  }

  for (unsigned i = 0 ; i < nip; ++i) {

    EvalDerivedShapeFunc(shapefunc, IntPoints, i);
    
    EvalJacobi(J, shapefunc, IntPoints, i);
    stCoord2globalCoord(J, shapefunc);
    SetBMatrix(shapefunc,B,FullIntegration);   

    for (unsigned j=0; j<((dim-1)*3); ++j) {epsilon[j]=0.0;}
    mult(B,u,epsilon);
    mult(D,epsilon,IntPointStress[i]);
  }

  if (nip==8) {
    for (unsigned i = 0 ; i < nip; ++i) {
      dummyM1(i,0) = 1;
      dummyM1(i,1) = IntPoints(i,0);
      dummyM1(i,2) = IntPoints(i,1);
      dummyM1(i,3) = IntPoints(i,2);
      dummyM1(i,4) = IntPoints(i,0) * IntPoints(i,1);
      dummyM1(i,5) = IntPoints(i,1) * IntPoints(i,2);
      dummyM1(i,6) = IntPoints(i,0) * IntPoints(i,2);
      dummyM1(i,7) = IntPoints(i,0) * IntPoints(i,1) * IntPoints(i,2);
    }
  }
  else if (nip==4 && dim==3) {
    for (unsigned i = 0 ; i < nip; ++i) {
      dummyM1(i,0) = 1;
      dummyM1(i,1) = IntPoints(i,0);
      dummyM1(i,2) = IntPoints(i,1);
      dummyM1(i,3) = IntPoints(i,2);
    }
 }
  else if (nip==4) {
    for (unsigned i = 0 ; i < nip; ++i) {
      dummyM1(i,0) = 1;
      dummyM1(i,1) = IntPoints(i,0);
      dummyM1(i,2) = IntPoints(i,1);
      dummyM1(i,3) = IntPoints(i,0) * IntPoints(i,1);
    }
 }
  else if (nip==3) {
    for (unsigned i = 0 ; i < nip; ++i) {
      dummyM1(i,0) = 1;
      dummyM1(i,1) = IntPoints(i,0);
      dummyM1(i,2) = IntPoints(i,1);
    }
 }
  lu_inversion(dummyM1, inv_dummyM1);

  if (nip==8) {
    for (unsigned i = 0 ; i < gnc; ++i) {
      dummyM2(i,0) = 1;
      dummyM2(i,1) = GetlNodeCoords()(i,0);
      dummyM2(i,2) = GetlNodeCoords()(i,1);
      dummyM2(i,3) = GetlNodeCoords()(i,2);
      dummyM2(i,4) = GetlNodeCoords()(i,0) * GetlNodeCoords()(i,1);
      dummyM2(i,5) = GetlNodeCoords()(i,1) * GetlNodeCoords()(i,2);
      dummyM2(i,6) = GetlNodeCoords()(i,0) * GetlNodeCoords()(i,2);
      dummyM2(i,7) = GetlNodeCoords()(i,0) * GetlNodeCoords()(i,1) * GetlNodeCoords()(i,2);
    }
  }
  else if (nip==4 && dim==3) {
    for (unsigned i = 0 ; i < gnc; ++i) {
      dummyM2(i,0) = 1;
      dummyM2(i,1) = GetlNodeCoords()(i,0);
      dummyM2(i,2) = GetlNodeCoords()(i,1);
      dummyM2(i,3) = GetlNodeCoords()(i,2);
    }
  }
  else if (nip==4) {
    for (unsigned i = 0 ; i < gnc; ++i) {
      dummyM2(i,0) = 1;
      dummyM2(i,1) = GetlNodeCoords()(i,0);
      dummyM2(i,2) = GetlNodeCoords()(i,1);
      dummyM2(i,3) = GetlNodeCoords()(i,0) * GetlNodeCoords()(i,1);
    }
  }
  else if (nip==3) {
    for (unsigned i = 0 ; i < gnc; ++i) {
      dummyM2(i,0) = 1;
      dummyM2(i,1) = GetlNodeCoords()(i,0);
      dummyM2(i,2) = GetlNodeCoords()(i,1);
    }
 }

  for (unsigned i=0; i<(dim-1)*3; i++)
    {
      mtl::set(p,0.0);
      mult(inv_dummyM1, columns(IntPointStress)[i], p);
      mult(dummyM2, p, columns(sigma_nodes)[i]);
    }

  for (unsigned k=0; k<gnc; ++k)
    {
      NodeVec[k]->SetStresses(sigma_nodes[k]);
    }

}
