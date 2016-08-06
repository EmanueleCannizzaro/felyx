//-----------------------------------------------------------------------------
// Darcy3D8.cc
//
// begin     : Mar 2 2004
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
#include "Darcy3D8.h"

using namespace fe_base;

////
//// Initialize static data members of class "Darcy3D8"
////
const int    Darcy3D8::NodeCount     	= 8;
const int    Darcy3D8::Id	       	= 70;
const string Darcy3D8::Name	       	= "LCM 3D Brick";
LcmDofSet    Darcy3D8::ElementDofSet 	= LcmDofSet("100000");
Dense_Matrix Darcy3D8::IP               = GetIntegrationPoints();

// EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to the standard
// coords L1, L2, L3 and L4 in a specific point.
// Its rows correspond to a specific shape function. The cols to the direction of derivation 
// EvalPoints is a Matrix of points, where the shape functions are to be evaluated. IP determines the specific point,
// means row of Evalpoints. The EvalPoints are given in volumetric coords
//----------------------------------------------------------------------------------------------------------------------

Dense_Matrix Darcy3D8::GetIntegrationPoints() {
  double coord = 1.0 / sqrt(3.0);

  Dense_Matrix IP(8,3); 
  
  IP(0,0) = -coord; IP(0,1) = -coord; IP(0,2) = -coord;
  IP(1,0) =  coord; IP(1,1) = -coord; IP(1,2) = -coord;
  IP(2,0) =  coord; IP(2,1) =  coord; IP(2,2) = -coord;
  IP(3,0) = -coord; IP(3,1) =  coord; IP(3,2) = -coord;
  IP(4,0) = -coord; IP(4,1) = -coord; IP(4,2) =  coord;
  IP(5,0) =  coord; IP(5,1) = -coord; IP(5,2) =  coord;
  IP(6,0) =  coord; IP(6,1) =  coord; IP(6,2) =  coord;
  IP(7,0) = -coord; IP(7,1) =  coord; IP(7,2) =  coord;
  return IP;
}

void Darcy3D8::EvalDerivedShapeFunc(Dense_Matrix& shapefunc, Dense_Vector& v)
{
  shapefunc(0,0) = -0.125 * (1 - v[1]) * (1 - v[2]);
  shapefunc(1,0) =  0.125 * (1 - v[1]) * (1 - v[2]);
  shapefunc(2,0) =  0.125 * (1 + v[1]) * (1 - v[2]);
  shapefunc(3,0) = -0.125 * (1 + v[1]) * (1 - v[2]);
  shapefunc(4,0) = -0.125 * (1 - v[1]) * (1 + v[2]);
  shapefunc(5,0) =  0.125 * (1 - v[1]) * (1 + v[2]);
  shapefunc(6,0) =  0.125 * (1 + v[1]) * (1 + v[2]);
  shapefunc(7,0) = -0.125 * (1 + v[1]) * (1 + v[2]);
  
  shapefunc(0,1) = -0.125 * (1 - v[0]) * (1 - v[2]);
  shapefunc(1,1) = -0.125 * (1 + v[0]) * (1 - v[2]);
  shapefunc(2,1) =  0.125 * (1 + v[0]) * (1 - v[2]);
  shapefunc(3,1) =  0.125 * (1 - v[0]) * (1 - v[2]);
  shapefunc(4,1) = -0.125 * (1 - v[0]) * (1 + v[2]);
  shapefunc(5,1) = -0.125 * (1 + v[0]) * (1 + v[2]);
  shapefunc(6,1) =  0.125 * (1 + v[0]) * (1 + v[2]);
  shapefunc(7,1) =  0.125 * (1 - v[0]) * (1 + v[2]);
  
  shapefunc(0,2) = -0.125 * (1 - v[0]) * (1 - v[1]);
  shapefunc(1,2) = -0.125 * (1 + v[0]) * (1 - v[1]);
  shapefunc(2,2) = -0.125 * (1 + v[0]) * (1 + v[1]);
  shapefunc(3,2) = -0.125 * (1 - v[0]) * (1 + v[1]);
  shapefunc(4,2) =  0.125 * (1 - v[0]) * (1 - v[1]);
  shapefunc(5,2) =  0.125 * (1 + v[0]) * (1 - v[1]);
  shapefunc(6,2) =  0.125 * (1 + v[0]) * (1 + v[1]);
  shapefunc(7,2) =  0.125 * (1 - v[0]) * (1 + v[1]);  
}

void Darcy3D8::Local2Global(Dense_Matrix& v)
{
  unsigned col=v.ncols();
  Dense_Matrix temp(3,col);
  for (unsigned j=0; j<col; ++j) {
    for (unsigned i=0; i<3; ++i) {
      temp(i,j) += 0.125*(1.0-v(0,j))*(1.0-v(1,j))*(1.0-v(2,j))*GetNodeIter(0)->GetCoords()[i];
      temp(i,j) += 0.125*(1.0+v(0,j))*(1.0-v(1,j))*(1.0-v(2,j))*GetNodeIter(1)->GetCoords()[i];
      temp(i,j) += 0.125*(1.0+v(0,j))*(1.0+v(1,j))*(1.0-v(2,j))*GetNodeIter(2)->GetCoords()[i];
      temp(i,j) += 0.125*(1.0-v(0,j))*(1.0+v(1,j))*(1.0-v(2,j))*GetNodeIter(3)->GetCoords()[i];
      temp(i,j) += 0.125*(1.0-v(0,j))*(1.0-v(1,j))*(1.0+v(2,j))*GetNodeIter(4)->GetCoords()[i];
      temp(i,j) += 0.125*(1.0+v(0,j))*(1.0-v(1,j))*(1.0+v(2,j))*GetNodeIter(5)->GetCoords()[i];
      temp(i,j) += 0.125*(1.0+v(0,j))*(1.0+v(1,j))*(1.0+v(2,j))*GetNodeIter(6)->GetCoords()[i];
      temp(i,j) += 0.125*(1.0-v(0,j))*(1.0+v(1,j))*(1.0+v(2,j))*GetNodeIter(7)->GetCoords()[i];
    }
  }
  mtl::copy(temp,v);
}

Dense_Matrix Darcy3D8::CalcEM() {

  double detJ;
  unsigned nSize = GetEMSize();
  unsigned dim = 3;

  //Initialization of some Matrices 
  Dense_Matrix 
    J(dim,dim),                     /* Jacobimatrix */
    K(nSize,nSize),                 /* "Stiffness" Matrix */
    shapefunc(GetNodeCount(), 3),   /* derived - shapefuncs */
    D,                              /* Permeabilty Matrix */
    dummy(GetNodeCount(),3);

  //Get the material data
  D.resize(3,3);
  
  for (unsigned i = 0 ; i < 8; i++)
    {
      mtl::set_value(D,0.0);
      if ( MaterialPtr -> ClassName() == "TransverseIsotropicMaterial23" ) {
	mtl::copy(K_Mat,D);
      }
      else {
	D(0,0) = MaterialPtr->Get("K");
	D(1,1) = MaterialPtr->Get("K");
	D(2,2) = MaterialPtr->Get("K");
      }
      
      Dense_Vector v(3);
      v[0]=IP(i,0); v[1]=IP(i,1); v[2]=IP(i,2);
      EvalDerivedShapeFunc(shapefunc, v);
      
      //Compute the Jacobian matrix
      EvalJacobi(J, shapefunc);
      
      //Transforme the shapefunctions to global x, y, z coordinates
      stCoord2globalCoord(J, shapefunc);
      
      detJ = determinant(J);
      
      scale(D,detJ);
      
      set_value(dummy,0.0);
      mult(shapefunc,D,dummy);
      mult(dummy,trans(shapefunc),K);
      
    }

  return K;
}


void Darcy3D8::EvalJacobi(Dense_Matrix& J_, Dense_Matrix& shapefunc_){

  //initialize the matrix to zero   
  mtl::set_value(J_,0.0);

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



void Darcy3D8::stCoord2globalCoord(const Dense_Matrix& Jacobi_,  Dense_Matrix& shapefunc_){

  //inversion of the Jacobian
  unsigned size = Jacobi_.nrows();
  Dense_Matrix invJacobi(size, size);

  inversion(Jacobi_, invJacobi);

    Dense_Vector dummyvec1(size), dummyvec2(size);
    for (unsigned k = 0 ; k < GetNodeCount() ; k++)
      {
	dummyvec1[0] = shapefunc_(k,0);
	dummyvec1[1] = shapefunc_(k,1);
	dummyvec1[2] = shapefunc_(k,2);
	
	mtl::set_value(dummyvec2,0.0);
	mult(invJacobi, dummyvec1, dummyvec2);
	
	shapefunc_(k,0) = dummyvec2[0];
	shapefunc_(k,1) = dummyvec2[1];
	shapefunc_(k,2) = dummyvec2[2];
      }

}

// compute pressure gradient dp (global coord.) at postition v (local coord.)
void Darcy3D8::CompPressGrad(Dense_Vector& v, Dense_Vector& dp) {

  Dense_Matrix shapefunc(GetNodeCount(),3), dummy(3,GetNodeCount());
  Dense_Matrix J(3,3), invJ(3,3);
  set_value(dp,0.0);
  EvalDerivedShapeFunc(shapefunc, v);
  EvalJacobi(J, shapefunc);
  inversion(J, invJ);
  mult(invJ,trans(shapefunc),dummy);
  for (unsigned j=0; j<3; ++j) {
    for (unsigned i=0; i<8; ++i) {
      dp[j]+=GetNodeIter(i)->Pressure*dummy(j,i);
    }
  }
}

void Darcy3D8::CompVolFlux(){
  // matrix storing local node coordinates of scv

  Dense_Matrix C(3,8);
  // scv 1
  C(0,0)=-1; C(1,0)=-1; C(2,0)=-1; C(0,1)= 0; C(1,1)=-1; C(2,1)=-1; 
  C(0,2)= 0; C(1,2)= 0; C(2,2)=-1; C(0,3)=-1; C(1,3)= 0; C(2,3)=-1; 
  C(0,4)=-1; C(1,4)=-1; C(2,4)= 0; C(0,5)= 0; C(1,5)=-1; C(2,5)= 0; 
  C(0,6)= 0; C(1,6)= 0; C(2,6)= 0; C(0,7)=-1; C(1,7)= 0; C(2,7)= 0; 
  unsigned FluxBCFG = 1, FluxEFGH = 4, FluxCDGH = 3, myself = 0;

  Flux(C,FluxBCFG,FluxEFGH,FluxCDGH,myself);
  // scv 2
  C(0,0)= 1; C(1,0)=-1; C(2,0)=-1; C(0,1)= 1; C(1,1)= 0; C(2,1)=-1; 
  C(0,2)= 0; C(1,2)= 0; C(2,2)=-1; C(0,3)= 0; C(1,3)=-1; C(2,3)=-1; 
  C(0,4)= 1; C(1,4)=-1; C(2,4)= 0; C(0,5)= 1; C(1,5)= 0; C(2,5)= 0; 
  C(0,6)= 0; C(1,6)= 0; C(2,6)= 0; C(0,7)= 0; C(1,7)=-1; C(2,7)= 0; 
  FluxBCFG = 2, FluxEFGH = 5, FluxCDGH = 0, myself = 1;
  Flux(C,FluxBCFG,FluxEFGH,FluxCDGH,myself);
  // scv 3
  C(0,0)= 1; C(1,0)= 1; C(2,0)=-1; C(0,1)= 0; C(1,1)= 1; C(2,1)=-1; 
  C(0,2)= 0; C(1,2)= 0; C(2,2)=-1; C(0,3)= 1; C(1,3)= 0; C(2,3)=-1; 
  C(0,4)= 1; C(1,4)= 1; C(2,4)= 0; C(0,5)= 0; C(1,5)= 1; C(2,5)= 0; 
  C(0,6)= 0; C(1,6)= 0; C(2,6)= 0; C(0,7)= 1; C(1,7)= 0; C(2,7)= 0; 
  FluxBCFG = 3, FluxEFGH = 6, FluxCDGH = 1, myself = 2;
  Flux(C,FluxBCFG,FluxEFGH,FluxCDGH,myself);
  // scv 4
  C(0,0)=-1; C(1,0)= 1; C(2,0)=-1; C(0,1)=-1; C(1,1)= 0; C(2,1)=-1; 
  C(0,2)= 0; C(1,2)= 0; C(2,2)=-1; C(0,3)= 0; C(1,3)= 1; C(2,3)=-1; 
  C(0,4)=-1; C(1,4)= 1; C(2,4)= 0; C(0,5)=-1; C(1,5)= 0; C(2,5)= 0; 
  C(0,6)= 0; C(1,6)= 0; C(2,6)= 0; C(0,7)= 0; C(1,7)= 1; C(2,7)= 0; 
  FluxBCFG = 0, FluxEFGH = 7, FluxCDGH = 2, myself = 3;
  Flux(C,FluxBCFG,FluxEFGH,FluxCDGH,myself);
  // scv 5
  C(0,0)=-1; C(1,0)=-1; C(2,0)= 1; C(0,1)=-1; C(1,1)= 0; C(2,1)= 1; 
  C(0,2)= 0; C(1,2)= 0; C(2,2)= 1; C(0,3)= 0; C(1,3)=-1; C(2,3)= 1; 
  C(0,4)=-1; C(1,4)=-1; C(2,4)= 0; C(0,5)=-1; C(1,5)= 0; C(2,5)= 0; 
  C(0,6)= 0; C(1,6)= 0; C(2,6)= 0; C(0,7)= 0; C(1,7)=-1; C(2,7)= 0; 
  FluxBCFG = 7, FluxEFGH = 0, FluxCDGH = 5, myself = 4;
  Flux(C,FluxBCFG,FluxEFGH,FluxCDGH,myself);
  // scv 6
   C(0,0)= 1; C(1,0)=-1; C(2,0)= 1; C(0,1)= 0; C(1,1)=-1; C(2,1)= 1; 
  C(0,2)= 0; C(1,2)= 0; C(2,2)= 1; C(0,3)= 1; C(1,3)= 0; C(2,3)= 1; 
  C(0,4)= 1; C(1,4)=-1; C(2,4)= 0; C(0,5)= 0; C(1,5)=-1; C(2,5)= 0; 
  C(0,6)= 0; C(1,6)= 0; C(2,6)= 0; C(0,7)= 1; C(1,7)= 0; C(2,7)= 0; 
  FluxBCFG = 4, FluxEFGH = 1, FluxCDGH = 6, myself = 5;
  Flux(C,FluxBCFG,FluxEFGH,FluxCDGH,myself);
  // scv 7
  C(0,0)= 1; C(1,0)= 1; C(2,0)= 1; C(0,1)= 1; C(1,1)= 0; C(2,1)= 1; 
  C(0,2)= 0; C(1,2)= 0; C(2,2)= 1; C(0,3)= 0; C(1,3)= 1; C(2,3)= 1; 
  C(0,4)= 1; C(1,4)= 1; C(2,4)= 0; C(0,5)= 1; C(1,5)= 0; C(2,5)= 0; 
  C(0,6)= 0; C(1,6)= 0; C(2,6)= 0; C(0,7)= 0; C(1,7)= 1; C(2,7)= 0; 
  FluxBCFG = 5, FluxEFGH = 2, FluxCDGH = 7, myself = 6;
  Flux(C,FluxBCFG,FluxEFGH,FluxCDGH,myself);
  // scv 8
  C(0,0)=-1; C(1,0)= 1; C(2,0)= 1; C(0,1)= 0; C(1,1)= 1; C(2,1)= 1; 
  C(0,2)= 0; C(1,2)= 0; C(2,2)= 1; C(0,3)=-1; C(1,3)= 0; C(2,3)= 1; 
  C(0,4)=-1; C(1,4)= 1; C(2,4)= 0; C(0,5)= 0; C(1,5)= 1; C(2,5)= 0; 
  C(0,6)= 0; C(1,6)= 0; C(2,6)= 0; C(0,7)=-1; C(1,7)= 0; C(2,7)= 0; 
  FluxBCFG = 6, FluxEFGH = 3, FluxCDGH = 4, myself = 7;
  Flux(C,FluxBCFG,FluxEFGH,FluxCDGH,myself);

}


void Darcy3D8::Flux(Dense_Matrix& C, unsigned BCFG, unsigned EFGH, unsigned CDGH, unsigned myself) {
  // get 3 pressure gradient vectors dp
  Dense_Vector dp1(3), dp2(3), dp3(3);
  Dense_Vector c1(3), c2(3), c3(3);
  c1[0]=(C(0,1)+C(0,2)+C(0,5)+C(0,6))/4.0;
  c1[1]=(C(1,1)+C(1,2)+C(1,5)+C(1,6))/4.0;
  c1[2]=(C(2,1)+C(2,2)+C(2,5)+C(2,6))/4.0;
  c2[0]=(C(0,4)+C(0,5)+C(0,6)+C(0,7))/4.0;
  c2[1]=(C(1,4)+C(1,5)+C(1,6)+C(1,7))/4.0;
  c2[2]=(C(2,4)+C(2,5)+C(2,6)+C(2,7))/4.0;
  c3[0]=(C(0,2)+C(0,3)+C(0,6)+C(0,7))/4.0;
  c3[1]=(C(1,2)+C(1,3)+C(1,6)+C(1,7))/4.0;
  c3[2]=(C(2,2)+C(2,3)+C(2,6)+C(2,7))/4.0;
  CompPressGrad(c1,dp1);
  CompPressGrad(c2,dp2);
  CompPressGrad(c3,dp3);

  // evaluating normals
  Local2Global(C);
  Dense_Vector n1(3), n2(3), n3(3);
  n1[0] = (-(C(1,5)-C(1,1))*(C(2,2)-C(2,1))+(C(2,5)-C(2,1))*(C(1,2)-C(1,1)) +
	   (C(1,5)-C(1,6))*(C(2,2)-C(2,6))-(C(2,5)-C(2,6))*(C(1,2)-C(1,6)))/2.0;
  n1[1] =-(-(C(0,5)-C(0,1))*(C(2,2)-C(2,1))+(C(2,5)-C(2,1))*(C(0,2)-C(0,1)) +
	   (C(0,5)-C(0,6))*(C(2,2)-C(2,6))-(C(2,5)-C(2,6))*(C(0,2)-C(0,6)))/2.0;
  n1[2] = (-(C(0,5)-C(0,1))*(C(1,2)-C(1,1))+(C(1,5)-C(1,1))*(C(0,2)-C(0,1)) +
	   (C(0,5)-C(0,6))*(C(1,2)-C(1,6))-(C(1,5)-C(1,6))*(C(0,2)-C(0,6)))/2.0;

  n2[0] = (-(C(1,4)-C(1,5))*(C(2,6)-C(2,5))+(C(2,4)-C(2,5))*(C(1,6)-C(1,5)) +
	   (C(1,4)-C(1,7))*(C(2,6)-C(2,7))-(C(2,4)-C(2,7))*(C(1,6)-C(1,7)))/2.0;
  n2[1] =-(-(C(0,4)-C(0,5))*(C(2,6)-C(2,5))+(C(2,4)-C(2,5))*(C(0,6)-C(0,5)) +
	   (C(0,4)-C(0,7))*(C(2,6)-C(2,7))-(C(2,4)-C(2,7))*(C(0,6)-C(0,7)))/2.0;
  n2[2] = (-(C(0,4)-C(0,5))*(C(1,6)-C(1,5))+(C(1,4)-C(1,5))*(C(0,6)-C(0,5)) +
	   (C(0,4)-C(0,7))*(C(1,6)-C(1,7))-(C(1,4)-C(1,7))*(C(0,6)-C(0,7)))/2.0;

  n3[0] = (-(C(1,6)-C(1,2))*(C(2,3)-C(2,2))+(C(2,6)-C(2,2))*(C(1,3)-C(1,2)) +
	   (C(1,6)-C(1,7))*(C(2,3)-C(2,7))-(C(2,6)-C(2,7))*(C(1,3)-C(1,7)))/2.0;
  n3[1] =-(-(C(0,6)-C(0,2))*(C(2,3)-C(2,2))+(C(2,6)-C(2,2))*(C(0,3)-C(0,2)) +
	   (C(0,6)-C(0,7))*(C(2,3)-C(2,7))-(C(2,6)-C(2,7))*(C(0,3)-C(0,7)))/2.0;
  n3[2] = (-(C(0,6)-C(0,2))*(C(1,3)-C(1,2))+(C(1,6)-C(1,2))*(C(0,3)-C(0,2)) +
	   (C(0,6)-C(0,7))*(C(1,3)-C(1,7))-(C(1,6)-C(1,7))*(C(0,3)-C(0,7)))/2.0;


  // compute velocities from dp
  Dense_Vector v1(3), v2(3), v3(3);
  //  print_all_matrix(K_Mat);
  mult(K_Mat,dp1,v1);
  scale(v1,-1.0/MaterialPtr->Get("visc"));
  mult(K_Mat,dp2,v2);
  scale(v2,-1.0/MaterialPtr->Get("visc"));
  mult(K_Mat,dp3,v3);
  scale(v3,-1.0/MaterialPtr->Get("visc"));

  // compute flux to neighbor cv's and apply 'em
  double d;
  d=scalar_prod_3d(v1,n1);
  if (d>0) { GetNodeIter(BCFG)->q += d*GetNodeIter(myself)->Fill;
  GetNodeIter(myself)->q -= d*GetNodeIter(myself)->Fill; }

  d=scalar_prod_3d(v2,n2);
  if (d>0) { GetNodeIter(EFGH)->q += d*GetNodeIter(myself)->Fill;
  GetNodeIter(myself)->q -= d*GetNodeIter(myself)->Fill; }

  d=scalar_prod_3d(v3,n3);
  if (d>0) { GetNodeIter(CDGH)->q += d*GetNodeIter(myself)->Fill;
  GetNodeIter(myself)->q -= d*GetNodeIter(myself)->Fill; }
}

void Darcy3D8::CompEleVelocity(){

}

void Darcy3D8::EvalCentroid(){
  
  K_Mat.resize(3,3); 
  
  if ( MaterialPtr -> ClassName() == "TransverseIsotropicMaterial23" ) {
    
    Dense_Matrix D(3,3);
    D(0,0) = MaterialPtr->Get("K1");
    D(1,1) = MaterialPtr->Get("K2");
    D(2,2) = MaterialPtr->Get("K3");
    if (EleCoordSysPtr==NULL){
      mtl::copy(D,K_Mat);}
    else{
      mult_AT_B_A_add(EleCoordSysPtr->GetRotMat(),D,K_Mat);}
  }
  else {
    K_Mat(0,0) = MaterialPtr->Get("K");
    K_Mat(1,1) = MaterialPtr->Get("K");
    K_Mat(2,2) = MaterialPtr->Get("K");
  }

}


void Darcy3D8::EvalCV() {
  
  Dense_Matrix J(3,3), shapefunc(8,3);
  Dense_Vector v(3);
  v[0]=1; v[1]=1; v[2]=1;
  EvalDerivedShapeFunc(shapefunc, v);
      
  EvalJacobi(J, shapefunc);
  double detJ = determinant(J);

  for (unsigned i=0; i<GetNodeCount(); ++i) {
    GetNodeIter(i)->Volume += detJ;
  }
  
}
