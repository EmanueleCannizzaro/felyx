//-----------------------------------------------------------------------------
// Darcy2D3.cc
//
// begin     : Jan 10 2003
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
   

#include "Darcy2D3.h"

using namespace fe_base;

////
//// Initialize static data members of class "Darcy2D3".
////

const int    Darcy2D3::NodeCount		= 3;
const int    Darcy2D3::Id			= 55;
const string Darcy2D3::Name	        	= "LCM 2D Membrane";
LcmDofSet Darcy2D3::ElementDofSet		= LcmDofSet("100000");


////
//! Implementation of element stiffness matrix computation in class "Darcy2D3".
//// 
Dense_Matrix Darcy2D3::CalcEM()
{
  if (!emExists) {
    emExists=true;
    Dense_Vector B(3,0.0), C(3,0.0);
    
    B[0] = NodeVec[1]->Cy - NodeVec[2]->Cy;
    B[1] = NodeVec[2]->Cy - NodeVec[0]->Cy;
    B[2] = NodeVec[0]->Cy - NodeVec[1]->Cy;
    
    C[0] = NodeVec[2]->Cx - NodeVec[1]->Cx;
    C[1] = NodeVec[0]->Cx - NodeVec[2]->Cx;
    C[2] = NodeVec[1]->Cx - NodeVec[0]->Cx;
    
    double D = NodeVec[0]->Cx * B[0] + NodeVec[1]->Cx * B[1] + NodeVec[2]->Cx * B[2];
    
    Dense_Matrix matB(3,3), matC(3,3), matD(3,3);
    A.resize(3,3);
    
    matB(0,0) = B[0] * B[0];
    matB(0,1) = B[0] * B[1];
    matB(0,2) = B[0] * B[2];
    matB(1,0) = B[1] * B[0];
    matB(1,1) = B[1] * B[1];
    matB(1,2) = B[1] * B[2];
    matB(2,0) = B[2] * B[0];
    matB(2,1) = B[2] * B[1];
    matB(2,2) = B[2] * B[2];
    
    matC(0,0) = C[0] * C[0];
    matC(0,1) = C[0] * C[1];
    matC(0,2) = C[0] * C[2];
    matC(1,0) = C[1] * C[0];
    matC(1,1) = C[1] * C[1];
    matC(1,2) = C[1] * C[2];
    matC(2,0) = C[2] * C[0];
    matC(2,1) = C[2] * C[1];
    matC(2,2) = C[2] * C[2];
    
    matD(0,0) = B[0] * C[0] + C[0] * B[0];
    matD(0,1) = B[0] * C[1] + C[0] * B[1];
    matD(0,2) = B[0] * C[2] + C[0] * B[2];
    matD(1,0) = B[1] * C[0] + C[1] * B[0];
    matD(1,1) = B[1] * C[1] + C[1] * B[1];
    matD(1,2) = B[1] * C[2] + C[1] * B[2];
    matD(2,0) = B[2] * C[0] + C[2] * B[0];
    matD(2,1) = B[2] * C[1] + C[2] * B[1];
    matD(2,2) = B[2] * C[2] + C[2] * B[2];
    
    mtl::set_value(A,0.0);
    if ( MaterialPtr -> ClassName() == "TransverseIsotropicMaterial23" )
      {
	if ( EleCoordSysPtr == NULL){
	  mtl::add(scaled(matB, MaterialPtr->Get("K1")),A);
	  mtl::add(scaled(matC, MaterialPtr->Get("K2")),A);
	}
	else {
	  double theta = EleCoordSysPtr -> Get("Thxy");
	  theta *= 2.0*3.14159265358/360.0;
	  double K1 = MaterialPtr->Get("K1");
	  double K2 = MaterialPtr->Get("K2");
	  double Kxx = K1*cos(theta)*cos(theta)+K2*sin(theta)*sin(theta);
	  double Kyy = K1*sin(theta)*sin(theta)+K2*cos(theta)*cos(theta);
	  double Kxy = -(K2-K1)*sin(theta)*cos(theta);
	  mtl::add(scaled(matB, Kxx),A);
	  mtl::add(scaled(matC, Kyy),A);
	  mtl::add(scaled(matD, Kxy),A);
	} 
      }
    else
      {
	mtl::add(scaled(matB, MaterialPtr->Get("K")),A);
	mtl::add(scaled(matC, MaterialPtr->Get("K")),A);
      }
    
    mtl::scale(A,0.5/D);
   // mtl::scale(A,-1.0);
  }
  return A;
  
}


void Darcy2D3::CompEleVelocity(){
  
  double D=GetNodeIter(0)->Cx*(GetNodeIter(1)->Cy-GetNodeIter(2)->Cy)
    +GetNodeIter(1)->Cx*(GetNodeIter(2)->Cy-GetNodeIter(0)->Cy)
    +GetNodeIter(2)->Cx*(GetNodeIter(0)->Cy-GetNodeIter(1)->Cy);
  vx=0; vy=0;
  for (unsigned i=0; i<GetNodeCount(); ++i) {
    vx += - GetNodeIter(i)->Pressure / D * (GetNodeIter((i+1)%3)->Cy - GetNodeIter((i+2)%3)->Cy);
    vy += - GetNodeIter(i)->Pressure / D * (GetNodeIter((i+2)%3)->Cx - GetNodeIter((i+1)%3)->Cx);
  }
  if ( MaterialPtr -> ClassName() == "TransverseIsotropicMaterial23" ){
    if ( EleCoordSysPtr == NULL){
      vx *= MaterialPtr->Get("K1")/MaterialPtr->Get("visc");
      vy *= MaterialPtr->Get("K2")/MaterialPtr->Get("visc");
    }
    else {
      double theta = EleCoordSysPtr -> Get("Thxy");
      theta *= 2.0*3.14159265358/360.0;
      double K1 = MaterialPtr->Get("K1");
      double K2 = MaterialPtr->Get("K2");
      double Kxx = K1*cos(theta)*cos(theta)+K2*sin(theta)*sin(theta);
      double Kyy = K1*sin(theta)*sin(theta)+K2*cos(theta)*cos(theta);
      double Kxy = -(K2-K1)*sin(theta)*cos(theta);
      double visc = MaterialPtr->Get("visc");
      double B = vx;
      double C = vy;
      vx = Kxx/visc*B + Kxy/visc*C;
      vy = Kxy/visc*B + Kyy/visc*C;
    }
  }
  else {
    vx *= MaterialPtr->Get("K")/MaterialPtr->Get("visc");
    vy *= MaterialPtr->Get("K")/MaterialPtr->Get("visc");
  }
}


void Darcy2D3::CompVolFlux() {
  for (unsigned i=0; i<3; ++i){
    double ny = 0.5*(GetNodeIter(i)->Cx + GetNodeIter((i+1)%3)->Cx) - sx;
    double nx =-0.5*(GetNodeIter(i)->Cy + GetNodeIter((i+1)%3)->Cy) + sy;
    double A = sqrt(nx*nx+ny*ny);
    nx /= A; ny /= A;
    if ( nx*vx + ny*vy > 0 ) {
      GetNodeIter(i)->q -= (nx*vx + ny*vy) * A * GetNodeIter(i)->Fill ;
      GetNodeIter((i+1)%3)->q += (nx*vx + ny*vy) * A * GetNodeIter(i)->Fill;
    }
    else{
      GetNodeIter(i)->q -= (nx*vx + ny*vy) * A * GetNodeIter((i+1)%3)->Fill ;
      GetNodeIter((i+1)%3)->q += (nx*vx + ny*vy) * A * GetNodeIter((i+1)%3)->Fill ;
    }    
  }
}


void Darcy2D3::EvalCV() { 
  Dense_Vector ab(3,0.0), ac(3,0.0), ad(3,0.0), temp(3,0.0);
  ab[0] = GetNodeIter(1)->Cx - GetNodeIter(0)->Cx;
  ab[1] = GetNodeIter(1)->Cy - GetNodeIter(0)->Cy;
  ab[2] = GetNodeIter(1)->Cz - GetNodeIter(0)->Cz;
  ac[0] = GetNodeIter(2)->Cx - GetNodeIter(0)->Cx;
  ac[1] = GetNodeIter(2)->Cy - GetNodeIter(0)->Cy;
  ac[2] = GetNodeIter(2)->Cz - GetNodeIter(0)->Cz;
  cross_prod_3d(ab,ac,temp);
  if (GetPropertiesPtr() != NULL) {
    GetNodeIter(0)->Volume += (two_norm(temp)/6.0) * GetPropertiesPtr()->GetDouble("Thickness");
    GetNodeIter(1)->Volume += (two_norm(temp)/6.0) * GetPropertiesPtr()->GetDouble("Thickness");
    GetNodeIter(2)->Volume += (two_norm(temp)/6.0) * GetPropertiesPtr()->GetDouble("Thickness"); }
  else {
    GetNodeIter(0)->Volume += (two_norm(temp)/6.0);
    GetNodeIter(1)->Volume += (two_norm(temp)/6.0);
    GetNodeIter(2)->Volume += (two_norm(temp)/6.0); }
}
