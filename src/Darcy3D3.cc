//-----------------------------------------------------------------------------
// Darcy3D3.cc
//
// begin     : Jun 16 2003
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
   

#include "Darcy3D3.h"

using namespace fe_base;

////
//// Initialize static data members of class "Darcy2D3".
////

const int    Darcy3D3::NodeCount		= 3;
const int    Darcy3D3::Id			= 57;
const string Darcy3D3::Name	        	= "LCM 3D Membrane";
LcmDofSet Darcy3D3::ElementDofSet		= LcmDofSet("100000");


////
//! Implementation of element stiffness matrix computation in class "Darcy3D3".
//// 
Dense_Matrix Darcy3D3::CalcEM()
{

  Dense_Vector B(3,0.0), C(3,0.0);

  B[0] = b[1] - c[1];
  B[1] = c[1] - a[1];
  B[2] = a[1] - b[1];

  C[0] = c[0] - b[0];
  C[1] = a[0] - c[0];
  C[2] = b[0] - a[0]; 

  double D = a[0] * B[0] + b[0] * B[1] + c[0] * B[2];
  Dense_Matrix matB(3,3), matC(3,3), A(3,3);

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

  mtl::set_value(A,0.0);

  if ( MaterialPtr -> ClassName() == "TransverseIsotropicMaterial23" )
    {
      mtl::add(scaled(matB, MaterialPtr->Get("K1")),A);
      mtl::add(scaled(matC, MaterialPtr->Get("K2")),A); }
  else
    {
      mtl::add(scaled(matB, MaterialPtr->Get("K")),A);
      mtl::add(scaled(matC, MaterialPtr->Get("K")),A);
    }
  
  mtl::scale(A,0.5/D);
  mtl::scale(A,-1.0);

  double t;
  if (GetPropertiesPtr() != NULL) {
    t = GetPropertiesPtr()->GetDouble("Thickness");
    mtl::scale(A,t);
  }

  //  cout << "nodes nr " << GetNode(0)->number << ", " << GetNode(1)->number << ", " << GetNode(2)->number << endl;
  //  print_all_matrix(A);

  return A;

}

void Darcy3D3::CompEleVelocity(){

  double D=a[0]*(b[1]-c[1]) + b[0]*(c[1]-a[1]) + c[0]*(a[1]-b[1]);
  vx=0; vy=0;
  
  vx += - GetNodeIter(0)->Pressure / D * (b[1] - c[1]);
  vy += - GetNodeIter(0)->Pressure / D * (c[0] - b[0]);
  
  vx += - GetNodeIter(1)->Pressure / D * (c[1] - a[1]);
  vy += - GetNodeIter(1)->Pressure / D * (a[0] - c[0]);
  
  vx += - GetNodeIter(2)->Pressure / D * (a[1] - b[1]);
  vy += - GetNodeIter(2)->Pressure / D * (b[0] - a[0]);
  
  if ( MaterialPtr -> ClassName() == "TransverseIsotropicMaterial23" ){
    vx *= MaterialPtr->Get("K1")/MaterialPtr->Get("visc");
    vy *= MaterialPtr->Get("K2")/MaterialPtr->Get("visc");
  }
  else {
    vx *= MaterialPtr->Get("K")/MaterialPtr->Get("visc");
    vy *= MaterialPtr->Get("K")/MaterialPtr->Get("visc");
  }
}


void Darcy3D3::CompVolFlux(){

  double t;
  if (GetPropertiesPtr() != NULL)
    t = GetPropertiesPtr()->GetDouble("Thickness");
  else t=1;

  // Flow between Nodes 0 and 1
  unsigned i=0;

  double ny = 0.5*(a[0] + b[0]) - s[0];

  double nx =-0.5*(a[1]+ b[1]) + s[1];

  double A = sqrt(nx*nx+ny*ny);

  nx /= A; ny /= A;
  if ( nx*vx + ny*vy > 0 ) {

    GetNodeIter(i)->q -= (nx*vx + ny*vy) * A * t * GetNodeIter(i)->Fill ;
    GetNodeIter((i+1)%3)->q += (nx*vx + ny*vy) * A * t * GetNodeIter(i)->Fill;
      }
  else{

    GetNodeIter(i)->q -= (nx*vx + ny*vy) * A * t * GetNodeIter((i+1)%3)->Fill ;
    GetNodeIter((i+1)%3)->q += (nx*vx + ny*vy) * A * t * GetNodeIter((i+1)%3)->Fill ;
  }

  // Flow between Nodes 1 and 2
  i=1;
  ny = 0.5*(b[0] + c[0]) - s[0];
  nx =-0.5*(b[1] + c[1]) + s[1];
  A = sqrt(nx*nx+ny*ny);
  nx /= A; ny /= A;

  if ( nx*vx + ny*vy > 0 ) {

    GetNodeIter(i)->q -= (nx*vx + ny*vy) * A * t * GetNodeIter(i)->Fill ;
    GetNodeIter((i+1)%3)->q += (nx*vx + ny*vy) * A * t * GetNodeIter(i)->Fill;
  }
  else{

    GetNodeIter(i)->q -= (nx*vx + ny*vy) * A * t * GetNodeIter((i+1)%3)->Fill ;
    GetNodeIter((i+1)%3)->q += (nx*vx + ny*vy) * A * t * GetNodeIter((i+1)%3)->Fill ;
  }

  // Flow between Nodes 2 and 0
  i=2;
  ny = 0.5*(c[0] + a[0]) - s[0];
  nx =-0.5*(c[1] + a[1]) + s[1];

  A = sqrt(nx*nx+ny*ny);
  nx /= A; ny /= A;
  if ( nx*vx + ny*vy > 0 ) {
    GetNodeIter(i)->q -= (nx*vx + ny*vy) * A * t * GetNodeIter(i)->Fill ;
    GetNodeIter((i+1)%3)->q += (nx*vx + ny*vy) * A * t * GetNodeIter(i)->Fill;

  }

  else{

    GetNodeIter(i)->q -= (nx*vx + ny*vy) * A * t * GetNodeIter((i+1)%3)->Fill ;
    GetNodeIter((i+1)%3)->q += (nx*vx + ny*vy) * A * t * GetNodeIter((i+1)%3)->Fill ;
  }
}


void Darcy3D3::EvalCV() { 
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



 if (!transformed) {

    Dense_Vector n1(3,0.0), n2(3,0.0), n3(3,0.0);
    Dense_Vector ab(3,0.0), bc(3,0.0);
    Dense_Matrix T(3,3), T_inv(3,3);
    
    a.resize(3);
    b.resize(3);
    c.resize(3);

    a = NodeVec[0]->Get();
    b = NodeVec[1]->Get();
    c = NodeVec[2]->Get();
    
    ab[0]=b[0]-a[0];
    ab[1]=b[1]-a[1];
    ab[2]=b[2]-a[2];
    
    bc[0]=c[0]-b[0];
    bc[1]=c[1]-b[1];
    bc[2]=c[2]-b[2];
    
    cross_prod_3d(ab,bc,n1);
    double temp = two_norm(n1);
    scale(n1,1.0/temp);

    if ( EleCoordSysPtr == NULL) {     
      mtl::copy(ab,n2);
    }
   
    else {
      Dense_Vector tmp(3,0.0);
      Dense_Matrix cs(3,3);
      // cs expresses the local-cs vectors in terms of the global cs 
      inversion(EleCoordSysPtr ->Euler2Mat(),cs);
      // cout << "euler2mat, dann cs" << endl;
      //  print_all_matrix(EleCoordSysPtr ->Euler2Mat());
      //   cout << EleCoordSysPtr -> Get("Thxy") << ", " << EleCoordSysPtr -> Get("Thzx") << ", " << EleCoordSysPtr -> Get("Thyz") << endl;
      //  print_all_matrix(cs);
      double dotpr = cs(0,0)*n1[0] + cs(1,0)*n1[1] + cs(2,0)*n1[2]; //dot product, angle btw local x-axis and n1
      //  cout << "dotpr = " << dotpr;
      if (abs(dotpr) < 0.7072) { // normal case, local x-axis is projected onto element
	// angle of 45 degree is used as in ansys
	add(scaled(n1,-dotpr),tmp);
	//	printv(tmp);
	add(columns(cs)[0],tmp);
	//	printv(tmp);
	add(tmp,n2); }
      else { // shit case, local x-axis (almost) normal to element -> take y-axis instead
	dotpr = cs(0,1)*n1[0] + cs(1,1)*n1[1] + cs(2,1)*n1[2]; //dot product, angle btw global y-axis and n1
	add(scaled(n1,-dotpr),tmp);
	add(columns(cs)[1],tmp);
	add(tmp,n2); }
    }

    temp = two_norm(n2);
    scale(n2,1.0/temp);
    cross_prod_3d(n1,n2,n3);
    temp = two_norm(n3);
    scale(n3,1.0/temp);
  
    mtl::copy(n2,columns(T)[0]);
    mtl::copy(n3,columns(T)[1]);
    mtl::copy(n1,columns(T)[2]);
  
    inversion(T,T_inv);
    
    mtl::set(a,0.0);
    mtl::set(b,0.0);
    mtl::set(c,0.0);
    mult(T_inv,ab,b);
    mult(T_inv,bc,c);
    add(b,c);
    
    s.resize(3);
    s[0] = (a[0]+b[0]+c[0])/3.0;
    s[1] = (a[1]+b[1]+c[1])/3.0;
    s[2] = (a[2]+b[2]+c[2])/3.0;
    
    transformed=true;
    //cout << "a = " << a << ", b = " << b << ", c = " << c << endl;
    //  cout << "a, b, c " << endl;
    //printv(a);
    //printv(b);
    //printv(c);
  }

}
