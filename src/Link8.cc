//-----------------------------------------------------------------------------
// Link8.cc
//
// begin     : Dec 20 2001
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
   

#include <cmath>
#include "Link8.h"

using namespace fe_base;

////
//// Initialize static data members of class "Link8"
////
const int    Link8::NodeCount		= 2;
const int    Link8::Id	    		= 8;
const string Link8::Name	    	= "Structural 3D Link";
StructDofSet Link8::ElementDofSet	= StructDofSet("111000");

////
//! Implementation of element stiffness matrix computation in class "Link8".
////
Dense_Matrix Link8::CalcEM(){

  // length of link
  double L   = sqrt( pow( NodeVec[1]->Cx - NodeVec[0]->Cx ,2) +
		     pow( NodeVec[1]->Cy - NodeVec[0]->Cy ,2) +
		     pow( NodeVec[1]->Cz - NodeVec[0]->Cz ,2) ) ;

  // Eval direction cosines in order to transform stiffness matrix
  // from local coords (link direction) to global coordinates
  double cx = ( NodeVec[1]->Cx - NodeVec[0]->Cx ) / L;
  double cy = ( NodeVec[1]->Cy - NodeVec[0]->Cy ) / L;
  double cz = ( NodeVec[1]->Cz - NodeVec[0]->Cz ) / L;

  // Transforamtion matrix
  Dense_Matrix T(2,6);
  T(0,0) = cx; T(0,1) = cy; T(0,2) = cz;
  T(1,3) = cx; T(1,4) = cy; T(1,5) = cz;

  // Ensure that Material and Properies ptr are set
  if (MaterialPtr == NULL || PropertiesPtr == NULL){
    cerr << endl << " ERROR: Link8::CalcEM(): " ;
    cerr << "MaterialPtr or PropertiesPtr are set to NULL"  << endl;
    exit(1);
  }

  // Stiffness matrix
  Dense_Matrix K(2,2);
  K(0,0) = 1; 		K(0,1) = -1;
  K(1,0) =-1;		K(1,1) =  1;

  // Multiply K with axial stiffness
  // K *= MaterialPtr->Get("E") * PropertiesPtr->GetDouble("Area") / L;
  scale(K, (MaterialPtr->Get("E") * PropertiesPtr->GetDouble("Area") /L));

  // Evaluate Transformation K = T^T * K * T
  
  //  K = K*T;
  // T.transpose();
  // K = T*K;
  
  Dense_Matrix dummyMatrix(6,6);
  mult_AT_B_A(T,K,dummyMatrix);
  
  // copy(dummyMatrix, K);  
  
  //return K;
  return dummyMatrix;
 
}

////
//! Implementation of element mass matrix computation in class "Link8".
////
Dense_Matrix Link8::CalcEmm(){

  // length of link
  double L   = sqrt( pow( NodeVec[1]->Cx - NodeVec[0]->Cx ,2) +
		     pow( NodeVec[1]->Cy - NodeVec[0]->Cy ,2) +
		     pow( NodeVec[1]->Cz - NodeVec[0]->Cz ,2) ) ;

  // Eval direction cosines in order to transform stiffness matrix
  // from local coords (link direction) to global coordinates
  double cxx = ( NodeVec[1]->Cx - NodeVec[0]->Cx ) / L;
  double cxy = ( NodeVec[1]->Cy - NodeVec[0]->Cy ) / L;
  double cxz = ( NodeVec[1]->Cz - NodeVec[0]->Cz ) / L;

  double L_ = sqrt(cxx*cxx+cxy*cxy);

  double cyx, cyy, cyz, czx, czy, czz;
  if(L_!=0.0)
    {
      cyx=-cxy/L_; cyy=cxx/L_; cyz=0;
    }
  else
    {
      cyx=0; cyy=1; cyz=0;
    }
  czx=cxy*cyz-cxz*cyy;
  czy=-cxx*cyz+cxz*cyx;
  czz=cxx*cyy-cxy*cyx;

  // Transforamtion matrix
  Dense_Matrix T(6,6);
  T(0,0) = cxx; T(0,1) = cxy; T(0,2) = cxz;
  T(1,0) = cyx; T(1,1) = cyy; T(1,2) = cyz;
  T(2,0) = czx; T(2,1) = czy; T(2,2) = czz;

  T(3,3) = cxx; T(3,4) = cxy; T(3,5) = cxz;
  T(4,3) = cyx; T(4,4) = cyy; T(4,5) = cyz;
  T(5,3) = czx; T(5,4) = czy; T(5,5) = czz;

  // Ensure that Material and Properies ptr are set
  if (MaterialPtr == NULL || PropertiesPtr == NULL){
    cerr << endl << " ERROR: Link8::CalcEsm(): " ;
    cerr << "MaterialPtr or PropertiesPtr are set to NULL"  << endl;
    exit(1);
  }

  // Mass matrix
  Dense_Matrix M(6,6);
  M(0,0) = 2; M(1,1) =  2; M(2,2) =  2;
  M(3,3) = 2; M(4,4) =  2; M(5,5) =  2;
  M(3,0) = 1; M(4,1) =  1; M(5,2) =  1;
  M(0,3) = 1; M(1,4) =  1; M(2,5) =  1;

  scale(M, (MaterialPtr->Get("rho") * PropertiesPtr->GetDouble("Area") * L / 6));
  
  Dense_Matrix dummyMatrix(6,6);
  mult_AT_B_A(T,M,dummyMatrix);
  
  return dummyMatrix;
}


////
// Eval length of beam
////
double Link8::EvalLength() const{
  
  return sqrt( pow( NodeVec[1]->Cx - NodeVec[0]->Cx ,2) +
         pow( NodeVec[1]->Cy - NodeVec[0]->Cy ,2) +
         pow( NodeVec[1]->Cz - NodeVec[0]->Cz ,2) ) ; 
}


////
// Evaluate volume of the element
////
double Link8::EvalVolume() const{

  return EvalLength() * PropertiesPtr->GetDouble("Area");
}

////
// Eval 3x3 transformation matrix from global to local coordinates
////
Dense_Matrix Link8::EvalTransformationMatrix(){

  // From local beam coordinate system as defined in ANSYS to global coordinates
  // ---------------------------------------------------------------------------
  Dense_Matrix T(3,3);
  double L = EvalLength();

  // Check if element length axis is within a 0.01 percent slope of global z-axis
  // using direction cosine cxz
  if ( abs( (NodeVec[1]->Cz - NodeVec[0]->Cz) / L )  >= 1/1.000000005 ){
    // 2a) If true, element y-axis is equal global y-axis 
    // mapping local x to global z
    //         local y to global y
    //         local z to -global x
    T(0,2) = 1;
    T(1,1) = 1;
    T(2,0) = -1;

    // Check in which direction on z-axis the element is oriented -> switch signs
    // accordingly ( THIS EFFECT WAS FOUND AFTER 4 DAYS OF DEBUGGING!!!)
    if ( NodeVec[0]->Cz > NodeVec[1]->Cz ) {
      T(0,2) =-1;
      T(2,0) = 1;
    }
  }
  else {
    // 2b) If false, element y-axis is defined to be in global  xy-plane
    //     local x is defined through the nodes
    // mapping to global coords through direction cosines...
    T(0,0) = ( NodeVec[1]->Cx - NodeVec[0]->Cx ) / L;
    T(0,1) = ( NodeVec[1]->Cy - NodeVec[0]->Cy ) / L;
    T(0,2) = ( NodeVec[1]->Cz - NodeVec[0]->Cz ) / L;
    
    // Using that the projection of x to global xy-plane is perpendicular to local y
    double Lx = sqrt( pow( T(0,0) ,2 ) + pow( T(0,1),2 ) );
    T(1,0) = - T(0,1) / Lx;
    T(1,1) = T(0,0) / Lx;
    T(1,2) = 0;
    
    // direction cosines of z direction are evaluated using cross product
    // of x and y vectors (in "direction cosine coordinates")
    T(2,0) = T(0,1) * T(1,2) - T(0,2) * T(1,1);
    T(2,1) = -( T(0,0) * T(1,2) - T(0,2) * T(1,0) );
    T(2,2) = T(0,0) * T(1,1) - T(0,1) * T(1,0);
  }

  return T;
  
}


////
// Eval nodal stresses of beam and store'em in vector< dense_vector > of Element
////
void Link8::EvalStresses(string StressType) {

  //   cout << "Print nodal deformations in global coordinates: " << endl;
  //   print_vector( NodePtr[0]->GetDeformations() );
  //   print_vector( NodePtr[1]->GetDeformations() );

  // Transform displacements from global to local coordinates
  // by multiplication of 3x3 subranges of ui with the transformation matrix
  // ------------------------------------------------------------------------
  Dense_Vector u1(6,0.0), u2(6,0.0);
  Dense_Matrix T = EvalTransformationMatrix();
  mult( T, NodeVec[0]->GetDeformations()(0,3), u1(0,3) );
  mult( T, NodeVec[0]->GetDeformations()(3,6), u1(3,6) );
  mult( T, NodeVec[1]->GetDeformations()(0,3), u2(0,3) );
  mult( T, NodeVec[1]->GetDeformations()(3,6), u2(3,6) );

  //   cout << "Print nodal deformations in local coordinates: " << endl;
  //   print_vector(u1);
  //   print_vector(u2);

  // Define the stress matrix with appropriate size: 2 nodes with 4 values 
  Stresses = Dense_Matrix(2,1);

  // Eval axial stresses of beam and store'em in sx of both nodes
  // --> Axial stress = E * eps_x = E * (u2x - u1x) / L
  // ------------------------------------------------------------------------
  float_type L = EvalLength();
  float_type eps = ( u2[0] - u1[0] ) / L;
  Stresses(0,0) = eps * MaterialPtr->Get("E");
  Stresses(1,0) = Stresses(0,0);
}



////
// Eval Euler buckling
////
double Link8::EvalEulerBuckling( double l_fac){
  EvalStresses("nodal");
  double buckling = 0;
  // get normal stress  
  double s_x = Stresses(0,0);
  // get inertia and area
  double area = GetPropertiesPtr()->GetDouble("Area");
  double pi = 3.14159265359;
  double radius = pow((area/pi),0.5);
  
  // the moments of inertia are calculated for circular links
  double i_y = pi*pow(radius,4)/4;
  
  // the moments of inertia are calculated for square links
  // double i_y  = pow(area,2)/12;
  
  // get length
  double length = EvalLength();
  // get youngs modulus
  double y_mod = GetMaterialPtr()->Get("E");
      
  // critical euler buckling load
  double euler = pow(pi,2)*y_mod/pow(l_fac*length/sqrt(i_y/area ),2);
  if (s_x < 0){
    buckling = euler / abs(s_x) - 1;
  }
  
  return buckling;
}
