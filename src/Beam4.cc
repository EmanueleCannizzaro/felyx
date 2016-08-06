//-----------------------------------------------------------------------------
// Beam4.cc
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
#include "Beam4.h"

extern const double PI = 3.14159265359;

using namespace fe_base;

////
//// Initialize static data members of class "Beam4".
////

const int    Beam4::NodeCount		= 2;
const int    Beam4::Id			= 4;
const string Beam4::Name	    	= "Structural 3D Beam";
StructDofSet Beam4::ElementDofSet	= StructDofSet("111111");

////
// Eval 3x3 transformation matrix from global to local coordinates
////
Dense_Matrix Beam4::EvalTransformationMatrix(){

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

  // If theta != 0: 
  // From general orientation defined through two nodes and Theta to 
  // orientation with x in direction of nodes, y in global xy-plane
  // --------------------------------------------------------------
  if ( PropertiesPtr->GetDouble("Theta") != 0 ){
    double co = cos( PropertiesPtr->GetDouble("Theta") * PI / 180);
    double si = sin( PropertiesPtr->GetDouble("Theta") * PI / 180);
    Dense_Matrix T1(3,3), multT(3,3);
    T1(0,0) = 1.0;
    T1(1,1) = co;  T1(1,2) = si;
    T1(2,1) = -si; T1(2,2) = co;
    
    //    T = T1*T;  
    mult(T1,T,multT);
    copy(multT,T);

  }

  return T;
  
}

////
//! Implementation of element stiffness matrix computation in class "Beam4".
//// 
Dense_Matrix Beam4::CalcEM()
{
  // Ensure that Material and Properies ptr are set
//  if (MaterialPtr == NULL || PropertiesPtr == NULL){
//    cerr << endl << " ERROR: Beam4::CalcEM(): ";
//   cerr << "MaterialPtr or PropertiesPtr are set to NULL" << endl;
//   exit(1);
//  }
  
  // Eval length of beam
  double L = EvalLength();
  
  // Evaluation of the stiffness matrix using local beam coordinates
  // ---------------------------------------------------------------
  // Ordering of coordinates:
  // u1x, u1y, u1z, phi1x, phi1y, phi1z, u2x, u2y, u2z, phi2x, phi2y, phi2z
  Dense_Matrix K(12,12);
  
  // Axial stiffness - x-direction
  double AS = MaterialPtr->Get("E") * PropertiesPtr->GetDouble("Area") / L;  //A*E/L
  double BS_Z = MaterialPtr->Get("E") * PropertiesPtr->GetDouble("Izz")  / pow(L,3); //E*Izz/L^3
  double BS_Y = MaterialPtr->Get("E") * PropertiesPtr->GetDouble("Iyy")  / pow(L,3); //E*Izz/L^3
    
  double ky = PropertiesPtr->GetDouble("ShearZ");
  double kz = PropertiesPtr->GetDouble("ShearY");
  double phi_y =  12. * ky * BS_Z * L / (PropertiesPtr->GetDouble("Area") * MaterialPtr->Get("G"));
  double phi_z =  12. * kz * BS_Y * L / (PropertiesPtr->GetDouble("Area") * MaterialPtr->Get("G"));
  
  double Y1 = 12. * BS_Z / (1. + phi_y);
  double Y2 = 6. * BS_Z * L / (1. + phi_y);
  double Y3 = (4. + phi_y) * BS_Z * pow(L,2.) / (1. + phi_y);
  double Y4 = (2. - phi_y) * BS_Z * pow(L,2.) / (1. + phi_y);
  
  double Z1 = 12. * BS_Y / (1. + phi_z);
  double Z2 = 6. * BS_Y * L / (1. + phi_z);
  double Z3 = (4. + phi_z) * BS_Y * pow(L,2.) / (1. + phi_z);
  double Z4 = (2. - phi_z) * BS_Y * pow(L,2.) / (1. + phi_z);
  
  K(0,0) = AS;  K(0,6) = -AS;
  K(6,0) = -AS; K(6,6) = AS;

  // Bending stiffnes - xy-plane
  K(1,1) = Y1;		K(1,5) = Y2; 		K(1,7) = -K(1,1); 	K(1,11) = K(1,5);
  K(5,1) = K(1,5);	K(5,5) = Y3; 		K(5,7) = -K(1,5); 	K(5,11) = Y4;
  K(7,1) = -K(1,1);	K(7,5) = -K(1,5); 	K(7,7) = K(1,1); 	K(7,11) = -K(1,5);
  K(11,1) = K(1,5);	K(11,5) = K(5,11); 	K(11,7) = K(7,11); 	K(11,11) = K(5,5);

  // Bending stiffnes - xz-plane
  // Pay attention to signs, different from bending stiffness in xy-plane
  K(2,2) = Z1;		K(2,4) = -Z2;	 	K(2,8) = -K(2,2); 	K(2,10) = K(2,4);
  K(4,2) = K(2,4);	K(4,4) = Z3;		K(4,8) = -K(2,4); 	K(4,10) = Z4;
  K(8,2) = -K(2,2);	K(8,4) = -K(2,4); 	K(8,8) = K(2,2); 	K(8,10) = -K(2,4);
  K(10,2) = K(2,4);	K(10,4) = K(4,10); 	K(10,8) = K(8,10); 	K(10,10) = K(4,4);

  // Torsional Stiffnes in Theta-x direction
  double Ipp = PropertiesPtr->GetDouble("Ip");
  if ( Ipp == 0 )
    Ipp = PropertiesPtr->GetDouble("Iyy") + PropertiesPtr->GetDouble("Izz");
  double S = MaterialPtr->Get("E") / ( 1+MaterialPtr->Get("nu") ) / 2 * Ipp / L;
  K(3,3) = S;  K(3,9) = -S;
  K(9,3) = -S; K(9,9) = S;
  
  // Eval coordinate transformation matrix
  Dense_Matrix Tl = EvalTransformationMatrix();

  // Build the final transformation matrix - a 12x12 matrix
  Dense_Matrix T(12,12);

  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      T(i,j) = Tl(i,j);
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      T(i+3,j+3) = Tl(i,j);
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      T(i+6,j+6) = Tl(i,j);
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      T(i+9,j+9) = Tl(i,j);

// Evaluate Transformation K = T^T * K * T
  Dense_Matrix Kglobal(12,12);
  mult_AT_B_A(T,K,Kglobal);

  return Kglobal; 
}

////
//! Implementation of element stiffness matrix computation in class "Beam4".
//// 
Dense_Matrix Beam4::CalcEmm()
{
  // Ensure that Material and Properies ptr are set
//  if (MaterialPtr == NULL || PropertiesPtr == NULL){
//    cerr << endl << " ERROR: Beam4::CalcEsm(): ";
//    cerr << "MaterialPtr or PropertiesPtr are set to NULL" << endl;
//    exit(1);
//  }
  
  // Eval length of beam
  double L = EvalLength();

  // Mass of Element
  double m = MaterialPtr->Get("rho") * PropertiesPtr->GetDouble("Area") * L;
  
  // Evaluation of the stiffness matrix using local beam coordinates
  // ---------------------------------------------------------------
  // Ordering of coordinates:
  // u1x, u1y, u1z, phi1x, phi1y, phi1z, u2x, u2y, u2z, phi2x, phi2y, phi2z
  Dense_Matrix M(12,12);
  
  M(0,0) = 1./3.;  M(0,6) = 1./6.;
  M(6,0) = 1./6.; M(6,6) = 1./3.;

  double s =  PropertiesPtr->GetDouble("Izz") /( PropertiesPtr->GetDouble("Area") * pow(L,2)); 
  M(1,1) = s*6./5.+13./35.;	M(1,5) = L*(s/10.+11./210.); 	M(1,7) = 9./70.-s*6./5.; 	M(1,11) = -L*(13./420.-s/10.);
  M(5,1) = M(1,5);	M(5,5) = pow(L,2)*(s*2./15.+1./105.);  M(5,7) = -M(1,11); 	M(5,11) = -pow(L,2)*(s/30.+1./140.);
  M(7,1) = M(1,7);	M(7,5) = M(5,7); 	M(7,7) = M(1,1); 	M(7,11) = -M(1,5);
  M(11,1) = M(1,11);	M(11,5) = M(5,11); 	M(11,7) = M(7,11); 	M(11,11) = M(5,5);

  s =  PropertiesPtr->GetDouble("Iyy") /( PropertiesPtr->GetDouble("Area") * pow(L,2)); 
  M(2,2) = s*6./5.+13./35.;	M(2,4) = -L*(s/10.+11./210.); 	M(2,8) = 9./70.-s*6./5.; 	M(2,10) = L*(13./420.-s/10.);
  M(4,2) = M(2,4);	M(4,4) = pow(L,2)*(s*2./15.+1./105.);  M(4,8) = -M(2,10); 	M(4,10) = -pow(L,2)*(s/30.+1./140.);
  M(8,2) = M(2,8);	M(8,4) = M(4,8); 	M(8,8) = M(2,2); 	M(8,10) = -M(2,4);
  M(10,2) = M(2,10);	M(10,4) = M(4,10); 	M(10,8) = M(8,10); 	M(10,10) = M(4,4);

  s = (PropertiesPtr->GetDouble("Iyy") + PropertiesPtr->GetDouble("Izz")) /  PropertiesPtr->GetDouble("Area");
  M(3,3) = s/3.;  M(3,9) = s/6.;
  M(9,3) = M(3,9); M(9,9) = M(3,3);

  scale(M,m);
  
  // Eval coordinate transformation matrix
  Dense_Matrix Tl = EvalTransformationMatrix();

  // Build the final transformation matrix - a 12x12 matrix
  Dense_Matrix T(12,12);

  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      T(i,j) = Tl(i,j);
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      T(i+3,j+3) = Tl(i,j);
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      T(i+6,j+6) = Tl(i,j);
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      T(i+9,j+9) = Tl(i,j);

// Evaluate Transformation M = T^T * M * T
  Dense_Matrix Mglobal(12,12);
  mult_AT_B_A(T,M,Mglobal);

  return Mglobal; 
}

////
// Evaluate volume of the element
////
double Beam4::EvalVolume() const{

  // Ensure that Properies ptr is set
//  if ( PropertiesPtr == NULL ){
//    cerr << endl << " ERROR: Beam4::EvalVolume(): ";
 //   cerr << "PropertiesPtr is set to NULL" << endl;
//    exit(1);
//  }

  return EvalLength() * PropertiesPtr->GetDouble("Area");
}

////
// Eval length of beam
////
double Beam4::EvalLength() const{
  
  return sqrt( pow( NodeVec[1]->Cx - NodeVec[0]->Cx ,2) +
	       pow( NodeVec[1]->Cy - NodeVec[0]->Cy ,2) +
	       pow( NodeVec[1]->Cz - NodeVec[0]->Cz ,2) ) ; 
}

////
// Eval Euler buckling
////
double Beam4::EvalEulerBuckling( double l_fac){
  EvalStresses("nodal");
  double buckling = 0;
  // get normal stress  
  double s_x = Stresses(0,0);
  // get inertia and area
  double area = GetPropertiesPtr()->GetDouble("Area");
  double i_y  = GetPropertiesPtr()->GetDouble("Iyy");
  double i_z  = GetPropertiesPtr()->GetDouble("Izz");
  // get length
  double length = EvalLength();
  // get youngs modulus
  double y_mod = GetMaterialPtr()->Get("E");
      
  // critical euler buckling load
  double euler = pow(PI,2)*y_mod/pow(l_fac*length/sqrt(min(i_y,i_z)/area ),2);
  if (s_x < 0){
    buckling = euler / abs(s_x) - 1;
  }
  
  return buckling;
}



////
// Eval nodal stresses of beam and store'em in vector< dense_vector > of Element
////
void Beam4::EvalStresses(string StressType) {

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
  Stresses = Dense_Matrix(2,4);

  // Eval axial stresses of beam and store'em in sx of both nodes
  // --> Axial stress = E * eps_x = E * (u2x - u1x) / L
  // ------------------------------------------------------------------------
  float_type L = EvalLength();
  float_type eps = ( u2[0] - u1[0] ) / L;
  Stresses(0,0) = eps * MaterialPtr->Get("E");
  Stresses(1,0) = Stresses(0,0);

  // Eval bending stresses of the beam at each node  
  // ------------------------------------------------------------------------
  // --> Bending stresses are evaluated using the analytically computed second 
  //     derivative of the cubic shape function of this beam.

  // --> Bending stress in xy-plane  
  /*     sigma_b_y = E * eps_b_y = E * (-1) * y * v''(x)
	 with y: max Y-Coord = ThicknessY/2 
	      v: Shape function for bending in xy-plane 
	      x: insert x=0 for the first node, x=L for the second node
  */
  // Bending stresses in xy-plane, y=max, x=0 (--> first node)
  eps = PropertiesPtr->GetDouble("ThicknessY") * ( -3* u2[1] + 3*u1[1] + 2*u1[5]*L + u2[5]*L  ) / pow(L,2);
  Stresses(0,1) = eps * MaterialPtr->Get("E");
  // Bending stresses in xy-plane, y=max, x=L (--> second node)
  eps = PropertiesPtr->GetDouble("ThicknessY") * -( -3* u2[1] + 3*u1[1] + u1[5]*L + 2*u2[5]*L  ) / pow(L,2);
  Stresses(1,1) = eps * MaterialPtr->Get("E");

  // --> Bending stress in xz-plane  
  /*     sigma_b_z = E * eps_b_z = E * (-1) * z * v''(x)
	 with z: max Z-Coord = ThicknessZ/2 
	      v: Shape function for bending in xy-plane 
	      x: insert x=0 for the first node, x=L for the second node
	 !! Watch the signs for the phi_y values !!
  */
  // Bending stresses in xz-plane, z=max, x=0 (--> first node)
  eps = PropertiesPtr->GetDouble("ThicknessZ") * ( -3* u2[2] + 3*u1[2] - 2*u1[4]*L - u2[4]*L  ) / pow(L,2);
  Stresses(0,2) = eps * MaterialPtr->Get("E");
  // Bending stresses in xz-plane, z=max, x=L (--> second node)
  eps = PropertiesPtr->GetDouble("ThicknessZ") * -( -3* u2[2] + 3*u1[2] - u1[4]*L - 2*u2[4]*L  ) / pow(L,2);
  Stresses(1,2) = eps * MaterialPtr->Get("E");

  // Eval maximum torsional shear stresses - only gives valid values for axysymmetric cross-sections
  // --> tau_phi_x = - G * r_max * (phi2 - phi1) / L
  // -----------------------------------------------------------------------------------------------
  double rmax = max( PropertiesPtr->GetDouble("ThicknessZ") , PropertiesPtr->GetDouble("ThicknessY") ) /2 ;
  eps = - rmax * ( u2[3] - u1[3] ) / L ;
  Stresses(0,3) = eps * MaterialPtr->Get("G");
  Stresses(1,3) = Stresses(0,3);

  // Print
  // cout << " Stress vectors of element: " << endl;
  // print_vector( Stresses[0] );
  // print_vector( Stresses[1] );

}
