//-----------------------------------------------------------------------------
// Beam3.cc
//
// begin     : Dec 6 2001
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
#include "Beam3.h"

using namespace fe_base;

////
//// Initialize static data members of class "Beam3"
////
const int    Beam3::NodeCount		= 2;
const int    Beam3::Id	    		= 3;
const string Beam3::Name	    	= "Structural 2D Beam";
StructDofSet Beam3::ElementDofSet	= StructDofSet("110001");

////
//! Implementation of element stiffness matrix computation in class "Beam3".
//// 
Dense_Matrix Beam3::CalcEM()
{
  // Evaluate length of beam, cos(phi) and sin(phi)
  double L   = sqrt( pow( NodeVec[1]->Cx - NodeVec[0]->Cx ,2) + pow( NodeVec[1]->Cy - NodeVec[0]->Cy ,2) );
  double sin = ( NodeVec[1]->Cy - NodeVec[0]->Cy )/L;
  double cos = ( NodeVec[1]->Cx - NodeVec[0]->Cx )/L;
  
  // Transformation matrix
  Dense_Matrix T(6,6);
  T(0,0) = cos;  T(0,1) = sin;
  T(1,0) = -sin; T(1,1) = cos;
  T(2,2) = 1;
  T(3,3) = cos;  T(3,4) = sin;
  T(4,3) = -sin; T(4,4) = cos;
  T(5,5) = 1;

  // Ensure that Material and Properies ptr are set
  if (MaterialPtr == NULL || PropertiesPtr == NULL){
   cerr << endl << " ERROR: Beam3::CalcEM(): ";
    cerr << "MaterialPtr or PropertiesPtr are set to NULL" << endl;
   exit(1);
  }
  
  // Axial and bending stiffness factors
  double AS = MaterialPtr->Get("E") * PropertiesPtr->GetDouble("Area") / L;  //A*E/L
  double BS = MaterialPtr->Get("E") * PropertiesPtr->GetDouble("Izz")  / pow(L,3); //E*Izz/L^3
  
  double ky = PropertiesPtr->GetDouble("ShearZ");
  double phi_y =  12. * ky * BS * L / (PropertiesPtr->GetDouble("Area") * MaterialPtr->Get("G"));
  double Y1 = 12. * BS / (1. + phi_y);
  double Y2 = 6. * BS * L / (1. + phi_y);
  double Y3 = (4. + phi_y) * BS * pow(L,2.) / (1. + phi_y);
  double Y4 = (2. - phi_y) * BS * pow(L,2.) / (1. + phi_y);

  // Stiffness matrix
  Dense_Matrix K(6,6);
  K(0,0) = AS; 		K(0,3) = -AS;
  K(1,1) = Y1;		K(1,2) = Y2;	 	K(1,4) = -K(1,1); 	K(1,5) = K(1,2);
  K(2,1) = K(1,2);	K(2,2) = Y3;		K(2,4) = -K(1,2); 	K(2,5) = Y4;
  K(3,0) = -AS;		K(3,3) = AS;
  K(4,1) = -K(1,1);	K(4,2) = -K(1,2); 	K(4,4) = K(1,1); 	K(4,5) = -K(1,2);
  K(5,1) = K(1,2);	K(5,2) = K(2,5); 	K(5,4) = -K(1,2); 	K(5,5) = K(2,2);

  // Evaluate Transformation K = T^T * K * T
  Dense_Matrix Kglobal(6,6);
  mult_AT_B_A(T,K,Kglobal);

  return Kglobal;
 
}


////
//! Implementation of element mass matrix computation in class "Beam3".
//// 
Dense_Matrix Beam3::CalcEmm()
{
  // Evaluate length of beam, cos(phi) and sin(phi)
  double L   = sqrt( pow( NodeVec[1]->Cx - NodeVec[0]->Cx ,2) + pow( NodeVec[1]->Cy - NodeVec[0]->Cy ,2) );
  double sin = ( NodeVec[1]->Cy - NodeVec[0]->Cy )/L;
  double cos = ( NodeVec[1]->Cx - NodeVec[0]->Cx )/L;
  
  // Transformation matrix
  Dense_Matrix T(6,6);
  T(0,0) = cos;  T(0,1) = sin;
  T(1,0) = -sin; T(1,1) = cos;
  T(2,2) = 1;
  T(3,3) = cos;  T(3,4) = sin;
  T(4,3) = -sin; T(4,4) = cos;
  T(5,5) = 1;

// Ensure that Material and Properies ptr are set
//  if (MaterialPtr == NULL || PropertiesPtr == NULL){
//    cerr << endl << " ERROR: Beam3::CalcEsm(): ";
//    cerr << "MaterialPtr or PropertiesPtr are set to NULL" << endl;
//    exit(1);
//  }
  
  // Mass of Element
  double m = MaterialPtr->Get("rho") * PropertiesPtr->GetDouble("Area") * L;

  double C = PropertiesPtr->GetDouble("Izz")  / ( pow(L,2) * PropertiesPtr->GetDouble("Area") );

  // Mass matrix
  Dense_Matrix M(6,6);
  M(0,0) = m/3.0;		M(0,3) = m/6.0;
  M(1,1) = (13.0/35.0+6.0/5.0*C)*m;  	M(1,2) = (11.0/210.0+1.0/10.0*C)*L*m ; 	M(1,4) = (9.0/70.0-6.0/5.0*C)*m ; 	M(1,5) = (-13.0/420.0+1.0/10.0*C)*L*m;
  M(2,1) = M(1,2);	M(2,2) = (1.0/105.0+2.0/15.0*C)*L*L*m;	M(2,4) = -M(1,5); 	M(2,5) = (-1.0/140.0-1.0/30.0*C)*L*L*m;
  M(3,0) = m/6.0;		M(3,3) = m/3.0;
  M(4,1) = M(1,4);	M(4,2) = M(2,4); 	M(4,4) = M(1,1); 	M(4,5) = -M(1,2);
  M(5,1) = M(1,5);	M(5,2) = M(2,5); 	M(5,4) = -M(1,2); 	M(5,5) = M(2,2);

  // Evaluate Transformation M = T^T * M * T
  Dense_Matrix Mglobal(6,6);
  mult_AT_B_A(T,M,Mglobal);

  return Mglobal;
 
}
