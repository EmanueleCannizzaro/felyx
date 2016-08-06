//-----------------------------------------------------------------------------
// Link1.cc
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
#include "Link1.h"

using namespace fe_base;

////
//// Initialize static data members of class "Link1"
////
const int    Link1::NodeCount		= 2;
const int    Link1::Id	    		= 1;
const string Link1::Name	    	= "Structural 2D Link";
StructDofSet Link1::ElementDofSet	= StructDofSet("110000");

////
//// Implementation of element stiffness matrix computation in class "Link1"
//// 
Dense_Matrix Link1::CalcEM(){
  
  // Evaluate length of link, cos(phi) and sin(phi)
  double L   = sqrt( pow( NodeVec[1]->Cx - NodeVec[0]->Cx ,2) + pow( NodeVec[1]->Cy - NodeVec[0]->Cy ,2) );
  double sin = ( NodeVec[1]->Cy - NodeVec[0]->Cy )/L;
  double cos = ( NodeVec[1]->Cx - NodeVec[0]->Cx )/L;
  
  // Transformation matrix
  Dense_Matrix T(4,4);
  T(0,0) = cos;  T(0,1) = sin;
  T(1,0) = -sin; T(1,1) = cos;
  T(2,2) = cos;  T(2,3) = sin;
  T(3,2) = -sin; T(3,3) = cos;

  // Ensure that Material and Properies ptr are set
  if (MaterialPtr == NULL || PropertiesPtr == NULL){
    cerr << endl << " ERROR: Link1::CalcEM(): " ;
    cerr << "MaterialPtr or PropertiesPtr are set to NULL"  << endl;
    exit(1);
  }
  
  // Axial stiffness
  double AS = MaterialPtr->Get("E") * PropertiesPtr->GetDouble("Area") / L;

  // Stiffness matrix
  Dense_Matrix K(4,4);
  K(0,0) = 1; 		K(0,2) = -1;
  K(2,0) =-1;		K(2,2) =  1;

  //  K *= AS;
  scale(K, AS);  

  // Evaluate Transformation K = T^T * K * T
  Dense_Matrix Kglobal(4,4);
  mult_AT_B_A(T,K,Kglobal);

  return Kglobal; 
}

////
//// Implementation of element stiffness matrix computation in class "Link1"
//// 
Dense_Matrix Link1::CalcEmm(){
  
  // Evaluate length of link, cos(phi) and sin(phi)
  double L   = sqrt( pow( NodeVec[1]->Cx - NodeVec[0]->Cx ,2) + pow( NodeVec[1]->Cy - NodeVec[0]->Cy ,2) );
  double sin = ( NodeVec[1]->Cy - NodeVec[0]->Cy )/L;
  double cos = ( NodeVec[1]->Cx - NodeVec[0]->Cx )/L;
  
  // Transformation matrix
  Dense_Matrix T(4,4);
  T(0,0) = cos;  T(0,1) = sin;
  T(1,0) = -sin; T(1,1) = cos;
  T(2,2) = cos;  T(2,3) = sin;
  T(3,2) = -sin; T(3,3) = cos;

  // Ensure that Material and Properies ptr are set
  if (MaterialPtr == NULL || PropertiesPtr == NULL){
    cerr << endl << " ERROR: Link1::CalcEsm(): " ;
    cerr << "MaterialPtr or PropertiesPtr are set to NULL"  << endl;
    exit(1);
  }
  
  double m = MaterialPtr->Get("rho") * PropertiesPtr->GetDouble("Area") * L / 6;

  // Mass matrix
  Dense_Matrix M(4,4);
  M(0,0) = 2; M(1,1) = 2; M(2,2) = 2; M(3,3) = 2; 
  M(0,2) = 1; M(1,3) = 1; M(2,0) = 1; M(3,1) = 1;

  scale(M, m);  

  Dense_Matrix Mglobal(4,4);
  mult_AT_B_A(T,M,Mglobal);

  return Mglobal; 
}
