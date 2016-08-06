//-----------------------------------------------------------------------------
// Solid185.cc
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
#include "Solid185.h"

using namespace fe_base;

////
//// Initialize static data members of class "Solid185"
////
const int    Solid185::NodeCount	= 8;
const int    Solid185::Id	       	= 185;
const string Solid185::Name	       	= "Structural 3D Solid";
StructDofSet Solid185::ElementDofSet 	= StructDofSet("111000");
Dense_Matrix Solid185::IntPoints  = GetIntegrationPoints(GaussQuadratic3D);

//! Implementation of element stiffness matrix computation in class "Solid185".
//---------------------------------------------------------------------------
Dense_Matrix Solid185::CalcEM()
{
  //Initialization of all variables
  double detJ, fdummy;

  //The size of the square ESM is given as follows
  unsigned  nSize = NodeCount*ElementDofSet.count();

  //Initialization of some Matrices 
  Dense_Matrix 
    J(3,3),                     /* 3D - Jacobimatrix */
    B(6,nSize),                 /* 3D - B-Matrix */
    K(nSize,nSize),             /* Stifnessmatrix */
    dummyMatr1(6,nSize),        /* 3D - Dummymatrix */
    dummyMatr2(nSize,nSize),
    shapefunc(NodeCount, ElementDofSet.count()); /* derived - shapefuncs */

  //Get the material data
  Dense_Matrix D = GetMaterialPtr() -> Get(Material::ThreeDSolid);

  //Calculate the function at every Gausspoint.
  for (unsigned i = 0 ; i < IntPoints.nrows(); i++)
  {
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
    //on the element (ansys uses bBar as default integration 
    //method with solid185)
    //----------------------------------------------------------
    SetBMatrix(shapefunc,B,bBarMethod);

    //----------------------------------------------------
    //Computation of the integrand at the ith integration
    //point and addition to the stifnessmatrix
    //----------------------------------------------------

    // detJ = J.pdet(); //compute the determinant of the jacobian
    detJ = lu_det(J); 

    //IntPoints(i,3) is the weight of the ith point.
    //It is included for to be able to change the 
    //degree of integration easely
    fdummy = detJ*IntPoints(i,3); 

     // dummyMatr1 = D*B;
    // B.transpose();
    // dummyMatr2 = B * dummyMatr1;
    // K += dummyMatr2*fdummy;
    // B.transpose(); //...to change the dimension of B 
    mult_AT_B_A(B,D,dummyMatr2);
    add(scaled(dummyMatr2,fdummy),K);


  }

  return K; //The stifnessmatrix of the element....
 
} //of CalcEM()

// EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to the standard 
// coords xi, eta, zeta in a specific point.
// Its rows correspond to a specific shape function. The cols to the direction of derivation 
// EvalPoints is a Matrix of points, where the shape functions are to be evaluated. IP determines the specific point,
// means row of Evalpoints
//----------------------------------------------------------------------------------------------------------------------
void Solid185::EvalDerivedShapeFunc(Dense_Matrix shapefunc,  const Dense_Matrix EvalPoints, int  IP) const
{
  shapefunc(0,0) = -0.125 * (1 - EvalPoints(IP,1)) * (1 - EvalPoints(IP,2));
  shapefunc(1,0) =  0.125 * (1 - EvalPoints(IP,1)) * (1 - EvalPoints(IP,2));
  shapefunc(2,0) =  0.125 * (1 + EvalPoints(IP,1)) * (1 - EvalPoints(IP,2));
  shapefunc(3,0) = -0.125 * (1 + EvalPoints(IP,1)) * (1 - EvalPoints(IP,2));
  shapefunc(4,0) = -0.125 * (1 - EvalPoints(IP,1)) * (1 + EvalPoints(IP,2));
  shapefunc(5,0) =  0.125 * (1 - EvalPoints(IP,1)) * (1 + EvalPoints(IP,2));
  shapefunc(6,0) =  0.125 * (1 + EvalPoints(IP,1)) * (1 + EvalPoints(IP,2));
  shapefunc(7,0) = -0.125 * (1 + EvalPoints(IP,1)) * (1 + EvalPoints(IP,2));
  
  shapefunc(0,1) = -0.125 * (1 - EvalPoints(IP,0)) * (1 - EvalPoints(IP,2));
  shapefunc(1,1) = -0.125 * (1 + EvalPoints(IP,0)) * (1 - EvalPoints(IP,2));
  shapefunc(2,1) =  0.125 * (1 + EvalPoints(IP,0)) * (1 - EvalPoints(IP,2));
  shapefunc(3,1) =  0.125 * (1 - EvalPoints(IP,0)) * (1 - EvalPoints(IP,2));
  shapefunc(4,1) = -0.125 * (1 - EvalPoints(IP,0)) * (1 + EvalPoints(IP,2));
  shapefunc(5,1) = -0.125 * (1 + EvalPoints(IP,0)) * (1 + EvalPoints(IP,2));
  shapefunc(6,1) =  0.125 * (1 + EvalPoints(IP,0)) * (1 + EvalPoints(IP,2));
  shapefunc(7,1) =  0.125 * (1 - EvalPoints(IP,0)) * (1 + EvalPoints(IP,2));
  
  shapefunc(0,2) = -0.125 * (1 - EvalPoints(IP,0)) * (1 - EvalPoints(IP,1));
  shapefunc(1,2) = -0.125 * (1 + EvalPoints(IP,0)) * (1 - EvalPoints(IP,1));
  shapefunc(2,2) = -0.125 * (1 + EvalPoints(IP,0)) * (1 + EvalPoints(IP,1));
  shapefunc(3,2) = -0.125 * (1 - EvalPoints(IP,0)) * (1 + EvalPoints(IP,1));
  shapefunc(4,2) =  0.125 * (1 - EvalPoints(IP,0)) * (1 - EvalPoints(IP,1));
  shapefunc(5,2) =  0.125 * (1 + EvalPoints(IP,0)) * (1 - EvalPoints(IP,1));
  shapefunc(6,2) =  0.125 * (1 + EvalPoints(IP,0)) * (1 + EvalPoints(IP,1));
  shapefunc(7,2) =  0.125 * (1 - EvalPoints(IP,0)) * (1 + EvalPoints(IP,1));  
}

void Solid185::EvalShapeFunc(Dense_Matrix shapefunc, const Dense_Matrix EvalPoints, int  IP)
{
      shapefunc(0,0) = 0.125*(1.0-EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0-EvalPoints(IP,2));
      shapefunc(1,0) = 0.125*(1.0+EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0-EvalPoints(IP,2));
      shapefunc(2,0) = 0.125*(1.0+EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0-EvalPoints(IP,2));
      shapefunc(3,0) = 0.125*(1.0-EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0-EvalPoints(IP,2));
      shapefunc(4,0) = 0.125*(1.0-EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0+EvalPoints(IP,2));
      shapefunc(5,0) = 0.125*(1.0+EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0+EvalPoints(IP,2));
      shapefunc(6,0) = 0.125*(1.0+EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0+EvalPoints(IP,2));
      shapefunc(7,0) = 0.125*(1.0-EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0+EvalPoints(IP,2));
}

Dense_Matrix Solid185::GetlNodeCoords()
{
  Dense_Matrix lC(8,3);
  lC(0,0) = -1.0; lC(0,1) = -1.0; lC(0,2) = -1.0;
  lC(1,0) =  1.0; lC(1,1) = -1.0; lC(1,2) = -1.0;
  lC(2,0) =  1.0; lC(2,1) =  1.0; lC(2,2) = -1.0;
  lC(3,0) = -1.0; lC(3,1) =  1.0; lC(3,2) = -1.0;
  lC(4,0) = -1.0; lC(4,1) = -1.0; lC(4,2) =  1.0;
  lC(5,0) =  1.0; lC(5,1) = -1.0; lC(5,2) =  1.0;
  lC(6,0) =  1.0; lC(6,1) =  1.0; lC(6,2) =  1.0;
  lC(7,0) = -1.0; lC(7,1) =  1.0; lC(7,2) =  1.0;
  return lC;
}
