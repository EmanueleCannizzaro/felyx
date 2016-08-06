//-----------------------------------------------------------------------------
// Solid186.cc
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
#include "Solid186.h"

using namespace fe_base;

////
//// Initialize static data members of class "Solid186"
////
const int    Solid186::NodeCount     	= 20;
const int    Solid186::Id	       	= 186;
const string Solid186::Name	       	= "Structural 3D Solid";
StructDofSet Solid186::ElementDofSet 	= StructDofSet("111000");
Dense_Matrix Solid186::IntPoints        = GetIntegrationPoints(GaussQuadratic3D);


// EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to the standard 
// coords xi, eta, zeta in a specific point.
// Its rows correspond to a specific shape function. The cols to the direction of derivation 
// EvalPoints is a Matrix of points, where the shape functions are to be evaluated. IP determines the specific point,
// means row of Evalpoints
//----------------------------------------------------------------------------------------------------------------------
void Solid186::EvalDerivedShapeFunc
(Dense_Matrix shapefunc, const Dense_Matrix EvalPoints, int  IP) const
{
  shapefunc(0,0)  = -0.125*(1-EvalPoints(IP,1))*(1-EvalPoints(IP,2))*(-2*EvalPoints(IP,0)-EvalPoints(IP,1)-EvalPoints(IP,2)-1);
  shapefunc(1,0)  =  0.125*(1-EvalPoints(IP,1))*(1-EvalPoints(IP,2))*( 2*EvalPoints(IP,0)-EvalPoints(IP,1)-EvalPoints(IP,2)-1);
  shapefunc(2,0)  =  0.125*(1+EvalPoints(IP,1))*(1-EvalPoints(IP,2))*( 2*EvalPoints(IP,0)+EvalPoints(IP,1)-EvalPoints(IP,2)-1);
  shapefunc(3,0)  = -0.125*(1+EvalPoints(IP,1))*(1-EvalPoints(IP,2))*(-2*EvalPoints(IP,0)+EvalPoints(IP,1)-EvalPoints(IP,2)-1);

  shapefunc(4,0)  = -0.125*(1-EvalPoints(IP,1))*(1+EvalPoints(IP,2))*(-2*EvalPoints(IP,0)-EvalPoints(IP,1)+EvalPoints(IP,2)-1);
  shapefunc(5,0)  =  0.125*(1-EvalPoints(IP,1))*(1+EvalPoints(IP,2))*( 2*EvalPoints(IP,0)-EvalPoints(IP,1)+EvalPoints(IP,2)-1);
  shapefunc(6,0)  =  0.125*(1+EvalPoints(IP,1))*(1+EvalPoints(IP,2))*( 2*EvalPoints(IP,0)+EvalPoints(IP,1)+EvalPoints(IP,2)-1);
  shapefunc(7,0)  = -0.125*(1+EvalPoints(IP,1))*(1+EvalPoints(IP,2))*(-2*EvalPoints(IP,0)+EvalPoints(IP,1)+EvalPoints(IP,2)-1);

  shapefunc(8,0)  = -0.5*EvalPoints(IP,0)*(1-EvalPoints(IP,1))*(1-EvalPoints(IP,2));
  shapefunc(9,0)  =  0.25*(1-EvalPoints(IP,1)*EvalPoints(IP,1))*(1-EvalPoints(IP,2));
  shapefunc(10,0) = -0.5*EvalPoints(IP,0)*(1+EvalPoints(IP,1))*(1-EvalPoints(IP,2));
  shapefunc(11,0) = -0.25*(1-EvalPoints(IP,1)*EvalPoints(IP,1))*(1-EvalPoints(IP,2));
  shapefunc(12,0) = -0.5*EvalPoints(IP,0)*(1-EvalPoints(IP,1))*(1+EvalPoints(IP,2));
  shapefunc(13,0) =  0.25*(1-EvalPoints(IP,1)*EvalPoints(IP,1))*(1+EvalPoints(IP,2));
  shapefunc(14,0) = -0.5*EvalPoints(IP,0)*(1+EvalPoints(IP,1))*(1+EvalPoints(IP,2));
  shapefunc(15,0) = -0.25*(1-EvalPoints(IP,1)*EvalPoints(IP,1))*(1+EvalPoints(IP,2));
  shapefunc(16,0) = -0.25*(1-EvalPoints(IP,1))*(1-EvalPoints(IP,2)*EvalPoints(IP,2));
  shapefunc(17,0) =  0.25*(1-EvalPoints(IP,1))*(1-EvalPoints(IP,2)*EvalPoints(IP,2));
  shapefunc(18,0) =  0.25*(1+EvalPoints(IP,1))*(1-EvalPoints(IP,2)*EvalPoints(IP,2));
  shapefunc(19,0) = -0.25*(1+EvalPoints(IP,1))*(1-EvalPoints(IP,2)*EvalPoints(IP,2));
  
  shapefunc(0,1)  = -0.125*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,2))*(-EvalPoints(IP,0)-2*EvalPoints(IP,1)-EvalPoints(IP,2)-1);
  shapefunc(1,1)  = -0.125*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,2))*(+EvalPoints(IP,0)-2*EvalPoints(IP,1)-EvalPoints(IP,2)-1);
  shapefunc(2,1)  =  0.125*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,2))*(+EvalPoints(IP,0)+2*EvalPoints(IP,1)-EvalPoints(IP,2)-1);
  shapefunc(3,1)  =  0.125*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,2))*(-EvalPoints(IP,0)+2*EvalPoints(IP,1)-EvalPoints(IP,2)-1);
  shapefunc(4,1)  = -0.125*(1-EvalPoints(IP,0))*(1+EvalPoints(IP,2))*(-EvalPoints(IP,0)-2*EvalPoints(IP,1)+EvalPoints(IP,2)-1);
  shapefunc(5,1)  = -0.125*(1+EvalPoints(IP,0))*(1+EvalPoints(IP,2))*(+EvalPoints(IP,0)-2*EvalPoints(IP,1)+EvalPoints(IP,2)-1);
  shapefunc(6,1)  =  0.125*(1+EvalPoints(IP,0))*(1+EvalPoints(IP,2))*(+EvalPoints(IP,0)+2*EvalPoints(IP,1)+EvalPoints(IP,2)-1);
  shapefunc(7,1)  =  0.125*(1-EvalPoints(IP,0))*(1+EvalPoints(IP,2))*(-EvalPoints(IP,0)+2*EvalPoints(IP,1)+EvalPoints(IP,2)-1);
  shapefunc(8,1)  =  -0.25*(1-EvalPoints(IP,0)*EvalPoints(IP,0))*(1-EvalPoints(IP,2));
  shapefunc(9,1)  =  -0.5*EvalPoints(IP,1)*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,2));
  shapefunc(10,1) =   0.25*(1-EvalPoints(IP,0)*EvalPoints(IP,0))*(1-EvalPoints(IP,2));
  shapefunc(11,1) =  -0.5*EvalPoints(IP,1)*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,2));
  shapefunc(12,1) =  -0.25*(1-EvalPoints(IP,0)*EvalPoints(IP,0))*(1+EvalPoints(IP,2));
  shapefunc(13,1) =  -0.5*EvalPoints(IP,1)*(1+EvalPoints(IP,0))*(1+EvalPoints(IP,2));
  shapefunc(14,1) =   0.25*(1-EvalPoints(IP,0)*EvalPoints(IP,0))*(1+EvalPoints(IP,2));
  shapefunc(15,1) =  -0.5*EvalPoints(IP,1)*(1-EvalPoints(IP,0))*(1+EvalPoints(IP,2));
  shapefunc(16,1) =  -0.25*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,2)*EvalPoints(IP,2));
  shapefunc(17,1) =  -0.25*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,2)*EvalPoints(IP,2));
  shapefunc(18,1) =   0.25*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,2)*EvalPoints(IP,2));
  shapefunc(19,1) =   0.25*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,2)*EvalPoints(IP,2));
  
  shapefunc(0,2)  = -0.125*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,1))*(-EvalPoints(IP,0)-EvalPoints(IP,1)-2*EvalPoints(IP,2)-1);
  shapefunc(1,2)  = -0.125*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,1))*(+EvalPoints(IP,0)-EvalPoints(IP,1)-2*EvalPoints(IP,2)-1);
  shapefunc(2,2)  = -0.125*(1+EvalPoints(IP,0))*(1+EvalPoints(IP,1))*(+EvalPoints(IP,0)+EvalPoints(IP,1)-2*EvalPoints(IP,2)-1);
  shapefunc(3,2)  = -0.125*(1-EvalPoints(IP,0))*(1+EvalPoints(IP,1))*(-EvalPoints(IP,0)+EvalPoints(IP,1)-2*EvalPoints(IP,2)-1);
  shapefunc(4,2)  =  0.125*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,1))*(-EvalPoints(IP,0)-EvalPoints(IP,1)+2*EvalPoints(IP,2)-1);
  shapefunc(5,2)  =  0.125*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,1))*(+EvalPoints(IP,0)-EvalPoints(IP,1)+2*EvalPoints(IP,2)-1);
  shapefunc(6,2)  =  0.125*(1+EvalPoints(IP,0))*(1+EvalPoints(IP,1))*(+EvalPoints(IP,0)+EvalPoints(IP,1)+2*EvalPoints(IP,2)-1);
  shapefunc(7,2)  =  0.125*(1-EvalPoints(IP,0))*(1+EvalPoints(IP,1))*(-EvalPoints(IP,0)+EvalPoints(IP,1)+2*EvalPoints(IP,2)-1);
  shapefunc(8,2)  = -0.25*(1-EvalPoints(IP,0)*EvalPoints(IP,0))*(1-EvalPoints(IP,1));
  shapefunc(9,2)  = -0.25*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,1)*EvalPoints(IP,1));
  shapefunc(10,2) = -0.25*(1-EvalPoints(IP,0)*EvalPoints(IP,0))*(1+EvalPoints(IP,1));
  shapefunc(11,2) = -0.25*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,1)*EvalPoints(IP,1));
  shapefunc(12,2) =  0.25*(1-EvalPoints(IP,0)*EvalPoints(IP,0))*(1-EvalPoints(IP,1));
  shapefunc(13,2) =  0.25*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,1)*EvalPoints(IP,1));
  shapefunc(14,2) =  0.25*(1-EvalPoints(IP,0)*EvalPoints(IP,0))*(1+EvalPoints(IP,1));
  shapefunc(15,2) =  0.25*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,1)*EvalPoints(IP,1));
  shapefunc(16,2) = -0.5*EvalPoints(IP,2)*(1-EvalPoints(IP,0))*(1-EvalPoints(IP,1));
  shapefunc(17,2) = -0.5*EvalPoints(IP,2)*(1+EvalPoints(IP,0))*(1-EvalPoints(IP,1));
  shapefunc(18,2) = -0.5*EvalPoints(IP,2)*(1+EvalPoints(IP,0))*(1+EvalPoints(IP,1));
  shapefunc(19,2) = -0.5*EvalPoints(IP,2)*(1-EvalPoints(IP,0))*(1+EvalPoints(IP,1));

}

void Solid186::EvalShapeFunc(Dense_Matrix shapefunc, const Dense_Matrix EvalPoints, int  IP)
{
      shapefunc(0,0) = 0.125*(1.0-EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0-EvalPoints(IP,2))*(-EvalPoints(IP,0)-EvalPoints(IP,1)-EvalPoints(IP,2)-2.0);
      shapefunc(1,0) = 0.125*(1.0+EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0-EvalPoints(IP,2))*(+EvalPoints(IP,0)-EvalPoints(IP,1)-EvalPoints(IP,2)-2.0);
      shapefunc(2,0) = 0.125*(1.0+EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0-EvalPoints(IP,2))*(+EvalPoints(IP,0)+EvalPoints(IP,1)-EvalPoints(IP,2)-2.0);
      shapefunc(3,0) = 0.125*(1.0-EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0-EvalPoints(IP,2))*(-EvalPoints(IP,0)+EvalPoints(IP,1)-EvalPoints(IP,2)-2.0);
      shapefunc(4,0) = 0.125*(1.0-EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0+EvalPoints(IP,2))*(-EvalPoints(IP,0)-EvalPoints(IP,1)+EvalPoints(IP,2)-2.0);
      shapefunc(5,0) = 0.125*(1.0+EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0+EvalPoints(IP,2))*(EvalPoints(IP,0)-EvalPoints(IP,1)+EvalPoints(IP,2)-2.0);
      shapefunc(6,0) = 0.125*(1.0+EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0+EvalPoints(IP,2))*(EvalPoints(IP,0)+EvalPoints(IP,1)+EvalPoints(IP,2)-2.0);
      shapefunc(7,0) = 0.125*(1.0-EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0+EvalPoints(IP,2))*(-EvalPoints(IP,0)+EvalPoints(IP,1)+EvalPoints(IP,2)-2.0);

      shapefunc(8,0) = 0.25*(1.0-EvalPoints(IP,0)*EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0-EvalPoints(IP,2));
      shapefunc(9,0) = 0.25*(1.0-EvalPoints(IP,1)*EvalPoints(IP,1))*(1.0+EvalPoints(IP,0))*(1.0-EvalPoints(IP,2));
      shapefunc(10,0) = 0.25*(1.0-EvalPoints(IP,0)*EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0-EvalPoints(IP,2));
      shapefunc(11,0) = 0.25*(1.0-EvalPoints(IP,1)*EvalPoints(IP,1))*(1.0-EvalPoints(IP,0))*(1.0-EvalPoints(IP,2));
      shapefunc(12,0) = 0.25*(1.0-EvalPoints(IP,0)*EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(1.0+EvalPoints(IP,2));
      shapefunc(13,0) = 0.25*(1.0-EvalPoints(IP,1)*EvalPoints(IP,1))*(1.0+EvalPoints(IP,0))*(1.0+EvalPoints(IP,2));
      shapefunc(14,0) = 0.25*(1.0-EvalPoints(IP,0)*EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(1.0+EvalPoints(IP,2));
      shapefunc(15,0) = 0.25*(1.0-EvalPoints(IP,1)*EvalPoints(IP,1))*(1.0-EvalPoints(IP,0))*(1.0+EvalPoints(IP,2));
      shapefunc(16,0) = 0.25*(1.0-EvalPoints(IP,2)*EvalPoints(IP,2))*(1.0-EvalPoints(IP,0))*(1.0-EvalPoints(IP,1));
      shapefunc(17,0) = 0.25*(1.0-EvalPoints(IP,2)*EvalPoints(IP,2))*(1.0+EvalPoints(IP,0))*(1.0-EvalPoints(IP,1));
      shapefunc(18,0) = 0.25*(1.0-EvalPoints(IP,2)*EvalPoints(IP,2))*(1.0+EvalPoints(IP,0))*(1.0+EvalPoints(IP,1));
      shapefunc(19,0) = 0.25*(1.0-EvalPoints(IP,2)*EvalPoints(IP,2))*(1.0-EvalPoints(IP,0))*(1.0+EvalPoints(IP,1));
}


Dense_Matrix Solid186::GetlNodeCoords()
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
