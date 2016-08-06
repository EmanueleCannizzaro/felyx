//-----------------------------------------------------------------------------
// Solid187.cc
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
#include "Solid187.h"

using namespace fe_base;

////
//// Initialize static data members of class "Solid187"
////
const int    Solid187::NodeCount     	= 10;
const int    Solid187::Id	       	= 187;
const string Solid187::Name	       	= "Structural 3D Tetraeder";
StructDofSet Solid187::ElementDofSet 	= StructDofSet("111000");
Dense_Matrix Solid187::IntPoints        = GetIntegrationPoints(VolumeCoords4Points);


// EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to the standard 
// coords xi, eta, zeta in a specific point.
// Its rows correspond to a specific shape function. The cols to the direction of derivation 
// EvalPoints is a Matrix of points, where the shape functions are to be evaluated. IP determines the specific point,
// means row of Evalpoints. The EvalPoints are given in volumetric coords
//----------------------------------------------------------------------------------------------------------------------
void Solid187::EvalDerivedShapeFunc(Dense_Matrix shapefunc, const Dense_Matrix EvalPoints, int  IP) const
{
      shapefunc(0,0) = 1.0 - 4.0*EvalPoints(IP,0);
      shapefunc(1,0) = 4.0*EvalPoints(IP,1) - 1.0;
      shapefunc(2,0) = 0.0;
      shapefunc(3,0) = 0.0;
      shapefunc(4,0) = 4.0*(EvalPoints(IP,0) - EvalPoints(IP,1));
      shapefunc(5,0) = 4.0*EvalPoints(IP,2);
      shapefunc(6,0) = -4.0*EvalPoints(IP,2);
      shapefunc(7,0) = -4.0*EvalPoints(IP,3);
      shapefunc(8,0) = 4.0*EvalPoints(IP,3);
      shapefunc(9,0) = 0.0;

      shapefunc(0,1) = 1.0 - 4.0*EvalPoints(IP,0);
      shapefunc(1,1) = 0.0;
      shapefunc(2,1) = 4.0*EvalPoints(IP,2) -1.0;
      shapefunc(3,1) = 0.0;
      shapefunc(4,1) = -4.0*EvalPoints(IP,1);
      shapefunc(5,1) = 4.0*EvalPoints(IP,1);
      shapefunc(6,1) = 4.0*(EvalPoints(IP,0) - EvalPoints(IP,2));
      shapefunc(7,1) = -4.0*EvalPoints(IP,3);
      shapefunc(8,1) = 0.0;
      shapefunc(9,1) = 4.0*EvalPoints(IP,3);

      shapefunc(0,2) = 1.0 - 4.0*EvalPoints(IP,0);
      shapefunc(1,2) = 0.0;
      shapefunc(2,2) = 0.0;
      shapefunc(3,2) = 4.0*EvalPoints(IP,3) -1.0;
      shapefunc(4,2) = -4.0*EvalPoints(IP,1);
      shapefunc(5,2) = 0.0;
      shapefunc(6,2) = -4.0*EvalPoints(IP,2);
      shapefunc(7,2) = 4.0*(EvalPoints(IP,0) - EvalPoints(IP,3));
      shapefunc(8,2) = 4.0*EvalPoints(IP,1);
      shapefunc(9,2) = 4.0*EvalPoints(IP,2);
}

void Solid187::EvalShapeFunc(Dense_Matrix shapefunc, const Dense_Matrix EvalPoints, int  IP)
{
      shapefunc(0,0) = (2.0*EvalPoints(IP,0)-1.0)*EvalPoints(IP,0);
      shapefunc(1,0) = (2.0*EvalPoints(IP,1)-1.0)*EvalPoints(IP,1);
      shapefunc(2,0) = (2.0*EvalPoints(IP,2)-1.0)*EvalPoints(IP,2);
      shapefunc(3,0) = (2.0*EvalPoints(IP,3)-1.0)*EvalPoints(IP,3);
      shapefunc(4,0) = 4.0*EvalPoints(IP,0)*EvalPoints(IP,1);
      shapefunc(5,0) = 4.0*EvalPoints(IP,1)*EvalPoints(IP,2);
      shapefunc(6,0) = 4.0*EvalPoints(IP,0)*EvalPoints(IP,2);
      shapefunc(7,0) = 4.0*EvalPoints(IP,0)*EvalPoints(IP,3);
      shapefunc(8,0) = 4.0*EvalPoints(IP,1)*EvalPoints(IP,3);
      shapefunc(9,0) = 4.0*EvalPoints(IP,2)*EvalPoints(IP,3);
}


Dense_Matrix Solid187::GetlNodeCoords()
{
  Dense_Matrix lC(10,4);
  lC(0,0) =  1.0; lC(0,1) =  0.0; lC(0,2) =  0.0; lC(0,3) =  0.0;
  lC(1,0) =  0.0; lC(1,1) =  1.0; lC(1,2) =  0.0; lC(1,3) =  0.0;
  lC(2,0) =  0.0; lC(2,1) =  0.0; lC(2,2) =  1.0; lC(2,3) =  0.0;
  lC(3,0) =  0.0; lC(3,1) =  0.0; lC(3,2) =  0.0; lC(3,3) =  1.0;
  lC(4,0) =  0.5; lC(4,1) =  0.5; lC(4,2) =  0.0; lC(4,3) =  0.0;
  lC(5,0) =  0.0; lC(5,1) =  0.5; lC(5,2) =  0.5; lC(5,3) =  0.0;
  lC(6,0) =  0.5; lC(6,1) =  0.0; lC(6,2) =  0.5; lC(6,3) =  0.0;
  lC(7,0) =  0.5; lC(7,1) =  0.0; lC(7,2) =  0.0; lC(7,3) =  0.5;
  lC(8,0) =  0.0; lC(8,1) =  0.5; lC(8,2) =  0.0; lC(8,3) =  0.5;
  lC(9,0) =  0.0; lC(9,1) =  0.0; lC(9,2) =  0.5; lC(9,3) =  0.5;
  return lC;
}
