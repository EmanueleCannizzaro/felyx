//-----------------------------------------------------------------------------
// Plane183.cc
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
#include "Plane183.h"

using namespace fe_base;

////
//// Initialize static data members of class "Plane183"
////
const int    Plane183::NodeCount  	= 8;
const int    Plane183::Id		= 183;
const string Plane183::Name	    	= "Structural 2D Plane";
StructDofSet Plane183::ElementDofSet 	= StructDofSet("110000");
Dense_Matrix Plane183::IntPoints        = GetIntegrationPoints(GaussQuadratic2D);



// EvalDerivedShapeFunc is a matrix of the values of the derived shape functions with respect to the standard 
// coords xi and eta in a specific point.
// Its rows correspond to a specific shape function. The cols to the direction of derivation 
// EvalPoints is a Matrix of points, where the shape functions are to be evaluated. IP determines the specific point,
// means row of Evalpoints
//----------------------------------------------------------------------------------------------------------------------
void Plane183::EvalDerivedShapeFunc(Dense_Matrix shapefunc, const Dense_Matrix EvalPoints, int  IP) const
{
  shapefunc(0,0) = -0.25*(1.0-EvalPoints(IP,1))*(-2.0*EvalPoints(IP,0)-EvalPoints(IP,1));
  shapefunc(1,0) =  0.25*(1.0-EvalPoints(IP,1))*( 2.0*EvalPoints(IP,0)-EvalPoints(IP,1));
  shapefunc(2,0) =  0.25*(1.0+EvalPoints(IP,1))*( 2.0*EvalPoints(IP,0)+EvalPoints(IP,1));
  shapefunc(3,0) = -0.25*(1.0+EvalPoints(IP,1))*(-2.0*EvalPoints(IP,0)+EvalPoints(IP,1));
  shapefunc(4,0) = -EvalPoints(IP,0)*(1.0-EvalPoints(IP,1));
  shapefunc(5,0) = 0.5*(1.0-EvalPoints(IP,1)*EvalPoints(IP,1));
  shapefunc(6,0) = -EvalPoints(IP,0)*(1.0+EvalPoints(IP,1));
  shapefunc(7,0) = -0.5*(1.0-EvalPoints(IP,1)*EvalPoints(IP,1));

  shapefunc(0,1) = -0.25*(1.0-EvalPoints(IP,0))*(-EvalPoints(IP,0)-2.0*EvalPoints(IP,1));
  shapefunc(1,1) = -0.25*(1.0+EvalPoints(IP,0))*( EvalPoints(IP,0)-2.0*EvalPoints(IP,1));
  shapefunc(2,1) =  0.25*(1.0+EvalPoints(IP,0))*( EvalPoints(IP,0)+2.0*EvalPoints(IP,1));
  shapefunc(3,1) =  0.25*(1.0-EvalPoints(IP,0))*(-EvalPoints(IP,0)+2.0*EvalPoints(IP,1));
  shapefunc(4,1) = -0.5*(1.0-EvalPoints(IP,0)*EvalPoints(IP,0));
  shapefunc(5,1) = -EvalPoints(IP,1)*(1.0+EvalPoints(IP,0));
  shapefunc(6,1) =  0.5*(1.0-EvalPoints(IP,0)*EvalPoints(IP,0));
  shapefunc(7,1) = -EvalPoints(IP,1)*(1.0-EvalPoints(IP,0));
}

void Plane183::EvalShapeFunc(Dense_Matrix shapefunc, const Dense_Matrix EvalPoints, int  IP)
{
  shapefunc(0,0) = 0.25*(1.0-EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(-EvalPoints(IP,0)-EvalPoints(IP,1)-1.0);
  shapefunc(1,0) = 0.25*(1.0+EvalPoints(IP,0))*(1.0-EvalPoints(IP,1))*(EvalPoints(IP,0)-EvalPoints(IP,1)-1.0);
  shapefunc(2,0) = 0.25*(1.0+EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(EvalPoints(IP,0)+EvalPoints(IP,1)-1.0);
  shapefunc(3,0) = 0.25*(1.0-EvalPoints(IP,0))*(1.0+EvalPoints(IP,1))*(-EvalPoints(IP,0)+EvalPoints(IP,1)-1.0);
  shapefunc(4,0) = 0.5*(1.0-EvalPoints(IP,0)*EvalPoints(IP,0))*(1.0-EvalPoints(IP,1));
  shapefunc(5,0) = 0.5*(1.0+EvalPoints(IP,0))*(1.0-EvalPoints(IP,1)*EvalPoints(IP,1));
  shapefunc(6,0) = 0.5*(1.0-EvalPoints(IP,0)*EvalPoints(IP,0))*(1.0+EvalPoints(IP,1));
  shapefunc(7,0) = 0.5*(1.0-EvalPoints(IP,0))*(1.0-EvalPoints(IP,1)*EvalPoints(IP,1));
}


Dense_Matrix Plane183::GetlNodeCoords()
{
  Dense_Matrix lC(8,3);
  lC(0,0) = -1.0; lC(0,1) = -1.0;
  lC(1,0) =  1.0; lC(1,1) = -1.0;
  lC(2,0) =  1.0; lC(2,1) =  1.0;
  lC(3,0) = -1.0; lC(3,1) =  1.0;
  lC(4,0) =  0.0; lC(4,1) = -1.0;
  lC(5,0) =  1.0; lC(5,1) =  0.0;
  lC(6,0) =  0.0; lC(6,1) =  1.0;
  lC(7,0) = -1.0; lC(7,1) =  0.0;
  return lC;
}
