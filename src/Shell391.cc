//-----------------------------------------------------------------------------
// Shell391.cc
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
   

#include "Shell391.h"
#include <iostream>

using namespace fe_base;

////
//// Initialize static data members of class "Shell391"
////
const unsigned    Shell391::NodeCount  =  6;
const unsigned    Shell391::Id	        = 391;
const string Shell391::Name	        = "Structural 3D Shell";
StructDofSet Shell391::ElementDofSet    = StructDofSet("111111");
Dense_Matrix Shell391::IntPoints        = GetIntegrationPoints(Gauss3x3);

// Returns the matrix of integration points for a specific layer number
Dense_Matrix Shell391::GetIntPoints(const unsigned layerNr ) const{
  
  FELYX_RUNTIME_ASSERT( layerNr >=0 && layerNr < GetLayerCount() && GetLayerThickness( layerNr ) >= 0.0 , "Shell91::GetIntPoints( const unsigned layerNr )" )

  //get the top and bottom interface of the current layer in global  dimensions
  double zTop = LaminatePtr->GetTopInterface(layerNr);
  double zBottom = LaminatePtr->GetBottomInterface(layerNr);

  //transforme them to natural coordinates
  double tLam = LaminatePtr->GetLaminateThickness();
  zTop *= 2.0/tLam;
  zBottom *= 2.0/tLam;

  //the center and the thickness of the layer in natural coords
  double zCenter = (zTop + zBottom)/2.0;
  double tLay = zTop - zBottom;

  // Create a new matrix for integration points of this layer
  Dense_Matrix IntPointsOfLayer(IntPoints.nrows(), IntPoints.ncols() );
  mtl::copy(IntPoints, IntPointsOfLayer);

  //adjust the matrix  of integration points zeta coordinate
  //to integrate only in the active layer...  
  for (unsigned l = 0 ; l < IntPointsOfLayer.nrows() ; ++l ){
    IntPointsOfLayer(l,3) *= tLay*0.5;
    IntPointsOfLayer(l,3) += zCenter;
    IntPointsOfLayer(l,4) *= tLay*0.5; 
  }

  return IntPointsOfLayer;
}


/*! EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to the standard 
  coords xi, eta, zeta in a specific point.
  Its rows correspond to a specific shape function. The cols to the direction of derivation
  The shape functions are independent of the thickness direction. Therefore Shell391 needs only 2 columns. The original,
  underived shapefunctions are stored in the third column. They are used for different transformations.
  EvalPoints is a Matrix of points, where the shape functions are to be evaluated. IP determines the specific point,
  means row of Evalpoints.
*/
//----------------------------------------------------------------------------------------------------------------------
Dense_Matrix Shell391::EvalDerivedShapeFunc( const Dense_Matrix::Row IntPoint ) const{
  
  float_type ip0 = IntPoint[0], ip1 = IntPoint[1], ip2 = IntPoint[2];
  
  FELYX_LOGIC_ASSERT( abs(ip0) <= 1.0 && abs(ip1) <= 1.0 && abs(ip2) <= 1.0,"Shell391::EvalDerivedShapeFunc()");
  
  Dense_Matrix DSF(6,3);    
  DSF(0,0) = 1.0 - 4.0*ip0;
  DSF(1,0) = 4.0*ip1 - 1.0;
  DSF(2,0) = 0;
  DSF(3,0) = 4.0*(ip0 - ip1);
  DSF(4,0) = 4.0*ip2;
  DSF(5,0) = -4.0*ip2;

  DSF(0,1) = 1.0 - 4.0*ip0;
  DSF(1,1) = 0;
  DSF(2,1) = 4.0*ip2 -1.0;
  DSF(3,1) = -4.0*ip1;
  DSF(4,1) = 4.0*ip1;
  DSF(5,1) = 4.0*(ip0 - ip2);

  DSF(0,2) = ip0 * (2*ip0 - 1.0);
  DSF(1,2) = ip1 * (2*ip1 - 1.0);
  DSF(2,2) = ip2 * (2*ip2 - 1.0);
  DSF(3,2) = 4*ip0*ip1;
  DSF(4,2) = 4*ip1*ip2;
  DSF(5,2) = 4*ip2*ip0;

  return DSF;
}


//this function evaluates a Matrix of deformations derived with respect to the local 
//coordinates. The Matrix would be multiplied with a vector of the nodal DOFs to 
//receive the derived deformations in natural coordinates
void Shell391::EvalDerivedDeformations(Dense_Matrix DerivedDeformations_, const Dense_Matrix TmpDirCos_,
				      const Dense_Matrix DerivedShapeFunc_, const Dense_Matrix  EvalPoints_,
				      unsigned IP){
  if ( abs(EvalPoints_(IP,0)) > 1.0 || abs(EvalPoints_(IP,1)) > 1.0 || abs(EvalPoints_(IP,2)) > 1.0 ) {
    cerr << endl << "****************************************************" 
	 << endl << "******   Error in Shell391::EvalDerivedShapeFunc*******" 
	 << endl << "******   One of the entered natural coords   *******" 
	 << endl << "******   Is grater than 1.0                  *******" 
	 << endl << "****************************************************" << endl;
    
    exit(0);
  }
  else{
    unsigned i;
    double zetaCoord = EvalPoints_(IP,3);
    
    mtl::set_value(DerivedDeformations_, 0.0);
        
    double t = GetShellThickness();

    for ( i = 0 ; i < NodeCount ; ++i ){
      
      DerivedDeformations_(0,0+5*i) = DerivedShapeFunc_(i, 0);
      DerivedDeformations_(0,3+5*i) = (-1)* zetaCoord*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(1+i*3,0)*0.5;
      DerivedDeformations_(0,4+5*i) = zetaCoord*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(0+i*3,0)*0.5;
      
      DerivedDeformations_(1,0+5*i) = DerivedShapeFunc_(i, 1);
      DerivedDeformations_(1,3+5*i) = (-1)* zetaCoord*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(1+i*3,0)*0.5;
      DerivedDeformations_(1,4+5*i) = zetaCoord*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(0+i*3,0)*0.5;
      
      DerivedDeformations_(2,3+5*i) = (-1)* t*DerivedShapeFunc_(i, 2)*TmpDirCos_(1+i*3,0)*0.5;
      DerivedDeformations_(2,4+5*i) = t*DerivedShapeFunc_(i, 2)*TmpDirCos_(0+i*3,0)*0.5;
      
      DerivedDeformations_(3,1+5*i) = DerivedShapeFunc_(i, 0);
      DerivedDeformations_(3,3+5*i) = (-1)* zetaCoord*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(1+i*3,1)*0.5;
      DerivedDeformations_(3,4+5*i) = zetaCoord*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(0+i*3,1)*0.5;
      
      DerivedDeformations_(4,1+5*i) = DerivedShapeFunc_(i, 1);
      DerivedDeformations_(4,3+5*i) = (-1)* zetaCoord*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(1+i*3,1)*0.5;
      DerivedDeformations_(4,4+5*i) = zetaCoord*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(0+i*3,1)*0.5;
      
      DerivedDeformations_(5,3+5*i) = (-1)* t*DerivedShapeFunc_(i, 2)*TmpDirCos_(1+i*3,1)*0.5;
      DerivedDeformations_(5,4+5*i) = t*DerivedShapeFunc_(i, 2)*TmpDirCos_(0+i*3,1)*0.5;
      
      DerivedDeformations_(6,2+5*i) = DerivedShapeFunc_(i, 0);
      DerivedDeformations_(6,3+5*i) = (-1)* zetaCoord*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(1+i*3,2)*0.5;
      DerivedDeformations_(6,4+5*i) = zetaCoord*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(0+i*3,2)*0.5;
      
      DerivedDeformations_(7,2+5*i) = DerivedShapeFunc_(i, 1);
      DerivedDeformations_(7,3+5*i) = (-1)* zetaCoord*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(1+i*3,2)*0.5;
      DerivedDeformations_(7,4+5*i) = zetaCoord*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(0+i*3,2)*0.5;
      
      DerivedDeformations_(8,3+5*i) = (-1)* t*DerivedShapeFunc_(i, 2)*TmpDirCos_(1+i*3,2)*0.5;
      DerivedDeformations_(8,4+5*i) = t*DerivedShapeFunc_(i, 2)*TmpDirCos_(0+i*3,2)*0.5;
    }
  }
} 


    
Dense_Matrix Shell391::GetlNodeCoords()
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


