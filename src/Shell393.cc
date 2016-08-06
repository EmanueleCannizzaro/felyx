//-----------------------------------------------------------------------------
// Shell393.cc
//
// begin     : Dec 9 2002
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


#include "Shell393.h"
#include <iostream>

using namespace fe_base;


////
//// Initialize static data members of class "Shell393"
////
const unsigned    Shell393::NodeCount  =  6;
const unsigned    Shell393::Id        = 393;
const string Shell393::Name           = "Structural 3D Shell";
StructDofSet Shell393::ElementDofSet   = StructDofSet("111111");
Dense_Matrix Shell393::IntPoints       = GetIntegrationPoints(Gauss3x2);



//! Implementation of element stiffness matrix computation in class "Shell393".
//---------------------------------------------------------------------------
Dense_Matrix Shell393::CalcEM() {

   double Angle;
   double PI = 4*atan(1.);

  if ( GetMaterialPtr()->ClassName() == "IsotropicMaterial" )
    Angle = 0.0;
  else
    Angle = (GetPropertiesPtr()->GetDouble("Theta"))*PI/180.; //transform to radiants...

  MaterialOrientation tmpOrient(Angle, 0.);

  //This is set before the ESM calcualtion so that also for layered shells
  //only one set of Transformation Matrices is generated
//  SetMaterialTransformationVector();

  return SingleLayerESM( GetMaterialPtr(), GetIntPoints(), tmpOrient, GetShellThickness() );
} //of CalcEM()

//! Implementation of element mass matrix computation in class "Shell393".
//---------------------------------------------------------------------------
Dense_Matrix Shell393::CalcEmm() {

   double Angle;
   double PI = 4*atan(1.);

  if ( GetMaterialPtr()->ClassName() == "IsotropicMaterial" )
    Angle = 0.0;
  else
    Angle = (GetPropertiesPtr()->GetDouble("Theta"))*PI/180.; //transform to radiants...

  MaterialOrientation tmpOrient(Angle, 0.);

  return SingleLayerEMM( GetMaterialPtr(), GetIntPoints(), tmpOrient );

} //of CalcEmm()


/*! EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to the standard
  coords xi, eta, zeta in a specific point. The point itself is given in area-coordinates.
  Its rows correspond to a specific shape function. The cols to the direction of derivation
  The shape functions are independent of the thickness direction. Therefore Shell393 needs only 2 columns. The original,
  underived shapefunctions are stored in the third column. They are used for different transformations.
  EvalPoints is a Matrix of points, where the shape functions are to be evaluated. IP determines the specific point,
  means row of Evalpoints.
*/
//----------------------------------------------------------------------------------------------------------------------
Dense_Matrix Shell393::EvalDerivedShapeFunc( const Dense_Matrix::Row IntPoint ) const{

  float_type ip0 = IntPoint[0], ip1 = IntPoint[1], ip2 = IntPoint[2];

  FELYX_LOGIC_ASSERT( abs(ip0) <= 1.0 && abs(ip1) <= 1.0 && abs(ip2) <= 1.0,"Shell393::EvalDerivedShapeFunc()");

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

  DSF(0,2) = ip0 * (2.0*ip0 - 1.0);
  DSF(1,2) = ip1 * (2.0*ip1 - 1.0);
  DSF(2,2) = ip2 * (2.0*ip2 - 1.0);
  DSF(3,2) = 4.0*ip0*ip1;
  DSF(4,2) = 4.0*ip1*ip2;
  DSF(5,2) = 4.0*ip2*ip0;

  return DSF;
}

void Shell393::EvalShapeFunc(Dense_Matrix DerivedShapeFunc_, const Dense_Matrix EvalPoints_, int  IP){
}

//this function evaluates a Matrix of deformations derived with respect to the local
//coordinates. The Matrix would be multiplied with a vector of the nodal DOFs to
//receive the derived deformations in natural coordinates
void Shell393::EvalDerivedDeformations(Dense_Matrix DerivedDeformations_, const Dense_Matrix TmpDirCos_,
                                      const Dense_Matrix DerivedShapeFunc_, const Dense_Matrix  EvalPoints_,
                                      unsigned IP){
  if ( abs(EvalPoints_(IP,0)) > 1.0 || abs(EvalPoints_(IP,1)) > 1.0 || abs(EvalPoints_(IP,2)) > 1.0 ) {
    cerr << endl << "****************************************************"
         << endl << "******   Error in Shell393::EvalDerivedShapeFunc     *******"
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




Dense_Matrix Shell393::GetlNodeCoords()
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


