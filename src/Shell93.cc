//-----------------------------------------------------------------------------
// Shell93.cc
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


#include "Shell93.h"
#include <iostream>

using namespace fe_base;


////
//// Initialize static data members of class "Shell93"
////
const unsigned    Shell93::NodeCount  =  8;
const unsigned    Shell93::Id         = 93;
const string Shell93::Name            = "Structural 3D Shell";
StructDofSet Shell93::ElementDofSet   = StructDofSet("111111");
Dense_Matrix Shell93::IntPoints       = GetIntegrationPoints(GaussQuadratic3D);



//! Implementation of element stiffness matrix computation in class "Shell93".
//---------------------------------------------------------------------------
Dense_Matrix Shell93::CalcEM() {

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

  return SingleLayerESM( GetMaterialPtr(), GetIntPoints(), tmpOrient, GetShellThickness());

} //of CalcEM()

//! Implementation of element mass matrix computation in class "Shell93".
//---------------------------------------------------------------------------
Dense_Matrix Shell93::CalcEmm() {

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
  coords xi, eta, zeta in a specific point.
  Its rows correspond to a specific shape function. The cols to the direction of derivation
  The shape functions are independent of the thickness direction. Therefore Shell93 needs only 2 columns. The original,
  underived shapefunctions are stored in the third column. They are used for different transformations.
*/
//----------------------------------------------------------------------------------------------------------------------
Dense_Matrix Shell93::EvalDerivedShapeFunc( const Dense_Matrix::Row IntPoint ) const{

    float_type ip0 = IntPoint[0], ip1 = IntPoint[1];

    FELYX_LOGIC_ASSERT( abs(ip0) <= 1.0 && abs(ip1) <= 1.0,"Shell93::EvalDerivedShapeFunc()");

    Dense_Matrix DSF(8,3);

    DSF(0,0) = -0.25*(1-ip1)*(-2*ip0-ip1);
    DSF(1,0) =  0.25*(1-ip1)*( 2*ip0-ip1);
    DSF(2,0) =  0.25*(1+ip1)*( 2*ip0+ip1);
    DSF(3,0) = -0.25*(1+ip1)*(-2*ip0+ip1);
    DSF(4,0) = -ip0*(1-ip1);
    DSF(5,0) =  0.5*(1-ip1*ip1);
    DSF(6,0) = -ip0*(1+ip1);
    DSF(7,0) = -0.5*(1-ip1*ip1);

    DSF(0,1) = -0.25*(1-ip0)*(-ip0-2*ip1);
    DSF(1,1) = -0.25*(1+ip0)*( ip0-2*ip1);
    DSF(2,1) =  0.25*(1+ip0)*( ip0+2*ip1);
    DSF(3,1) =  0.25*(1-ip0)*(-ip0+2*ip1);
    DSF(4,1) = -0.5*(1-ip0*ip0);
    DSF(5,1) = -ip1*(1+ip0);
    DSF(6,1) =  0.5*(1-ip0*ip0);
    DSF(7,1) = -ip1*(1-ip0);

    DSF(0,2) = 0.25*(1-ip0)*(1-ip1)*(-ip0-ip1-1);
    DSF(1,2) = 0.25*(1+ip0)*(1-ip1)*(ip0-ip1-1);
    DSF(2,2) = 0.25*(1+ip0)*(1+ip1)*(ip0+ip1-1);
    DSF(3,2) = 0.25*(1-ip0)*(1+ip1)*(-ip0+ip1-1);
    DSF(4,2) = 0.5*(1-ip1)*(1-ip0*ip0);
    DSF(5,2) = 0.5*(1+ip0)*(1-ip1*ip1);
    DSF(6,2) = 0.5*(1+ip1)*(1-ip0*ip0);
    DSF(7,2) = 0.5*(1-ip0)*(1-ip1*ip1);

    return DSF;
}

//this function evaluates a Matrix of deformations derived with respect to the local
//coordinates. The Matrix would be multiplied with a vector of the nodal DOFs to
//receive the derived deformations in natural coordinates
void Shell93::EvalDerivedDeformations(Dense_Matrix DerivedDeformations_, const Dense_Matrix TmpDirCos_,
                                      const Dense_Matrix DerivedShapeFunc_, const Dense_Matrix  EvalPoints_,
                                      unsigned IP){
  if ( abs(EvalPoints_(IP,0)) > 1.0 || abs(EvalPoints_(IP,1)) > 1.0 || abs(EvalPoints_(IP,2)) > 1.0 ) {
    cerr << endl << "****************************************************"
         << endl << "******   Error in Shell93::EvalDerivedShapeFunc     *******"
         << endl << "******   One of the entered natural coords   *******"
         << endl << "******   Is grater than 1.0                  *******"
         << endl << "****************************************************" << endl;

    exit(0);
  }

  else{
    unsigned i;

    mtl::set_value(DerivedDeformations_, 0.0);

    //    double t = PropertiesPtr->GetDouble("Thickness");
    double t = GetShellThickness();

    for ( i = 0 ; i < NodeCount ; ++i ){

      DerivedDeformations_(0,0+5*i) = DerivedShapeFunc_(i, 0);
      DerivedDeformations_(0,3+5*i) = (-1)* EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(1+i*3,0)*0.5;
      DerivedDeformations_(0,4+5*i) = EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(0+i*3,0)*0.5;

      DerivedDeformations_(1,0+5*i) = DerivedShapeFunc_(i, 1);
      DerivedDeformations_(1,3+5*i) = (-1)* EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(1+i*3,0)*0.5;
      DerivedDeformations_(1,4+5*i) = EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(0+i*3,0)*0.5;

      DerivedDeformations_(2,3+5*i) = (-1)* t*DerivedShapeFunc_(i, 2)*TmpDirCos_(1+i*3,0)*0.5;
      DerivedDeformations_(2,4+5*i) = t*DerivedShapeFunc_(i, 2)*TmpDirCos_(0+i*3,0)*0.5;

      DerivedDeformations_(3,1+5*i) = DerivedShapeFunc_(i, 0);
      DerivedDeformations_(3,3+5*i) = (-1)* EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(1+i*3,1)*0.5;
      DerivedDeformations_(3,4+5*i) = EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(0+i*3,1)*0.5;

      DerivedDeformations_(4,1+5*i) = DerivedShapeFunc_(i, 1);
      DerivedDeformations_(4,3+5*i) = (-1)* EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(1+i*3,1)*0.5;
      DerivedDeformations_(4,4+5*i) = EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(0+i*3,1)*0.5;

      DerivedDeformations_(5,3+5*i) = (-1)* t*DerivedShapeFunc_(i, 2)*TmpDirCos_(1+i*3,1)*0.5;
      DerivedDeformations_(5,4+5*i) = t*DerivedShapeFunc_(i, 2)*TmpDirCos_(0+i*3,1)*0.5;

      DerivedDeformations_(6,2+5*i) = DerivedShapeFunc_(i, 0);
      DerivedDeformations_(6,3+5*i) = (-1)* EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(1+i*3,2)*0.5;
      DerivedDeformations_(6,4+5*i) = EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 0)*TmpDirCos_(0+i*3,2)*0.5;

      DerivedDeformations_(7,2+5*i) = DerivedShapeFunc_(i, 1);
      DerivedDeformations_(7,3+5*i) = (-1)* EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(1+i*3,2)*0.5;
      DerivedDeformations_(7,4+5*i) = EvalPoints_(IP, 2)*t*DerivedShapeFunc_(i, 1)*TmpDirCos_(0+i*3,2)*0.5;

      DerivedDeformations_(8,3+5*i) = (-1)* t*DerivedShapeFunc_(i, 2)*TmpDirCos_(1+i*3,2)*0.5;
      DerivedDeformations_(8,4+5*i) = t*DerivedShapeFunc_(i, 2)*TmpDirCos_(0+i*3,2)*0.5;
    }
  }
}

Dense_Matrix Shell93::GetlNodeCoords()
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



