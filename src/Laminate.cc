//-----------------------------------------------------------------------------
// Laminate.cc
//
// begin     : July 2004
// copyright : (c) 2004 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
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
#include "Laminate.h"
#include "FelyxException.h"
#include <typeinfo>

namespace fe_base {

double Laminate::GetLaminateThickness() const {
  double t = 0.0;
  for (unsigned i = 0 ; i < size() ; ++i )
    t += operator[](i)->GetThickness();
  return t;
}

double Laminate::GetTopInterface( const unsigned l ) const{

  FELYX_LOGIC_ASSERT( l < size() , " Laminate::GetTopInterface( l ) " );

  double inter = -(GetLaminateThickness()/2.0);
  for (unsigned i = 0 ; i <= l ; ++i )
    inter += operator[](i)->GetThickness();

  return inter;
}

double Laminate::GetBottomInterface( const unsigned l ) const{

  FELYX_LOGIC_ASSERT( l < size() , " Laminate::GetBottomInterface( l ) " );

  double inter = -(GetLaminateThickness()/2.0);
  for (unsigned i = 0 ; i < l ; ++i )
      inter += operator[](i)->GetThickness();

  return inter;
}


void Laminate::printOn(ostream& stream) const{
  stream << "Laminate: nlayer= " << size() << " total thickness = " << GetLaminateThickness() << std::endl;
  for ( unsigned i = 0 ; i < size() ; ++i ){
    stream.width(10); stream << i << " ";
    operator[](i)->printOn( stream );
    stream << std::endl;
  }
}

mtl::Dense_Matrix Laminate::GetAMatrix( std::vector<MaterialOrientation> ) const
{   //! returns the normalized laminate stiffness matrix A*
    mtl::Dense_Matrix A(6,6),Q(6,6);
    mtl::Dense_Matrix rotation(6,6), rotation_T(6,6), rotation_T_inv(6,6) , Q_global(6,6);
    double t = 0, h_layer, factor, PI = 4*atan(1.);

    for (unsigned i = 0; i < size(); ++i) {
    //  get material properties
      mtl::set_value(Q_global,0.);
      OrthotropicMaterial *ortho  = dynamic_cast<OrthotropicMaterial*>(operator[](i)->GetMaterialPtr());
      double E1 = ortho->OrthotropicMaterial::Get("E1"),
             E2 = ortho->OrthotropicMaterial::Get("E2"),
             nu12 = ortho->OrthotropicMaterial::Get("nu12");
      double G12 = ortho->OrthotropicMaterial::Get("G12"),
             G23 = ortho->OrthotropicMaterial::Get("G23")*5/6,
             G31 = ortho->OrthotropicMaterial::Get("G31")*5/6;

    // stiffnesmatrix Q in local coordinates
      h_layer = operator[](i)->GetThickness();
      t += h_layer;
      factor = h_layer/(1. - E2/E1*nu12*nu12);
      Q(0,0) = E1*factor; Q(1,1) = E2*factor;
      Q(0,1) = Q(1,0) = nu12*E2*factor;
      Q(3,3) = G12*h_layer;  //with shear change to Q(3,3)
      Q(4,4) = G23*h_layer; Q(5,5) = G31*h_layer;

    // transformation: local into global
      //rotation_T_inv = GetPlaneMaterialTransformationMatrix((operator[](i)->GetAngle())*PI/180.);
      rotation = Get3DMaterialTransformationMatrix((operator[](i)->GetAngle())*PI/180.);
      mtl::transpose(rotation, rotation_T);
      mtl::lu_inversion(rotation_T, rotation_T_inv);
      mult_AT_B_A(rotation_T_inv,Q,Q_global);

    // thickness
      mtl::add(Q_global,A);     // multiply stiffness matrix with layer thickness
      }

      mtl::scale(A,1./t);
    return A;   // units of the A Matrix are N/mm
}

mtl::Dense_Matrix Laminate::GetABDMatrix( std::vector<MaterialOrientation>) const
{   //! returns the laminate stiffness matrix ABD
    mtl::Dense_Matrix ABD(6,6), A(3,3), B(3,3), D(3,3), Q(3,3), rotation(3,3) , Q_global(3,3);
    double PI = 4*atan(1.);

    for (unsigned i = 0; i < size(); ++i) {
      // get material properties
      mtl::set_value(Q_global,0.);

      OrthotropicMaterial *ortho  = dynamic_cast<OrthotropicMaterial*>(operator[](i)->GetMaterialPtr());

      double E1 = ortho->OrthotropicMaterial::Get("E1"), E2 = 0., nu12 = 0., G12 = 0.;
             E2 = ortho->OrthotropicMaterial::Get("E2"),
             nu12 = ortho->OrthotropicMaterial::Get("nu12"),
             G12 = ortho->OrthotropicMaterial::Get("G12");
      double factor = 1./(1. - E2/E1*nu12*nu12);

      // stiffnessmatrix Q in local coordinates
      double h_layer = operator[](i)->GetThickness();
      double h_top = GetTopInterface(i);
      double h_bottom = GetBottomInterface(i);
      Q(0,0) = E1*factor; Q(1,1) = E2*factor;
      Q(0,1) = Q(1,0) = nu12*E2*factor;
      Q(2,2) = G12;

      // transform Q into global coordinates
      rotation = GetPlaneMaterialTransformationMatrix((operator[](i)->GetAngle())*PI/180.);
      mult_AT_B_A(rotation,Q,Q_global);

      // A, B and D matrix
      mtl::add(scaled(Q_global, h_layer), A);
      mtl::add(scaled(Q_global, 0.5*(h_top*h_top - h_bottom*h_bottom)), B);
      mtl::add(scaled(Q_global, 1./3.*(h_top*h_top*h_top - h_bottom*h_bottom*h_bottom)), D);
      }

    ABD(0,0) = A(0,0); ABD(0,1) = ABD(1,0) = A(0,1); ABD(0,2) = ABD(2,0) = A(0,2);
    ABD(1,1) = A(1,1); ABD(1,2) = ABD(2,1) = A(1,2); ABD(2,2) = A(2,2);
    ABD(0,3) = ABD(3,0) = B(0,0);
    ABD(0,4) = ABD(1,3) = ABD(3,1) = ABD(4,0) = B(0,1);
    ABD(0,5) = ABD(2,3) = ABD(3,2) = ABD(5,0) = B(0,2);
    ABD(1,4) = ABD(4,1) = B(1,1);
    ABD(1,5) = ABD(2,4) = ABD(4,2) = ABD(5,1) = B(1,2);
    ABD(2,5) = ABD(5,2) = B(2,2);
    ABD(3,3) = D(0,0); ABD(3,4) = ABD(4,3) = D(0,1); ABD(3,5) = ABD(5,3) = D(0,2);
    ABD(4,4) = D(1,1); ABD(4,5) = ABD(5,4) = D(1,2); ABD(5,5) = D(2,2);

    return ABD;         // units of the ABD Matrix are N/mm, N and mm/N
}

ostream& fe_base::operator<< (ostream& stream, const Laminate& L){
  L.print(stream);

  return stream;
}

void Laminate::print(std::ostream& stream) const
{
  for (unsigned i = 0 ; i < size() ; ++i )
  {
    stream << std::endl << "Layer " << i << " Thickness: " << operator[](i)->GetThickness()
           << " Angle [deg]: " << operator[](i)->GetAngle()
           << " MaterialType : " << operator[](i)->GetMaterialPtr()-> ClassName() << std::endl;
  }
}

};
