//-----------------------------------------------------------------------------
// StructNode.cc
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
   n
   You should have received a copy of the GNU General Public License
   along with FELyX; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
//-----------------------------------------------------------------------------


#include "StructNode.h"

using namespace fe_base;

////
//// Constructors
////

// Default constructor.
StructNode::StructNode()
  : NodeCoordSysPtr(NULL), BoundConPtr(NULL), counter(0)
{ NodeDofSet.reset(); }

// Specialized constructor.
StructNode::StructNode(Dense_Vector v)
  :NodeCoordSysPtr(NULL), BoundConPtr(NULL), counter(0)
{
  set(v);
  NodeDofSet.reset();
}

// Specialized constructor
StructNode::StructNode(double x, double y, double z)
  :NodeCoordSysPtr(NULL), BoundConPtr(NULL), counter(0)
{
  set(x, y, z);
  NodeDofSet.reset();
}

// Copy constructor
StructNode::StructNode( const StructNode& n )
  : Cx(n.Cx),
    Cy(n.Cy),
    Cz(n.Cz),
    NodeCoordSysPtr(n.NodeCoordSysPtr),
    BoundConPtr(n.BoundConPtr),
    NodeDofSet(n.NodeDofSet),
    Idx2GM(n.Idx2GM)
{
  // In order to make a "hard" copy of this MTL vector
  copy( n.Deformations, Deformations);
}



////
//// Sets
////
void StructNode::set(Dense_Vector a){
  Cx = a[0];
  Cy = a[1];
  Cz = a[2];
}

//! Setting node coordinates by three doubles.
void StructNode::set(double x, double y, double z){
    Cx = x;
    Cy = y;
    Cz = z;
}

//! Setting node coordinates and node coordinate system pointer.
void StructNode::set(double x, double y, double z, CoordSys* nc){
    set(x,y,z);
    NodeCoordSysPtr = nc;
}

//! Setting node coordinate system pointer.
void StructNode::set(CoordSys* nc){
    NodeCoordSysPtr = nc;
}

//! Setting boundary condition pointer.
void StructNode::set(StructBoundCon* ptr){
    BoundConPtr = ptr;
}

Dense_Vector StructNode::Get(){

  Dense_Vector bla(3,0.0);
  bla[0] = Cx;
  bla[1] = Cy;
  bla[2] = Cz;

  return bla;
}



////
//// Queries
////

//! Returns true if NodeCoordSysPtr != NULL.
bool StructNode::ExistNodeCoordSys() const{
    return ( NodeCoordSysPtr != NULL );
}

//! Returns true if BoundConPtr != NULL.
bool StructNode::ExistBoundCon() const{
    return ( BoundConPtr != NULL );
}

//! Get node coordinate system if exists, else give back (0,0,0).
bool StructNode::GetNodeCoordSysEuler312( Dense_Vector& a) const{

  bool b = ExistNodeCoordSys();
  if ( b )
    a = NodeCoordSysPtr->Get( CoordSys::Euler312 );
  else
    a[0] = a[1] = a[2] =0.;

  return b;
}

//! Get node coordinate system if exists, else give back (0,0,0).
bool StructNode::GetNodeCoordSysEuler313( Dense_Vector& a) const{

  bool b = ExistNodeCoordSys();
  if ( b )
    a = NodeCoordSysPtr->Get( CoordSys::Euler313 );
  else
    a[0] = a[1] = a[2] =0.;

  return b;
}

//! Return active Dof's evaluated from NodeDofSet and homogeneous BC's.
StructDofSet  StructNode::GetDofSet() const{
  if (BoundConPtr != NULL)
    return NodeDofSet & (BoundConPtr->GetDofSet());
  else
    return NodeDofSet;
}

// GetGMIndex
// Input:       certain DOF as int: 0-5 -> x,y,z,rx,ry,rz
// Output:      global index in GM, if DOF is active
//              else function returns -1
int StructNode::GetGMIndex(unsigned dof) const{
  if ( dof < 6 ){

    StructDofSet ActiveSet = GetDofSet();

    if ( ActiveSet.test(5-dof) )               // is this dof active?
      // Return Idx2GM + number of bits set til "dof" - 1
      return Idx2GM + ActiveSet.count(dof)-1;
    else
      return -1;
  }
  else
    return -1;
}

// Get coordinates of nodes in MTL-vector
Dense_Vector StructNode::GetCoords() const {
  Dense_Vector a(3,0);
  a[0] = Cx;
  a[1] = Cy;
  a[2] = Cz;
  return a;
}

Dense_Vector StructNode::GetTranslations() const
{ Dense_Vector Transl(3,0.);
  Transl[0] = Deformations[0];
  Transl[1] = Deformations[1];
  Transl[2] = Deformations[2];
  return Transl;
}

////
// Overload operator == and != in order to apply STL-finds...
////
bool StructNode::operator==(const StructNode& N ) const{
  return ( N.Cx == Cx && N.Cy == Cy && N.Cz == Cz );
}
bool StructNode::operator!=(const StructNode& N ) const{
  return ( N.Cx != Cx || N.Cy != Cy || N.Cz != Cz );
}

// overloaded '<<'-operator to easily print single nodes
ostream& fe_base::operator<< (ostream& stream, const StructNode& a){

  stream.precision(9);
  stream.setf(ios::right, ios::adjustfield);

  stream.width(16); stream << a.Cx <<" ";
  stream.width(16); stream << a.Cy <<" ";
  stream.width(16); stream << a.Cz <<" ";


  //Print the pointer to the Nodal Coordsys
  //stream << endl << "The pointer to the nodal cs: "
  //       << a.NodeCoordSysPtr << endl;

  // Print nodal coordinate systems, if any
  if ( a.NodeCoordSysPtr != NULL ){
    stream << *a.NodeCoordSysPtr;
  }

  stream.unsetf(ios::left);

  return stream;
}


float_type StructNode::GetVMStress() const {
  float_type VMStress=0;
  if (Stresses.size()==4)
    VMStress = sqrt(0.5*((Stresses[0]-Stresses[1])*(Stresses[0]-Stresses[1]) + Stresses[0]*Stresses[0] +  Stresses[1]*Stresses[1] + 6.0*Stresses[2]*Stresses[2] ));

  if (Stresses.size()==7)
    VMStress = sqrt(0.5*((Stresses[0]-Stresses[1])*(Stresses[0]-Stresses[1]) + (Stresses[0]-Stresses[2])*(Stresses[0]-Stresses[2]) + (Stresses[2]-Stresses[1])*(Stresses[2]-Stresses[1]) + 6.0*(Stresses[3]*Stresses[3] + Stresses[4]*Stresses[4] + Stresses[5]*Stresses[5]  )));

  return VMStress;
}
