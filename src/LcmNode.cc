//-----------------------------------------------------------------------------
// LcmNode.cc
//
// begin     : Jan 12 2003
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


#include "LcmNode.h"

using namespace fe_base;

////
//// Constructors
////

// Default constructor.
LcmNode::LcmNode()
  : NodeCoordSysPtr(NULL), OrigBoundConPtr(NULL), BoundConPtr(NULL), Fill(0), Pressure(0), Volume(0), q(0), time(0), entrapment(0)
{ NodeDofSet.reset();}

// Specialized constructor.
LcmNode::LcmNode(Dense_Vector v)
  :NodeCoordSysPtr(NULL), OrigBoundConPtr(NULL), BoundConPtr(NULL), Fill(0), Pressure(0), Volume(0), q(0), time(0), entrapment(0)
{
  set(v);
  NodeDofSet.reset();
}

// Specialized constructor
LcmNode::LcmNode(double x, double y, double z)
  :NodeCoordSysPtr(NULL), OrigBoundConPtr(NULL), BoundConPtr(NULL), Fill(0), Pressure(0), Volume(0), q(0), time(0), entrapment(0)
{
  set(x, y, z);
  NodeDofSet.reset();
}

// Copy constructor
LcmNode::LcmNode( const LcmNode& n )
  : Cx(n.Cx),
    Cy(n.Cy),
    Cz(n.Cz),
    NodeCoordSysPtr(n.NodeCoordSysPtr),
    OrigBoundConPtr(n.OrigBoundConPtr),
    BoundConPtr(n.BoundConPtr),
    NodeDofSet(n.NodeDofSet),
    Idx2GM(n.Idx2GM),
    Fill(n.Fill)
{
  // In order to make a "hard" copy of this MTL vector
  Pressure = n.Pressure;
  Volume = n.Volume;
  q = n.q;
  entrapment = n.entrapment;
  time = n.time;
}



////
//// Sets
////
void LcmNode::set(Dense_Vector a){
  Cx = a[0];
  Cy = a[1];
  Cz = a[2];
}

//! Setting node coordinates by three doubles.
void LcmNode::set(double x, double y, double z){
    Cx = x;
    Cy = y;
    Cz = z;
}

//! Setting node coordinates and node coordinate system pointer.
void LcmNode::set(double x, double y, double z, CoordSys* nc){
    set(x,y,z);
    NodeCoordSysPtr = nc;
}

//! Setting node coordinate system pointer.
void LcmNode::set(CoordSys* nc){
    NodeCoordSysPtr = nc;
}

//! Setting boundary condition pointer.
void LcmNode::set(LcmBoundCon* ptr){
  Dense_Vector NodalLoads(12,0);
  ptr->GetActiveLoads(NodalLoads);
  BoundConPtr = ptr;
  LcmBoundCon* ptr2;
  ptr2 = new LcmBoundCon;
  *ptr2 = *ptr;
  OrigBoundConPtr = ptr2;
  //OrigBoundConPtr->set(NodalLoads);
}

//! Setting boundary condition pointer.
void LcmNode::setTempBCptr(LcmBoundCon* ptr){
    BoundConPtr = ptr;
}

void LcmNode::setFill(double x){
    Fill = x;
}

Dense_Vector LcmNode::Get(){

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
bool LcmNode::ExistNodeCoordSys() const{
    return ( NodeCoordSysPtr != NULL );
}

//! Returns true if BoundConPtr != NULL.
bool LcmNode::ExistOrigBoundCon() const{
    return ( OrigBoundConPtr != NULL );
}

//! Returns true if BoundConPtr != NULL.
bool LcmNode::ExistBoundCon() const{
    return ( BoundConPtr != NULL );
}

//! Return active Dof's evaluated from NodeDofSet and homogeneous BC's.
LcmDofSet  LcmNode::GetDofSet() const{
  if (BoundConPtr != NULL)
    return NodeDofSet & BoundConPtr -> GetDofSet();
  else
    return NodeDofSet;
}

// GetGMIndex
// Input:       certain DOF as int: 0-5 -> x,y,z,rx,ry,rz
// Output:      global index in GM, if DOF is active
//              else function returns -1
int LcmNode::GetGMIndex(unsigned dof) const{
  if ( dof < 6 ){

    LcmDofSet ActiveSet = GetDofSet();

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
Dense_Vector LcmNode::GetCoords() const {
  Dense_Vector a(3,0);
  a[0] = Cx;
  a[1] = Cy;
  a[2] = Cz;
  return a;
}

//! Get node coordinate system if exists, else give back (0,0,0).
bool LcmNode::GetNodeCoordSysEuler312( Dense_Vector& a) const{

  bool b = ExistNodeCoordSys();
  if ( b )
    a = NodeCoordSysPtr->Get( CoordSys::Euler312 );
  else
    a[0] = a[1] = a[2] =0.;

  return b;
}

//! Get node coordinate system if exists, else give back (0,0,0).
bool LcmNode::GetNodeCoordSysEuler313( Dense_Vector& a) const{

  bool b = ExistNodeCoordSys();
  if ( b )
    a = NodeCoordSysPtr->Get( CoordSys::Euler313 );
  else
    a[0] = a[1] = a[2] =0.;

  return b;
}

////
// Overload operator == and != in order to apply STL-finds...
////
bool LcmNode::operator==(const LcmNode& N ) const{
  return ( N.Cx == Cx && N.Cy == Cy && N.Cz == Cz );
}
bool LcmNode::operator!=(const LcmNode& N ) const{
  return ( N.Cx != Cx || N.Cy != Cy || N.Cz != Cz );
}

// overloaded '<<'-operator to easily print single nodes
// not member of class LcmNode
ostream& fe_base::operator<< (ostream& stream, const LcmNode& a){

  stream.precision(9);
  stream.setf(ios::right, ios::adjustfield);

  stream.width(16); stream << a.Cx <<" ";
  stream.width(16); stream << a.Cy <<" ";
  stream.width(16); stream << a.Cz <<" ";

  // Print nodal coordinate systems, if any
  if ( a.NodeCoordSysPtr != NULL ){
    stream << *a.NodeCoordSysPtr;
  }

  stream.unsetf(ios::left);

  return stream;
}


