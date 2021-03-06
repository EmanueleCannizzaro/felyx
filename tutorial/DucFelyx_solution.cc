///-----------------------------------------------------------------------------
// DucFelyx.cc
//
// begin     : Jan 2005
// copyright : (c) 2005 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {koenig,wintermantel,zehnder}@even-ag.ch
// www       : even-ag.ch
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

#include "DucFelyx_solution.h"
#include <stdexcept>

////
// Eval mass of frame in kg
////
double DucFelyxObject::EvalMass() const {

  double mass = 0.0;
  PtrVector<element_type*>::const_iterator eleit;
  for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ) {

    // only take elements with reference number 1 or 2
    if ( ( *eleit ) ->GetRefNumber() >= 1 && ( *eleit ) ->GetRefNumber() <= 2 )
      mass += ( *eleit ) ->EvalMass();
  }

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed()
    << " >-- " << "Mass of frame [kg] : " << mass * 1000 << endl;

  // Return mass of frame in kg
  return mass * 1000;
}

////
// Delete all forces and moments and apply torsion loadcase
////
void DucFelyxObject::ApplyTorsionLoadcase() {

  // Delete all consisting forces and moments
  for ( vector<bc_type>::iterator it = BoundaryConditions.begin();
        it != BoundaryConditions.end(); ++it )
    it->deleteForces();

  // Apply torsional moment
  vector<node_type>::iterator nit = find ( Nodes.begin(), Nodes.end(), node_type( 785.0 , 0.0 , 135.0 ) );
  if ( nit->ExistBoundCon() )
    nit->BoundConPtr->set( Mx, 700000 );
  else
    throw std::logic_error( "ERROR in DucFelyxObject::ApplyTorsionLoadcase : BoundConPtr not set!" );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed()
    << " >-- " << " Applied torsion loadcase to node : " << distance( Nodes.begin(), nit ) << endl;

}

////
// Delete all forces and moments and apply braking loadcase
////
void DucFelyxObject::ApplyBrakingLoadcase() {

  // Delete all consisting forces and moments
  for ( vector<bc_type>::iterator it = BoundaryConditions.begin();
        it != BoundaryConditions.end(); ++it )
    it->deleteForces();

  // Apply Forces that introduce the braking moment
  double Force = 4256.0e3 / 182.5;
  vector<node_type>::iterator nit = find ( Nodes.begin(), Nodes.end(), node_type( 785.0 , 0.0 , 44.0 ) );
  // Set boundary condition pointer
  nit->set(new fe_base::StructBoundCon() );
  
  if ( nit->ExistBoundCon() )
    nit->BoundConPtr->set( Fx, -Force );
  else
    throw std::logic_error( "ERROR in DucFelyxObject::ApplyBrakingLoadcase : BoundConPtr not set1!" );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed()
    << " >-- " << " Applied braking loadcase to node : " << distance( Nodes.begin(), nit ) << endl;

  nit = find ( Nodes.begin(), Nodes.end(), node_type( 785.0 , 0.0 , 226.5 ) );
  nit->set(new fe_base::StructBoundCon() );
  if ( nit->ExistBoundCon() )
    nit->BoundConPtr->set( Fx, Force );
  else
    throw std::logic_error( "ERROR in DucFelyxObject::ApplyBrakingLoadcase : BoundConPtr not set2!" );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed()
    << " >-- " << " and node : " << distance( Nodes.begin(), nit ) << endl;

}




////
// Eval torsion stiffness of frame, including engine in Nm/degree
////
double DucFelyxObject::EvalTorsionStiffness() const {

  vector<node_type>::const_iterator nodeit1, nodeit2, nodeit;
  mtl::dense1D<float_type> v1( 3, 0 ), v2( 3, 0 );

  // Find iterators for the end nodes of the rear fork side of the frame
  nodeit1 = find ( Nodes.begin(), Nodes.end(), node_type( 0.0 , 134.0 , 0.0 ) );
  nodeit2 = find ( Nodes.begin(), Nodes.end(), node_type( 0.0 , -134.0 , 0.0 ) );

  // Eval the vector giving the directon of rear fork axis in deformed shape
  // v1 = C2 + U2 - (C1 + U1)
  mtl::add
    ( nodeit2->GetCoords(),
        nodeit2->GetDeformations(),
        v1 );

  mtl::add
    ( mtl::scaled( nodeit1->GetCoords(), -1 ),
        mtl::scaled( nodeit1->GetDeformations(), -1 ),
        v1,
        v1 );

  // Find iterators for the end nodes of the steering head of the frame
  nodeit1 = find ( Nodes.begin(), Nodes.end(), node_type( 785.0 , 0.0 , 44.0 ) );
  nodeit2 = find ( Nodes.begin(), Nodes.end(), node_type( 785.0 , 0.0 , 226.5 ) );

  // Eval the vector giving the directon of rear fork axis in deformed shape
  // v2 = C2 + U2 - (C1 + U1)
  mtl::add
    ( nodeit2->GetCoords(),
        nodeit2->GetDeformations(),
        v2 );

  mtl::add
    ( mtl::scaled( nodeit1->GetCoords(), -1 ),
        mtl::scaled( nodeit1->GetDeformations(), -1 ),
        v2,
        v2 );

  // Twist of frame in degrees for actual load
  double twist = acos ( mtl::dot ( v1, v2 ) / ( mtl::two_norm( v1 ) * two_norm( v2 ) ) );

  // Twist minus the 90 degree of undeformed shape, in degrees
  twist = abs( PI / 2 - twist ) * 180 / PI;

  // Eval applied torsion moment on frame
  // ------------------------------------
  // Get iterator to node where moment is applied
  nodeit = find ( Nodes.begin(), Nodes.end(), node_type( 785.0 , 0.0 , 135.0 ) );

  if ( !nodeit->ExistBoundCon() )
    throw std::logic_error( "ERROR in DucFelyxObject::EvalTorsionStiffnes : BoundaryCondition set used to evaluate torsional stiffness is missing!" );

  if ( !nodeit ->BoundConPtr->getStatus( Mx ) )
    throw std::logic_error( "ERROR in DucFelyxObject::EvalTorsionStiffnes : Mx not defined for the BCset used to evaluate torsional stiffness!" );

  double Moment = nodeit->BoundConPtr->getValue( Mx );
  double TorsionStiffness = Moment / 1000.0 / twist;

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed()
    << " >-- " << "Torsion stiffness of frame [ Nm/degree ] : " << TorsionStiffness << endl;

  // Eval torsion stiffness in Nm / degree
  return TorsionStiffness;
}

////
// Eval max stress of frame, in a range of element reference numbers
////
double DucFelyxObject::EvalMaxStress() {

  // Eval stress vectors of elements
  EvalStresses();

  unsigned n = 0;
  float_type sb_max = 0.0, sax_max = 0.0, sn_max = 0.0, s_mises = 0.0, s_mises_max = 0.0;

  // Loop through all elements
  PtrVector<StructElement*>::const_iterator eleit_begin, eleit;
  eleit_begin = Elements.begin();
  for ( eleit = eleit_begin; eleit != Elements.end(); ++eleit ) {

    // Only look at elements in certain range of element reference numbers
    if ( ( *eleit ) ->GetRefNumber() == 1 || ( *eleit ) ->GetRefNumber() == 2 ) {


      // Loop through stress vectors of element
      for ( n = 0; n < ( *eleit ) ->Stresses.nrows(); ++n ) {

        // Eval axial stress -> using abs(s_ax) to always get a max stress
        sax_max = abs ( ( *eleit ) ->Stresses( n, 0 ) );

        // Eval max bending stress in model
        // equals vector addition of s_b_y_max and s_b_z_max
        sb_max = sqrt( pow( ( *eleit ) ->Stresses( n, 1 ), 2 ) +
                       pow( ( *eleit ) ->Stresses( n, 2 ), 2 ) );

        // Gives a total normal stress
        sn_max = sax_max + sb_max;

        // Eval equivalent maximal stress after von Mises by including torsional shear stress
        s_mises = sqrt( pow( sn_max, 2 ) + 3 * pow( ( *eleit ) ->Stresses( n, 3 ) , 2 ) );

        if ( s_mises > s_mises_max )
          s_mises_max = s_mises;

        if ( noise > 1 ) {
          OutStream << "Stresses at node " << setw( 3 ) << n + 1
          << " of elem " << setw( 5 ) << distance( eleit_begin, eleit ) + 1
          << " : sax_max=" << setw( 13 ) << sax_max
          << " sb_max=" << setw( 13 ) << sb_max
          << " sn_max=" << setw( 13 ) << sn_max
          << " s_mises=" << setw( 13 ) << s_mises << endl;
        }
      }
    }
  }
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed()
    << " >-- " << "Maximum von Mises stress in frame [ N/mm2 ] : " << s_mises_max << endl;

  return s_mises_max;
}

////
// Eval and set Beam properties for a certain Ri and t
////
void DucFelyxObject::SetTubeProperties( unsigned tube,
                                        double Ri,
                                        double t,
                                        bool final ) {
  // Map tube index to iterator
  vector<PropertySet>::iterator propit = Properties.begin();
  if ( tube < nVarTubes )
    propit += tube;
  else
    throw std::logic_error( "ERROR in DucFelyxObject::SetTubeProperties : Asked for the modification of a tube that is not within allowed set!" );

  // eval area's of beams and apply'em to FE-model
  propit->Set( "Area", PI * ( pow( Ri + t, 2 ) - pow( Ri, 2 ) ) );

  // eval moment's of inertia for this beams
  double I = PI * ( pow( Ri + t, 4 ) - pow( Ri, 4 ) ) / 4;

  propit->Set( "Izz", I );
  propit->Set( "Iyy", I );
  propit->Set( "Ip", 2 * I );

  // store outer radius in ThicknessZ and ThicknessY - needed for bending stress evaluations
  propit->Set( "ThicknessZ", Ri + t );
  propit->Set( "ThicknessY", Ri + t );

}
