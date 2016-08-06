/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// PostProcessing.h
//
// begin     : July 8 2002
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
   

#ifndef PostProcessing_h
#define PostProcessing_h PostProcessing_h

#include "ElementHeaders.h"
#include "mtl/mtl_felyx_utils.h"
#include "PtrVector.h"
#include <vector>

namespace fe_base{

  //! Store deformations of solution vector in nodes
  template<class SolVector>
  void StoreNodalDeformations( vector<StructNode>& Nodes_, const SolVector& NodalDeformations_ ){
    vector<StructNode>::iterator nodeit;
    for (nodeit = Nodes_.begin(); nodeit != Nodes_.end(); ++nodeit )
      nodeit->SetDeformations( NodalDeformations_ );    
  }
  
  //! Get maximum deformation
  float_type GetMaximumDeformation( const vector<StructNode>& Nodes_ ); 

  //! Compute compliance of structure based on nodal loads and deformations
  float_type EvalCompliance( const std::vector<StructNode>& Nodes );


  //! Eval element stress vectors
  void EvalElementStresses(PtrVector<StructElement*>& Elements, string StressType="");

  //! Eval element strains
  void EvalStrains(PtrVector<StructElement*>& Elements);
  
  //! Eval element strains of layered elements for specific layer
  void EvalStrains( PtrVector<StructElement*>& Elements, StructElement::LayerEnum layer );
  
  //! Eval max von Mises Stresses for Nodes
  float_type EvalMaxVonMisesStress( const vector<StructNode>& Nodes_ );   
  
  //! Maximum von Mises strain of the elements
  float_type evalMaxVonMisesStrain( const PtrVector<StructElement*>& Elements );
  
  /*! Eval max failure criteria. Available types: "MaximumStress", "Tsai-Wu", "Tsai-Hill",
        and "Hashin".
  */
  float_type EvalMaxFailureCriteria( std::string type, const PtrVector<StructElement*>& Elements );
  
  //! Eval total mass of active elements in the actual model
  float_type EvalMass(const PtrVector<StructElement*>& Elements);
  
  //! Eval von Mises value for 2D; size of vector v must be 3
  template <class Vector3>
  float_type vonMises2D( const Vector3 v){
    FELYX_LOGIC_ASSERT( v.size() == 3, "vonMises2D( Vector3 v )" )
    
    return sqrt( pow(v[0],2) + pow(v[1],2) - v[0]*v[1] + 3*pow(v[2],2) );
  }
  
  //! Eval von Mises value for 3D; size of vector v must be 6
  template <class Vector6>
  float_type vonMises3D( const Vector6 v){
    FELYX_LOGIC_ASSERT( v.size() == 6, "vonMises3D( Vector6 v )" )
    return 1/sqrt(2.0) * sqrt( pow(v[0]-v[1],2) + pow(v[1]-v[2],2) + pow(v[2]-v[0],2) + 6*( pow(v[3],2) + pow(v[4],2) + pow(v[5],2) ) );
  }
  
  //! Eval vonMises value for any vector v of size 3 or 6 (not needed at the moment)
  /*
  template <class aVector>
  float_type vonMises( const aVector v ){
    FELYX_LOGIC_ASSERT( v.size() == 6 || v.size() == 3 , "vonMises( aVector v )" )
    if ( v.size() == 3 )
      return vonMises2D( v );
    
    return vonMises3D( v );
  }
  */
  
} // of namespace


#endif
