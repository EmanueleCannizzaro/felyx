/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// PostProcessing.cc
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


#include "PostProcessing.h"

#include "ElementHeaders.h"
#include "mtl/mtl_felyx_utils.h"
#include <vector>

using namespace std;
namespace fe_base{


  //! Get maximum deformation
  float_type GetMaximumDeformation( const vector<StructNode>& Nodes_ ) {
    float_type maxdef=0.0, def=0.0;
    std::vector<StructNode>::const_iterator nodeit;
   for (nodeit = Nodes_.begin(); nodeit != Nodes_.end(); ++nodeit ){
      def = mtl::two_norm( (nodeit->GetDeformations())(0,3) );
      if (def > maxdef)
        maxdef = def;
    }
   return maxdef;
  }

  //! Compute compliance, based on nodal loads!
  float_type EvalCompliance( const std::vector<StructNode>& Nodes ){
    float_type compliance=0.0;
    Dense_Vector forces(6,0.0);

    std::vector<StructNode>::const_iterator nodeit;
    for (nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
      if ( nodeit->ExistBoundCon() ){
        nodeit->BoundConPtr->GetActiveLoads( forces,6,12);
        compliance+= mtl::dot( nodeit->GetDeformations(), forces );
      }
    }
    return compliance;
  }

  //! Eval element stress vectors
  void EvalElementStresses(PtrVector<StructElement*>& Elements, string StressType){
    for (PtrVector<StructElement*>::iterator it=Elements.begin(); it !=Elements.end(); ++it )
      { (*it)->EvalStresses(StressType);}
  }

  //! Maximum nodal von Mises stress
  float_type EvalMaxVonMisesStress( const std::vector<StructNode>& Nodes ){
    float_type max=0.0, s=0.0;
    for ( std::vector<StructNode>::const_iterator nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
      s = nodeit->GetVMStress();
      if ( s > max )
        max = s;
    }
    return max;
  }


  //! Maximum von Mises strain of the elements
  float_type evalMaxVonMisesStrain( const PtrVector<StructElement*>& Elements ){
    float_type max=0.0, e=0.0;
    // Loop over elements
    for (PtrVector<StructElement*>::const_iterator it=Elements.begin(); it !=Elements.end(); ++it ){
      e = (*it)->evalMaxVonMisesStrain();
      if ( e > max ){
        max = e;
      }
    }
    return max;
  }

  //! Maximum failure criteria
  float_type EvalMaxFailureCriteria( std::string type, const PtrVector<StructElement*>& Elements){
    float_type max=0.0, e=0.0;

    // Loop over elements
    for (PtrVector<StructElement*>::const_iterator it=Elements.begin(); it !=Elements.end(); ++it ){
      e = (*it)->evalMaxFailureCriteria( type );
      if ( e > max ){
        max = e;
      }
    }
    return max;
    
  }
  
  //! Eval total mass of active elements in the actual model
  float_type EvalMass(const PtrVector<StructElement*>& Elements){
    double mass = 0.0;
    PtrVector<StructElement*>::const_iterator eleit;
    for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ){
      // Only take active elements
      if ( (*eleit)->GetStatus() )
        mass += (*eleit)->EvalMass();
    }
    return mass;
  }

} // of namespace
