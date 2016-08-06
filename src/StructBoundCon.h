//-----------------------------------------------------------------------------
// StructBoundCon.h
//
// begin     : Jan 11 2003
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


#ifndef StructBoundCon_h
#define StructBoundCon_h StructBoundCon_h

#include <iostream>
#include <vector>

#include "BoundCon.h"
#include "StructDofSet.h"

#include "mtl/mtl_felyx_utils.h"

namespace fe_base {	//begin of namespace fe_base

  namespace bc {
    //! Enumeration for the different boundary conditions stored in StructBoundCon.BC[].
    enum StructBCtype {Dx, Dy, Dz, Rx, Ry, Rz, Fx, Fy, Fz, Mx, My, Mz};
  }

  //Class StructBoundCon
  class StructBoundCon: public BoundCon<bc::StructBCtype, StructDofSet> {

    public:

      StructBoundCon() : BoundCon<bc::StructBCtype, StructDofSet>( 12 ) {}

      friend ostream& operator<<( ostream&, const StructBoundCon& ); // overloaded << operator

      //! Delete all force boundary conditions
      void deleteForces();

    protected:
      void print( ostream& ) const;		//print-function

  };

}	// end of namespace fe_base

#endif
