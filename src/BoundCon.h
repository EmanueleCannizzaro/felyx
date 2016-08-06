//-----------------------------------------------------------------------------
// BoundCon.h
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

#ifndef BoundCon_h
#define BoundCon_h BoundCon_h

#include <iostream>
#include <vector>

#include "mtl/mtl_felyx_utils.h"
#include "FelyxException.h"

namespace fe_base {

  //! Abstract boundary condition base class
  template <class BCtype, class dof_type>
  class BoundCon {

    public:

      typedef BCtype enumType;

      //! Base constructor initializing member data
      BoundCon( const size_t& BCcount )
          : BCstatus( BCcount, false ), BCvalue( BCcount, 0.0 ) {}

      virtual ~BoundCon() {}

      //!Interface that assigns a value to a certain boundary condition.
      void set
        ( BCtype type, double value ) {
        BCstatus[ type ] = true;
        BCvalue[ type ] = value;
      }

      //!Interface that assigns a value to a certain boundary condition.
      void set
        ( mtl::Dense_Vector vec ) {
        FELYX_RUNTIME_ASSERT( vec.size() <= BCvalue.size(), "BoundCon::set(Dense_Vector vec): vector sizes do not match" );
        for ( int i = 0; i < vec.size(); i++ ) {
          BCvalue[ i ] = vec[ i ];
          BCstatus[ 0 ] = true;
        }
      }

      //! Returns the value of a certain boundary condition, doesn't check if the BC is used.
      double getValue( BCtype type ) {
        return BCvalue[ type ];
      }

      //! Returns the status of a boundary condition, used or not.
      bool getStatus( size_t i ) {
        return BCstatus[ i ];
      }

      /*! Returns the complete BC. If the BC is used the function returns true and the double
        is overwritten, if the BC is not used the function returns false and double is left unchanged. */
      bool get
        ( BCtype type, double& value ) {
        if ( BCstatus[ type ] ) {
          value = BCvalue[ type ];
          return true;
        } else
          return false;
      }

      /*! Returns active Dofs from point of view of BC's
        (from BC's only homogeneous Dofs are inactive). */
      dof_type GetDofSet() {
        dof_type ActiveSet;
        ActiveSet.set();		// set all bits to 1

        for ( int i = 0; i < 6; i++ ) {
          if ( ( BCstatus[ i ] == 1 ) && ( BCvalue[ i ] == 0 ) )
            ActiveSet[ 5 - i ] = 0;
        }
        return ActiveSet;
      }

      //! Returns active Loads, they consist of non-zero displacement constraints and forces
      virtual void GetActiveLoads( mtl::Dense_Vector NodalLoads, unsigned minRange = 0, unsigned maxRange = 12 ) {
        FELYX_RUNTIME_ASSERT( maxRange <= BCstatus.size(), "BoundCon::GetActiveLoads() out of range index access");

        for ( unsigned i = minRange; i < maxRange; i++ ) {
          if ( BCstatus[ i ] == 1 && BCvalue[ i ] != 0 )
            NodalLoads[ i - minRange ] = BCvalue[ i ];
          else
            NodalLoads[ i - minRange ] = 0;
        }
      }

      //! Delete a certain boundary condition
      void deleteBC( BCtype type ) {
        BCstatus[ type ] = false;
        BCvalue[ type ] = 0.0;
      }


      //! Overload == operator
      bool operator==( const BoundCon<BCtype, dof_type>& BC ) {
        return ( BCstatus == BC.BCstatus && BCvalue == BC.BCvalue );
      }

    protected:
      virtual void print( ostream& ) const = 0;

      //!Determines which boundary conditions are active.
      std::vector<bool> BCstatus;
      //!Values of the boundary conditions.
      std::vector<double> BCvalue;

    private:
      BoundCon(){}

  };

}	// end of namespace fe_base

#endif
