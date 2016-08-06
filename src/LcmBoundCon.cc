//-----------------------------------------------------------------------------
// LcmBoundCon.cc
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


#include "LcmBoundCon.h"

namespace fe_base {


  std::ostream& fe_base::operator<< ( std::ostream& stream, const LcmBoundCon& p ) {
    p.print( stream );
    return stream;
  }

  //! Print-function for boundary-conditions.
  void LcmBoundCon::print( std::ostream& stream ) const {
    int count = 0;
    for ( size_t i = 0; i < BCstatus.size(); i++ )
      if ( BCstatus[ i ] ) {
        count++;
      }

    stream.width( 3 );
    stream << count;

    stream.setf( ios::left, ios::adjustfield );
    //stream.precision(9);
    if ( BCstatus[ 0 ] ) {
      stream << "  P: " ;
      stream.width( 6 );
      stream << BCvalue[ 0 ];
    }
    if ( BCstatus[ 1 ] ) {
      stream << "  V: " ;
      stream.width( 6 );
      stream << BCvalue[ 1 ];
    }
  }

  /*! Returns active Loads that consist of non-zero displacement
    constraints and forces. */
  void LcmBoundCon::GetActiveLoads( mtl::Dense_Vector NodalLoads, unsigned minRange, unsigned maxRange ) {

    if ( BCstatus[ 0 ] == 1 && BCvalue[ 0 ] != 0 )
      NodalLoads[ 0 ] = BCvalue[ 0 ];
    else
      NodalLoads[ 0 ] = 0;

    if ( BCstatus[ 1 ] == 1 && BCvalue[ 1 ] != 0 )
      NodalLoads[ 6 ] = BCvalue[ 1 ];
    else
      NodalLoads[ 6 ] = 0;

  }

}
