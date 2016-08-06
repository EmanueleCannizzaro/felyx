//-----------------------------------------------------------------------------
// StructBoundCon.cc
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


#include "StructBoundCon.h"

namespace fe_base {

  void StructBoundCon::deleteForces() {
    for ( int i = 5; i < 12; ++i ) {
      BCstatus[ i ] = false;
      BCvalue[ i ] = 0.0;
    }
  }

  std::ostream& fe_base::operator<< ( std::ostream& stream, const StructBoundCon& p ) {
    p.print( stream );
    return stream;
  }

  //! Print-function for boundary-conditions.
  void StructBoundCon::print( std::ostream& stream ) const {
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
      stream << "  Dx: " ;
      stream.width( 6 );
      stream << BCvalue[ 0 ];
    }
    if ( BCstatus[ 1 ] ) {
      stream << "  Dy: " ;
      stream.width( 6 );
      stream << BCvalue[ 1 ];
    }
    if ( BCstatus[ 2 ] ) {
      stream << "  Dz: " ;
      stream.width( 6 );
      stream << BCvalue[ 2 ];
    }
    if ( BCstatus[ 3 ] ) {
      stream << "  Rx: " ;
      stream.width( 2 );
      stream << BCvalue[ 3 ];
    }
    if ( BCstatus[ 4 ] ) {
      stream << "  Ry: " ;
      stream.width( 6 );
      stream << BCvalue[ 4 ];
    }
    if ( BCstatus[ 5 ] ) {
      stream << "  Rz: " ;
      stream.width( 6 );
      stream << BCvalue[ 5 ];
    }
    if ( BCstatus[ 6 ] ) {
      stream << "  Fx: " ;
      stream.width( 6 );
      stream << BCvalue[ 6 ];
    }
    if ( BCstatus[ 7 ] ) {
      stream << "  Fy: " ;
      stream.width( 6 );
      stream << BCvalue[ 7 ];
    }
    if ( BCstatus[ 8 ] ) {
      stream << "  Fz: " ;
      stream.width( 6 );
      stream << BCvalue[ 8 ];
    }
    if ( BCstatus[ 9 ] ) {
      stream << "  Mx: " ;
      stream.width( 6 );
      stream << BCvalue[ 9 ];
    }
    if ( BCstatus[ 10 ] ) {
      stream << "  My: " ;
      stream.width( 6 );
      stream << BCvalue[ 10 ];
    }
    if ( BCstatus[ 11 ] ) {
      stream << "  Mz: " ;
      stream.width( 6 );
      stream << BCvalue[ 11 ];
    }
  }

}
