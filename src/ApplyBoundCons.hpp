/***************************************************************************
*   Copyright (C) 2005 by Oliver Koenig                                   *
*   Email: <koenig@even-ag.ch>                                            *
*                                                                         *
*   This file is part of FELyX (Finite Element Library eXperiment).       *
*   See library home page at http://felyx.sourceforge.net/                *
*                                                                         *
*   FELyX is free software; you can redistribute it and/or modify         *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   FELyX is distributed in the hope that it will be useful,              *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along FELyX; if not, write to the                                     *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/

#ifndef fe_base_ApplyBoundCons_h
#define fe_base_ApplyBoundCons_h

#include <algorithm>

namespace fe_base {

  //! Apply inhomogeneous boundary conditions
  /*! Within this function, the boundary conditions
      defined in BoundCon lists are applied to the global equation system Ax=B.
      Two different types of BC's are handled:
      - Given values for the right hand side vector B.
      - Prescribed values for the solution vector x.
      Compared to the existing solution in Preprocessing.inl, this solution si magnitudes faster!
  */
  template <class nodesT, class matrixT, class vectorT >
  void ApplyBoundCons( const nodesT& nodes, matrixT A, vectorT& B ) {

    // Init right hand side vector of equation system to zero
    std::fill( B.begin(), B.end(), 0.0 );

    // First loop through all boundary conditions
    // (For now this means to loop through all nodes, bad design...)
    
    Dense_Vector NodalLoads( 12, 0 );
    int gsmindex_int;   // This one is needed, since if GetGMIndex() fails, it returns -1 !!
    unsigned gsmindex;
    
    // Temporarly put prescribed BC's into a vector
    std::vector<float_type> prescribed_values( B.size(), 0.0 );
    
    for ( typename nodesT::const_iterator nit = nodes.begin(); nit != nodes.end(); ++nit ) {

      // If there is a BC set for this node
      if ( nit->ExistBoundCon() ) {

        // Get active load for that node
        nit->BoundConPtr->GetActiveLoads( NodalLoads );

        // Loop through DOF's of node
        for ( size_t i = 0; i < 6; ++i ) {

          // Get global index of this DOF
          gsmindex_int = nit->GetGMIndex( i );
          // Check if this DOF is active:
          if ( gsmindex_int >= 0 ) {

            gsmindex = ( unsigned ) gsmindex_int;

            // If there are, write RHS BC's to global RHS vector
            if ( NodalLoads[ i + 6 ] != 0 ) {
              B[ gsmindex ] += NodalLoads[ i + 6 ];
            }
            // If there are prescribed BC's different from zero
            if ( NodalLoads[ i ] != 0 ) {
              //Store prescribed prescribed_values temporarly in a vector
              //std::cout << "Prescribed displ: ind="<<gsmindex << " val="<<NodalLoads[i] <<"\n";
              prescribed_values[ gsmindex ] = NodalLoads[ i ];
            }

          }


        }
      }

    }

    // In a second loop we apply prescribed values directly by a single loop through
    // the global matrix
    typename matrixT::iterator cit = A.begin();
    typename matrixT::Column::iterator rit;
    typename matrixT::size_type cind, rind;
    while ( cit != A.end() ) {
      cind = cit.index();

      // Handle diagonal value
      // ---------------------
      rit = ( *cit ).begin();
      if ( prescribed_values[ cind ] != 0.0 ) {
        ( *rit ) = 1.0;
        //std::cout << "Diag: c="<<cind << " r="<<rit.index() << " val="<<prescribed_values[cind]<< "\n";
        B[ cind ] = prescribed_values[ cind ];
      }

      // Handle other values
      // -------------------
      ++rit;
      // Check if actual column has prescribed prescribed_values
      if ( prescribed_values[ cind ] != 0.0 ) {
        while ( rit != ( *cit ).end() ) {
          //std::cout << "colwise: c="<<cind << " r="<<rit.index() << " val="<<prescribed_values[cind]<< "\n";
          B[ rit.index() ] -= prescribed_values[ cind ] * ( *rit );
          ( *rit ) = 0.0;
          ++rit;
        }
      }
      // Else check if we have to apply prescribed displacement values for actual row
      else {

        while ( rit != ( *cit ).end() ) {
          rind = rit.index();
          if ( prescribed_values[ rind ] != 0.0 ) {
            //std::cout << "rowwise: c="<<cind << " r="<<rit.index() << " val="<<prescribed_values[rind]<< "\n";
            B[ cind ] -= prescribed_values[ rind ] * ( *rit );
            ( *rit ) = 0.0;
          }
          ++rit;
        }
      }

      ++cit;
    }

  }


}
#endif
