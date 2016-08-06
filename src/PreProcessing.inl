//-----------------------------------------------------------------------------
// PreProcessing.inl
//
// begin     : Mai 8 2002
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

////
// Evaluation of envelope and profile of the GM
////
template <class Array, class element_type>
unsigned EvalEnvelope( Array& Envelope , const PtrVector<element_type*>& Elements ) {
  // initialize envelope vector to zero
  fill( Envelope.begin(), Envelope.end(), 0 );

  // loop through elements and call envelope member function of the element
  typename PtrVector<element_type*>::const_iterator eleit;
  for ( eleit = Elements.begin(); eleit < Elements.end(); ++eleit ) {
    if ( ( *eleit ) ->GetStatus() )                                	// check if element is active
      ( *eleit ) ->EvalEnvelope( Envelope );
  }

  // envelope of first row is per definition equal 1
  Envelope[ 0 ] = 1;

  // return sum of envelope vector, equals profile
  return accumulate( Envelope.begin(), Envelope.end(), 0 );
}

////
// Link nodes to global matrix and give back its size.
////
template <class element_type, class node_type>
unsigned LinkNodes2Gsm( vector<node_type>& Nodes, PtrVector<element_type*>& Elements ) {

  // Reset NodeDofSets
  typename vector<node_type>::iterator nodeit;
  for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ) {
    nodeit -> NodeDofSet.reset();
    nodeit -> Idx2GM = 0;
  }

  // Set node DOF vectors based on DOF's of elements
  typename PtrVector<element_type*>::iterator eleit;
  for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ) {
    if ( ( *eleit ) ->GetStatus() )                                	// check if element is active
      ( *eleit ) ->SetNodeDofSet();
  }

  // Set index for each node, pointing to first DOF in the GM
  // The active Dofs of each node are determined through GetActiveDof(),
  // which takes into account the NodeDofVec and the homogeneous BC's of a node
  unsigned Idx = 0;
  //DofSet ds;

  for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ) {
    nodeit -> Idx2GM = Idx;
    Idx += nodeit -> GetDofSet().count();
  }

  // Give back Size of GM to be created later on
  return Nodes.back().Idx2GM + Nodes.back().GetDofSet().count();
}

////
// Assembly: Evaluate element matrix and fill GM.
////
template <class Matrix, class element_type>
void AssembleGM( PtrVector<element_type*>& Elements, Matrix& globMat ) {

  typename PtrVector<element_type*>::iterator eleit;
  for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ) {
    if ( ( *eleit ) ->GetStatus() )                                 		// check if element is active
      ( *eleit ) -> AssembleElement2GM( globMat );	//Add the EM to the GM
  }
}

////
// Apply inhomogeonous boundary conditions:
////
template <class Matrix, class Array, class node_type, class envelope>
void ApplyLoads( vector<node_type>& Nodes, Matrix& globMat, Array& Forces, envelope& Envelope ) {

  // Set Force vector to zero
  fill( Forces.begin(), Forces.end(), 0.0 );

  Dense_Vector NodalLoads( 12, 0 );
  int gsmindex_int;   // This one is needed, since if GetGMIndex() fails, it returns -1 !!
  unsigned gsmindex;

  typename vector<node_type>::iterator nodeit;
  for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ) {

    if ( nodeit->ExistBoundCon() ) {			// If there is a BC set for this node:
      nodeit->BoundConPtr->GetActiveLoads( NodalLoads );	// Get active loads for that node

      for ( unsigned i = 0; i < 6; ++i ) {			// Loop over the 6 DOF's
        gsmindex_int = nodeit->GetGMIndex( i );		// Get global index of this DOF

        if ( gsmindex_int >= 0 ) {				// Check if this DOF is active:

          gsmindex = ( unsigned ) gsmindex_int;

          // If there are, write nodal force to global load vector
          if ( NodalLoads[ i + 6 ] != 0 )
            Forces[ gsmindex ] += NodalLoads[ i + 6 ];

          // If there is a displacement constraint different from zero, apply it
          if ( NodalLoads[ i ] != 0 ) {
            for ( unsigned j = 0; j < Forces.size(); ++j ) {	// loop over all DOF's
              if ( j < gsmindex && Envelope[ gsmindex ] > gsmindex - j )       // upper triangle
              {
                Forces[ j ] -= NodalLoads[ i ] * globMat( gsmindex, j );	// for all force values:
                globMat( gsmindex, j ) = 0;
              } else if ( Envelope[ j ] > j - gsmindex ) {              // lower triangle
                Forces[ j ] -= NodalLoads[ i ] * globMat( j, gsmindex );	// for all force values:
                // for actual row and column of GM
                globMat( j, gsmindex ) = 0;
              }
            }
            Forces[ gsmindex ] = NodalLoads[ i ];		// for the actual force value:
            globMat( gsmindex, gsmindex ) = 1;		// for globMat(gsmindex, gsmindex):
          }
        }
      }
    }
  }
}


