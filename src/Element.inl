// -*- c++ -*-
//-----------------------------------------------------------------------------
// Element.inl
//
// begin     : May 8 2002
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
//// Constructors
////

//! Copy constructor
template<class node_type, class dof_type>
Element<node_type, dof_type>::Element( const Element& e):
NodeVec(e.NodeVec), Active(e.Active)
{

}

////
//// Create and destroy elements
////

//! Function that allocates memory for n pointer to nodes
template<class node_type, class dof_type> 
void Element<node_type, dof_type>::create(unsigned n)
{
  NodeVec.resize(n);
  
  	// Set array of pointers to NULL
 	// for (unsigned i=0; i < n; i++)
  	//  NodePtr[i] = NULL;
}

//! Function that dealloctes the memory of an element reserved for the pointer to the nodes
template<class node_type, class dof_type> 
void Element<node_type, dof_type>::destroy()
{
  if(!NodeVec.empty()){
     NodeVec.clear();
  }
}

////
//// Set members
////

//! Set Node ptr of i-th node, returning true if range check for i succeeded.
template<class node_type, class dof_type> 
bool Element<node_type, dof_type>::SetNodeIter(unsigned i, typename std::vector<node_type>::iterator Iter){

  if ( i >= 0 && i < NodeVec.size() ){
    NodeVec[i] = Iter;
    return 1;
  }
  else
    return 0;
}

////
// Implementation of function transform element matrices, if element has any 
// nodes with rotated coordinate systems
////
template<class node_type, class dof_type> 
void Element<node_type, dof_type>::TransformEM( Dense_Matrix EM ){
  
  // Check if any transformations of EM have to be done
  // by finding out, if any of the nodes have nodal coord sys different from global coord sys
  vector<unsigned> RotatedNodes;

  for ( unsigned nodei=0; nodei < GetNodeCount(); ++nodei)
    if ( NodeVec[nodei]->ExistNodeCoordSys() )
      RotatedNodes.push_back(nodei);

  // Transform element matrix using nodal coord systems of rotated nodes found
  if ( RotatedNodes.size() > 0 ){
    
    const unsigned pr=0;
    if (pr) cout << endl << "START function Element::TransformEM( Dense_Matrix EM ): " << endl;
    unsigned nodei, dofi, dofj, posT, DofCount1, DofCount2;
    unsigned sizeEM = EM.nrows();
    
    // Create transformation matrix T, first as identity matrix
    Dense_Matrix T( sizeEM, sizeEM );
    mtl::set_value(T,0.0);
    mtl::set_diagonal(T,1.0);

    for ( unsigned i=0; i < RotatedNodes.size(); ++i ){
      nodei = RotatedNodes[i];

      // Check for not allowed coordinate system rotations for 2D models
      if ( GetDofSet().ElementDimension() < 3  && 
	   ( NodeVec[nodei]->NodeCoordSysPtr->Get("Thyz") != 0 ||  NodeVec[nodei]->NodeCoordSysPtr->Get("Thzx") != 0  ) ){
	cerr << "ERROR in TransformeM( Dense_Matrix EM ) " << endl;
	cerr << "This a 2D FE calculation: Thyz || Thzx != 0 not allowed! Exit program... " << endl;
	exit(1);
      }

      // Get rotation matrix of coord sys of this node
      Dense_Matrix T_node = NodeVec[nodei]->NodeCoordSysPtr->GetRotMat();
      
      // Inverse T_node
      Dense_Matrix T_node_inv( T_node.nrows(), T_node.ncols() );
      lu_inversion(T_node, T_node_inv );

      // Eval first index of T to be modified
      posT = nodei * GetDofSet().count();

      // Replace vals of T with vals from T_node_inv at appropriate places
      DofCount1 = GetDofSet().count();			// number of DOF's of this element
      DofCount2 = DofCount1 <= 3 ? DofCount1 : 3;	// number of DOF's for "first three DOF entries" of this element
     
      // Replace vals for first three DOF's of actual node
      for (dofi =0; dofi < DofCount2; ++dofi )
	for (dofj =0; dofj < DofCount2; ++dofj )
	  T( posT + dofi, posT+ dofj )= T_node_inv( dofi, dofj );

      // Replace vals for second three DOF's of actual node
      for ( dofi = 3; dofi < DofCount1; ++dofi )
	for ( dofj = 3; dofj < DofCount1; ++dofj )
	  T( posT + dofi, posT+ dofj )= T_node_inv( dofi-3, dofj-3 );

      if(pr){
	cout << " --> Node coord sys of node " << nodei <<" is used " <<  *( NodeVec[nodei]->NodeCoordSysPtr )  << endl;
	cout << "     - # DOF's of this element: " << DofCount1 << endl;
	cout << "     - Rotation matrix of node:" << endl;
	mtl::print_all_matrix( T_node_inv );
      }

    } 

    if(pr){
      cout << " --> Rotation matrix for this element: " << endl;
      mtl::print_all_matrix(T);
      cout << "END function Element::TransformEM( Dense_Matrix EM ): " << endl << endl;
    }

    Dense_Matrix temp(EM.nrows(),EM.ncols());
    mtl::copy(EM,temp);

    // Eval transformation
    mult_AT_B_A(T,temp,EM);
  
  }  

}

////
// Implementation of function to evaluate influence of an element on envelope and profile
////
template<class node_type, class dof_type>
template<class Array>
void Element<node_type, dof_type>::EvalEnvelope( Array& Envelope ) const {
  
  typename Array::value_type linebandwidth;
  
  unsigned nodei, nodej;
  int activedofi, activedofj,dof, gsmi, gsmj;
  bool swapping = 0;
  
  // Loop through stiffness matrix of element, on node level
  for ( nodei=0; nodei < GetNodeCount(); ++nodei){
    
    // number of active dofs of this node
    activedofi = NodeVec[nodei]->GetDofSet().count();
    
    // Do only something if this node has any active DOF's
    if ( activedofi > 0 ){
      
      // first gsm index of this node
      gsmi = NodeVec[nodei]->Idx2GM;
      
      // Loop through lower half stiffness matrix of element, on node level
      for ( nodej=0; nodej <= nodei; ++nodej ){
	
	// number of active dofs of this node
	activedofj = NodeVec[nodej]->GetDofSet().count();
	
	// Do only something if this node has any active DOF's
	if ( activedofj > 0 ){
	  
	  // first gsm index of this node
	  gsmj = NodeVec[nodej]->Idx2GM;
	  
	  // if considered node is in the upper half of symmetric GSM swap indices
	  if ( gsmj > gsmi ){
	    swapping = 1;
	    swap( gsmi, gsmj);
	    swap( activedofi, activedofj );
	  }
	  
	  // Check bandwidthes for all lines in GSM that get influenced from Node i
	  for (dof = 0; dof < activedofi; ++dof){
	      linebandwidth = gsmi - gsmj + 1 + dof;
	      if ( linebandwidth > Envelope[ gsmi + dof ] )
		Envelope[ gsmi + dof ] = linebandwidth;
	  }

	  // Swap indices back for further use
	  if (swapping){
	    swap( gsmi, gsmj);
	    swap( activedofi, activedofj );
	    swapping = 0;
	  }
	}
      }
    }
  }
}

////
// Implementation of Assembly of one element into global structure
////
template<class node_type, class dof_type>
template<class mtlMatrix>
void Element<node_type, dof_type>::AssembleElement2GM( mtlMatrix& globMat ){ 
 // Evaluate stiffness matrix of element
  Dense_Matrix K = CalcEM();

  // Transform EM, if there are any rotated nodal coordinate systems
  TransformEM(K);

  unsigned i, j;
  int globMatRow, globMatCol, nodenr, dofnr, dofpos;
  bool swapping = 0;

  for ( i=0; i < K.nrows(); ++i ){			// loop over lower triangle of symmetric K
    
    //    cout << endl << " i = " << i << endl;
    nodenr = i / GetDofSet().count();			// nodenr in this element
    //    cout << " nodenr = " << nodenr << endl;
    dofnr  = i % GetDofSet().count() + 1;		// dofnr-th DOF is looked at right now
    //    cout << " dofnr = " << dofnr << endl;
    dofpos = GetDofSet().GetDofPos( dofnr );		// this corresponds to dofpos  (x,y,z,rx,ry,rz)
    //    cout << " dofpos = " << dofpos << endl;
    globMatRow = NodeVec[nodenr]->GetGMIndex(dofpos);	// Get appropriate globMat row
    //    cout << " globMatRow = " << globMatRow << endl;

    if (globMatRow >= 0 ){					// If considered DOF is active (else globMatRow = -1)
      
      for ( j = 0; j <= i ; ++j ){			// loop over lower triangle of symmetric K
	
	nodenr = j / GetDofSet().count();		// nodenr in this element
	dofnr  = j % GetDofSet().count() + 1;		// dofnr-th DOF is looked at right now
	dofpos = GetDofSet().GetDofPos( dofnr );	// this corresponds to dofpos  (x,y,z,rx,ry,rz)
	globMatCol = NodeVec[nodenr]->GetGMIndex(dofpos);	// Get appropriate globMat col

	if (globMatCol >= 0 ){				// If considered DOF is active (else globMatCol = -1)
	  
	  // if considered dof is located in upper half of globMat, swap indices
	  if ( globMatCol > globMatRow ){
	    swapping = 1;
	    swap( globMatRow, globMatCol );
	  }
	  
	  globMat( globMatRow, globMatCol ) += K( i,j );	// Add K(i,j) to globMat

	  // Swap indices back for further use
	  if (swapping){
	    swap( globMatRow, globMatCol );
	    swapping = 0;
	  }
	}
      }
      
    }
  }
}

