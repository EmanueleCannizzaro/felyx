//-----------------------------------------------------------------------------
// PreProcessing.h
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

#ifndef PreProcessing_h
#define PreProcessing_h PreProcessing_h

// include FELyX headers
#include "ElementHeaders.h"
#include "PtrVector.h"
#include <algorithm>

using namespace std;

namespace fe_base {

  //! Evaluation of envelope and profile of the GSM
  template <class Array, class element_type>
  unsigned EvalEnvelope( Array& , const PtrVector<element_type*>& );

  //! Link nodes to GSM (Global Stiffness Matrix) and give back its size.
  /*! Evaluate Node DOF vector and appropriate index for each node which
      links nodes to its appropriate positions in the GSM. 
      Evaluation is based on the DOFs of each node and the homogeneous boundary conditions. 
      Member function gives back the size of the GSM to be created. */
  template <class element_type, class node_type>
  unsigned LinkNodes2Gsm( vector<node_type>&, PtrVector<element_type*>& );

  //   template<class element_type, class node_type>
  //   unsigned LinkNodes2Gmm(vector<node_type>&, PtrVector<element_type*>& );

  //! Assembly: Evaluate element stiffnesses and fill GSM.
  template <class mtlMatrix, class element_type>
  void AssembleGM( PtrVector<element_type*>&, mtlMatrix& );

  //! Apply inhomogeonous boundary conditions:
  /*! Forces and displacement constraints not equal 0.
      This function is very slow, and should be replaced by a faster function as
      done for the sparse matrix formats (used with the pardiso sparse solver) 
      in function ApplyBoundCons 
  */
  template <class mtlMatrix, class Array, class node_type, class envelope>
  void ApplyLoads( vector<node_type>&, mtlMatrix&, Array&, envelope& );

  // Include implementation of the above functions
#include "PreProcessing.inl"

} // of namespace

#endif
