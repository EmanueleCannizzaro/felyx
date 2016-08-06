//-----------------------------------------------------------------------------
// LcmElement.cc
//
// begin     : Jan 3 2003
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
   
#include "LcmElement.h"

using namespace fe_base;


double LcmElement::Get(string myString_) const {

  if ( myString_ == "K1") 
    return K1;
  else if ( myString_ == "K2" )
    return K2;
  else {
    cerr << endl << "ERROR in LcmElement::Get(string) " 
	 << endl << "Wrong string asked " << endl <<  endl; 
    double dummy = 0.0;
    return dummy;
  } 

}



// /*! Compare the existing NodeDofVec of every node with the DofVec
//   write the more free one back to NodeDofVec */
void LcmElement::SetNodeDofSet() const{
  for (unsigned i = 0 ; i < GetNodeCount() ; i++){
    if (!NodeVec.empty())
      NodeVec[i]->NodeDofSet |= GetDofSet();     //bitwise or
  }  
}

void LcmElement::EvalCentroid(){
  unsigned n = GetNodeCount();
  sx=0; sy=0; sz=0;
  for (unsigned i=0; i<n; ++i) {
    sx += GetNodeIter(i)->Cx;
    sy += GetNodeIter(i)->Cy;
    sz += GetNodeIter(i)->Cz;
  }
  sx /= n; sy /= n; sz /= n;
}


unsigned LcmElement::GetDimension(){
  LcmDofSet myDofSet_ = GetDofSet();
  unsigned Dim = myDofSet_[5]+myDofSet_[4]+myDofSet_[3];
  return Dim;
}


// Material* LcmElement::GetMaterialPtr() const {
//   return NULL;
// }


PropertySet* LcmElement::GetPropertiesPtr() const {
  return NULL;
}


CoordSys* LcmElement::GetEleCoordSysPtr() const {
  return EleCoordSysPtr;
}

Laminate* LcmElement::GetLaminatePtr() const {
  return NULL;
}
