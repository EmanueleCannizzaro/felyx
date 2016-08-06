//-----------------------------------------------------------------------------
// StructDofSet.cc
//
// begin     : Jan 14 2003
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
   

#include "StructDofSet.h"

using namespace fe_base;

// Special access, giving back pos, at which the bitcount-th bit (set bit) occurs
// If bitcount-th bit is not set at all, give back -1
int StructDofSet::GetDofPos( unsigned bitcount ){
  
  unsigned b=0;
  int i=-1;
  
  while (b < bitcount){
    i++;
    if ( i == int( size() ) )	// If bitcount-th bit that is set does not occur
      return -1;
    if ( test(5-i) )		// If bit i is set
	  b++;
  }
  return i;
}

// Count number of set bits until pos
size_t StructDofSet::count(unsigned pos){
  
  if (pos >= size() )
    return bitset<6>::count();
  
  size_t count=0;
  for(unsigned i=0; i <= pos; i++){
    if ( test(5-i) ) 
      count++;				
  }
  return count;
}

// Overload some extern operators
StructDofSet fe_base::operator& (const StructDofSet& ds1, const StructDofSet& ds2){
  StructDofSet result(ds1);
  result &= ds2;
  return result;
};

//!Check if Element is a 3D Solid
bool StructDofSet::IsElement3DSolid(){
  return ( test(3) && count()==3 );
}

//!Check if Element is a 3D Shell
bool StructDofSet::IsElement3DShell(){
  return ( count()==6 );
}

//Check the dimension of the element
unsigned StructDofSet::ElementDimension() const {
  int i = test(5) + test(4) + test(3);
  return i;
}
