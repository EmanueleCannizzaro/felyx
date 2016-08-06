//-----------------------------------------------------------------------------
// StructDofSet.h
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
   

#ifndef StructDofSet_h
#define StructDofSet_h StructDofSet_h

#include <iostream>
#include <bitset>
#include <vector>

using namespace std;
namespace fe_base{		// Put classes into namespace fe_base

  //! Class StructDofSet is derived from bitset<6> to handle all dof stuff
  //! in nodes and elements
  class StructDofSet : public bitset<6>{
    
  public:

    // Constructors
    // ------------
    StructDofSet() : bitset<6>() {}			 // Empty constructor
    StructDofSet(const string str) : bitset<6>(str,0,6) {} // "string" constructor

    // Access to single bits
    // ---------------------
    // The following overloading did not work in gcc 2.95.x 
    // In the following three member functions, access is straight-forward
    // bit[0], or test(0) accesses first bit in set, not last as defined in bitset
    // bool test(size_t idx) const			{ return bitset<6>::test( size()-1-idx); }
    // reference operator[] (size_t idx){ return reference(*this, size()-1-idx );}
    // bool operator[] (size_t idx) const 		{ return test(idx); }

    // Special access, giving back pos, at which the bitcount-th bit (set bit) occurs
    // If bitcount-th bit is not set at all, give back -1
    int GetDofPos( unsigned );

    //!Returns true if the third bit is set (displacement in z direction) 
    bool IsElement3DSolid();
    //!Returns true if all bits are set to 1
    bool IsElement3DShell();
    //!Returns the count of the first 3 bits
    unsigned ElementDimension() const;
    
    // Overload "count()"
    // ------------------
    size_t count() { return bitset<6>::count(); }	// Count bits in DofVec
    size_t count(unsigned);			    	// Count bits in DofVec from begin to pos
    
  };

  // Overload operators not belonging to any class
  // ---------------------------------------------
  StructDofSet operator& (const StructDofSet&, const StructDofSet&);
  
  
} // of namespace

  
#endif
  
