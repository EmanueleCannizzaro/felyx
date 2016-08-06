//-----------------------------------------------------------------------------
// PtrVector.h
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
   

#ifndef PtrVector_h
#define PtrVector_h PtrVector_h

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <typeinfo>

using namespace std;
namespace fe_base{		// Put classes into namespace fe_base

    ////
    //// DEFINITION OF TEMPLATE CLASS PTRVECTOR
    //// ======================================
    ////
    
    // Derive a class PtrVector from STL-vector to handle heterogeneous 
    // classes; ptrT must be a pointer to a base class
    template<class ptrT > 
    class PtrVector : public vector<ptrT>{
    public:
	// Constructors
	// ------------
	PtrVector() : vector<ptrT>() {}
	PtrVector( const PtrVector<ptrT>& );
	~PtrVector();
	
	// Adding elements
	// ---------------
	void push_back(const ptrT& );	// define push_back() new

    // Clear/Delete
    // ------------
    void clear();
    
	// Iterators
	// ---------
	typedef typename vector<ptrT>::iterator iterator;
	typedef typename vector<ptrT>::const_iterator const_iterator;
	
	iterator begin()  	{ return vector<ptrT>::begin(); }
	iterator end()  	{ return vector<ptrT>::end(); }
	const_iterator begin() const  { return vector<ptrT>::begin(); }
	const_iterator end()   const 	{ return vector<ptrT>::end(); }
    };

    ////
    //// IMPLEMENTATION OF TEMPLATE CLASS PTRVECTOR
    //// ==========================================
    ////
    
    // Copy constructor
    // use "virtual constructor" idiom -> member function "clone"
    // of class "ptrT" or a derived class of "ptrT"
    template<class ptrT>
    PtrVector<ptrT>::PtrVector(const PtrVector<ptrT>& Other)
	: vector<ptrT>( Other )
    {
	for (unsigned i=0; i < size(); i++)
	    (*this)[i] = Other[i]->Clone();
    }
    
    
    // Destructor which deletes every element once created with new
    template<class ptrT>
    PtrVector<ptrT>::~PtrVector()
    {
	// delete all stuff created with new
	for (unsigned i=0; i < size(); i++)
      if ( (*this)[i] != NULL )
	    delete (*this)[i];
    }

    // push_back function which simply creates new element for the new
    // pointer using "new"
    template<class ptrT>
    void PtrVector<ptrT>::push_back(const ptrT& X)
    {
	ptrT Y = X->Clone();	// clone the element the pointer points to
	insert(end(), Y);	// add it at the end (no memory realloc,
				// as long as capazity() > size()
    }
    
    // Overload clear function of std vector
    template<class ptrT>
    void PtrVector<ptrT>::clear(){
      // delete all stuff created with new
      for (unsigned i=0; i < size(); i++){
        if ( (*this)[i]!=NULL)
          delete (*this)[i];
      }
      
      std::vector<ptrT>::clear();
    }
  
    
    
} // of namespace

#endif
