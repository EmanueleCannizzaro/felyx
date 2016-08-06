/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// Element.h
//
// begin     : Jun 17 2004
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
   

#ifndef Handle_h
#define Handle_h Handle_h

//#include "Element.h"

using namespace std;
//8using namespace mtl;

namespace fe_base{		// Put classes into namespace fe_base

  template<class T>
  class Handle{
    
    public:
      Handle(): refptr(new size_t(1)), p(0) {}
      //Handle(const Handle& s): p(0) {if (s.p) p=s.p->Clone();}
      Handle(const Handle& s): refptr(s.refptr), p(s.p) {++*refptr;}
      Handle& operator=(const Handle&);
      ~Handle();
      
      Handle(T* t): refptr(new size_t(1)), p(t) {}
      
      operator bool() const {return p;}
      T& operator*() const;
      T* operator->() const;
    
      void make_unique() {
        if(*refptr != 1) {
          --*refptr;
          refptr = new size_t(1);
          p = p? p->Clone(): 0; }
      }
      
    private:
      size_t* refptr;
      T* p;
  };
  
  template<class T>
  Handle<T>::~Handle() {
    if(--*refptr == 0) {
      delete refptr;
      delete p; }
  }
  
  template<class T>
  Handle<T>& Handle<T>::operator=(const Handle& rhs) {
    ++*rhs.refptr;
    if (--*refptr == 0) {
      delete refptr;
      delete p; }
    refptr = rhs.refptr;
    p = rhs.p;
    return *this;
  }
  
  template<class T>
  T& Handle<T>::operator*() const {
    if(p)return *p;
    throw runtime_error("unbound Handle");
  }
  
  template<class T>
  T* Handle<T>::operator->() const {
    if(p)return p;
    throw runtime_error("unbound Handle");
  }
      
} // of namespace

#endif

