//-----------------------------------------------------------------------------
// BaseBeam.h
//
// begin     : Nov 11 2002
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
   

#ifndef BaseBeam_h
#define BaseBeam_h BaseBeam_h

#include "StructElement.h"

using namespace std;


namespace fe_base{		// Put classes into namespace fe_base

  class BaseBeam:public StructElement{

  public:
    BaseBeam()            : StructElement(), MaterialPtr(NULL), PropertiesPtr(NULL) {};
    BaseBeam(unsigned n_) : StructElement(n_), MaterialPtr(NULL), PropertiesPtr(NULL) {};

    virtual ~BaseBeam() {};

    virtual BaseBeam* Clone() const = 0;

    //CLASS FUNCTIONS
    //--------------

    //! Set MaterialPtr
    virtual void SetMaterialPtr( Material* Ptr ) { MaterialPtr = Ptr; };    
    //! Get MaterialPtr
    virtual Material* GetMaterialPtr() const { return MaterialPtr; };    

    /*! Set PropertiesPtr, only if element allows propertysets, otherwise
      PropertiesPtr is set to NULL --> therefore a virtual function */
    virtual void SetPropertiesPtr( PropertySet* Ptr ) { PropertiesPtr = Ptr; };
    //! Return Properties Ptr
    virtual PropertySet* GetPropertiesPtr() const { return PropertiesPtr; };

    virtual double EvalMass() const;
    
    virtual double EvalEulerBuckling( double length_factor );

    //DATA MEMBERS
    //------------
    //! Pointer to material  
    Material* MaterialPtr;
    //! Pointer to properties
    PropertySet* PropertiesPtr;
   };
} //end of namespace

#endif
