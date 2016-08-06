//-----------------------------------------------------------------------------
// Darcy3D3.h
//
// begin     : Jun 16 2003
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
   

#ifndef Darcy3D3_h
#define Darcy3D3_h Darcy3D3_h

#include "LcmElement.h"

using namespace std;

namespace fe_base{		// Put classes into namespace fe_base
    
//! Darcy3D3 	==> 3-Node LCM membrane

  class Darcy3D3:public LcmElement{
  public:
    // Constructors
    // ------------
    //! Standard constructor for Darcy3D3 element
    //Darcy3D3()		: LcmElement(NodeCount) {}
    Darcy3D3()		: LcmElement(NodeCount), transformed(false), PropertiesPtr(NULL) {}
    //! Standard destructor for Darcy3D3 element    
    virtual ~Darcy3D3()	{}
    //! Virtual clone function, needed in PtrVector
    virtual Darcy3D3* Clone() const { return new Darcy3D3(*this); }
    
    bool transformed;

    // Get functions
    // -------------
    virtual void SetPropertiesPtr( PropertySet* Ptr ) {PropertiesPtr = Ptr;} ;
    virtual PropertySet* GetPropertiesPtr() const { return PropertiesPtr; };
    //! Returns the number of nodes of an element.    
    virtual unsigned	GetNodeCount()	const { return NodeCount;}
    //! Returns the ID-number of the element.    
    virtual unsigned    GetId()        	const { return Id; }
     //!Returns the name of the element.    
    virtual string	GetName()	const { return Name; }
     //! Returns a set which contains the degrees of freedom of the element.
    virtual LcmDofSet   GetDofSet()	const { return ElementDofSet; }
    //! Returns the size of the element stiffnes matrix.
    virtual unsigned   	GetEMSize()	const { return ElementDofSet.count()*NodeCount; }

    //! Eval element stiffness matrix
    virtual Dense_Matrix CalcEM();   

    //! Compute element velocity
    virtual void CompEleVelocity();

    //! Compute flux between nodes of element
    virtual void CompVolFlux();

    virtual void EvalCV();

    // Constant static data members
    // -------------------
    const static int      NodeCount;
    const static int      Id ; 
    const static string   Name;
    static LcmDofSet      ElementDofSet;

    PropertySet* PropertiesPtr;

  };  //of class Darcy3D3
  
} // of namespace

#endif
