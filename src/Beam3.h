//-----------------------------------------------------------------------------
// Beam3.h
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
   

#ifndef Beam3_h
#define Beam3_h Beam3_h

#include "BaseBeam.h"

using namespace std;

namespace fe_base{		// Put classes into namespace fe_base
    
//! Beam3 	==> 2-Node structural Beam
/*! Implemented as described in:
  Methode der Finiten Elemente,
  H.R.Schwarz, Teubner Studienbuecher, p72f 
  Stiffness Matrix implemented as described in:
  Concepts and Applications of Finite Element Analysis 
  R.D.Cook, D.S.Malkus, p26f */
  class Beam3:public BaseBeam{
  public:
    // Constructors
    // ------------
    //! Standard constructor for Beam3 element
    Beam3()		: BaseBeam(NodeCount) {}
    //! Standard destructor for Beam3 element    
    virtual ~Beam3()	{}
    //! Virtual clone function, needed in PtrVector
    virtual Beam3* Clone() const { return new Beam3(*this); }
    
    // Get functions
    // -------------
    //! Returns the number of nodes of an element.    
    virtual unsigned	GetNodeCount()	const { return NodeCount;}
    //! Returns the ID-number of the element.    
    virtual unsigned    GetId()        	const { return Id; }
     //!Returns the name of the element.    
    virtual string	GetName()	const { return Name; }
     //! Returns a set which contains the degrees of freedom of the element.
    virtual StructDofSet GetDofSet()	const { return ElementDofSet; }
    //! Returns the size of the element stiffnes matrix.
    virtual unsigned   	GetEMSize()	const { return ElementDofSet.count()*NodeCount; }

    //! Eval element stiffness matrix
    /*! Includes transverse shear deformation. If no shear deflection constant (ShearZ) is set,
	the displacement considers bending deformation only. 
	Shear deflection constants for several common sections are as follows: rectangle (6/5), solid circle (10/9),
	hollow (thin-walled) circle (2), hollow (thin-walled) square (12/5). Shear deflection constants for other
	cross-sections can be found in structural handbooks.*/
    virtual Dense_Matrix CalcEM();   
    //! Eval element mass matrix
    virtual Dense_Matrix CalcEmm();  

    // Constant static data members
    // -------------------
    const static int      NodeCount;
    const static int      Id ; 
    const static string   Name;
    static StructDofSet   ElementDofSet;
    
  };  //of class Beam3
  
} // of namespace

#endif


