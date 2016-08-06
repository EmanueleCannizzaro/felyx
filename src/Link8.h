//-----------------------------------------------------------------------------
// Link8.h
//
// begin     : Dec 20 2001
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
   

#ifndef Link8_h
#define Link8_h Link8_h

#include "BaseBeam.h"

using namespace std;
	// Using namespace from mat_base

namespace fe_base{		// Put classes into namespace fe_base
    
//! Link8 ==> 2-Node structural Beam
/*! Implemented as described in:
  Methode der Finiten Elemente,
  H.R.Schwarz, Teubner Studienbuecher, p71f */
  class Link8:public BaseBeam{
  public:
    // Constructors
    
    // ------------
    //! Standard constructor for Link8 element
    Link8()		: BaseBeam(NodeCount) {}
    //! Standard destructor for Link8 element
    virtual ~Link8()	{}
    //! Virtual clone function, needed in PtrVector
    virtual Link8* Clone() const { return new Link8(*this); }
    
    // Get functions
    // -------------
    //! Returns the number of nodes of an element.
    virtual unsigned	GetNodeCount()	const { return NodeCount;}
    //! Returns the ID-number of the element.
    virtual unsigned    GetId()        	const { return Id; }
    //!Returns the name of the element.    
    virtual string	GetName()	const { return Name; }
    //! Returns a set which contains the degrees of freedom of the element.
    virtual StructDofSet    	GetDofSet()	const { return ElementDofSet; }
    //! Returns the size of the element stiffnes matrix.
    virtual unsigned 	GetEMSize()	const { return ElementDofSet.count()*NodeCount; }


    //! Eval element stiffness matrix
    virtual Dense_Matrix CalcEM();   

    //! Eval element mass matrix
    virtual Dense_Matrix CalcEmm();  
    
    //! Eval length of beam - non virtual function, makes no sense for most of the elements
    double EvalLength() const;
    
    //! Returns the volume of this element
    virtual double EvalVolume() const; 
    
    //! Eval 3x3 transformation matrix from global to local coordinates - non virtual function
    Dense_Matrix EvalTransformationMatrix();
    
    //! Eval stress vectors at nodes
    /*! For link elements, stresses are evaluated in local coordinates, x in beam direction! 
        The stress matrix of the element is filled as follows:
    For each node, 1 value is stored:
    - (0) Axial normal stress
          ( Only results in valid values for axisymmetric crossection!)
    */
    virtual void EvalStresses(std::string ="" );
    
    
    /*! Eval Euler Buckling of beam - the length_factor defines the kind of support.
        Typical factors are:  - 2     for one free (and loaded) end, the other end clamped
                              - 1     both ends simply supported
                              - 0.707 one end simply supported, the other end clamped
                              - 0.5   both ends clamped
        The result is a safety margin: mos = (sigma_euler / s_x) - 1
        Thus, buckling occurs for values below 0!
        Needless to say, the beam must be subjected to compressive loads. The result for
        beams subject to tensile stresses is always 0! (no buckling)
    */
    virtual double EvalEulerBuckling( double length_factor );

    // Static data members
    // -------------------
    const static int      NodeCount;
    const static int      Id;
    const static string   Name;
    static StructDofSet   ElementDofSet;

  };  //of class Link8  

} // of namespace

#endif


