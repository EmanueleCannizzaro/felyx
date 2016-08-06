/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// Beam4.h
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
   

#ifndef Beam4_h
#define Beam4_h Beam4_h

#include "BaseBeam.h"

using namespace std;

namespace fe_base{		// Put classes into namespace fe_base
    
//! Beam4 ==> 2-Node structural Beam
/*! Implemented as described in:
  Methode der Finiten Elemente,
  H.R.Schwarz, Teubner Studienbuecher, p72f
  CalEM implemented ad described in R.D.Cook, D.S.Malkus, p26f 
  Stiffness Matrix implemented as described in:
  Concepts and Applications of Finite Element Analysis 
  R.D.Cook, D.S.Malkus, p26f */

  class Beam4:public BaseBeam{
  public:

    // Constructors
    // ------------  
    //! Standard constructor for Beam4 element
    Beam4()		: BaseBeam(NodeCount) {}
    //! Standard destructor for Beam4 element
    virtual ~Beam4()	{}
    //! Virtual clone function, needed in PtrVector
    virtual Beam4* Clone() const { return new Beam4(*this); }
    
    // Get functions
    // -------------
    //! Returns the number of nodes of an element.
    virtual unsigned 	GetNodeCount()	const { return NodeCount;}
    //! Returns the ID-number of the element.
    virtual unsigned    GetId()        	const { return Id; }
    //!Returns the name of the element.    
    virtual string	GetName()	const { return Name; }
    //! Returns a set which contains the degrees of freedom of the element.
    virtual StructDofSet    	GetDofSet()	const { return ElementDofSet; }
    //! Returns the size of the element stiffnes matrix.
    virtual unsigned	GetEMSize()	const { return ElementDofSet.count()*NodeCount; }

    // Eval functions
    // --------------

    //! Eval element stiffness matrix
    /*! Includes transverse shear deformation. If no shear deflection constant (ShearZ and ShearY) is set,
    	the displacement considers bending deformation only. 
	Shear deflection constants for several common sections are as follows: rectangle (6/5), solid circle (10/9),
	hollow (thin-walled) circle (2), hollow (thin-walled) square (12/5). Shear deflection constants for other
	cross-sections can be found in structural handbooks.*/
    virtual Dense_Matrix CalcEM();   

    //! Eval element mass matrix
    virtual Dense_Matrix CalcEmm();   

    //! Returns the volume of this element
    virtual double EvalVolume() const;

    //! Eval stress vectors at nodes
    /*! For beam elements, stresses are evaluated in local coordinates, x in beam direction! 
        The stress matrix of the element is filled as follows:
	For each node, 4 values are stored:
	  - (0) Axial normal stress
	  - (1) Maximal bending stress in xy-plane: sigma_b_y(ThicknessY/2)
	  - (2) Maximal bending stress in xz-plane: sigma_b_z(ThicknessZ/2)
	  - (3) Maximal torsional shear stress    : tau_t(max( ThicknessY, ThicknessZ)/2 )
	        ( Only results in valid values for axisymmetric crossection!)
    */
    virtual void EvalStresses(std::string ="" );

    //! Eval 3x3 transformation matrix from global to local coordinates - non virtual function
    Dense_Matrix EvalTransformationMatrix();

    //! Eval length of beam - non virtual function, makes no sense for most of the elements
    double EvalLength() const;
    
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

  };  //of class Beam4
  
} // of namespace

#endif


