//-----------------------------------------------------------------------------
// Solid185.h
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
#ifndef Solid185_h
#define Solid185_h Solid185_h

#include "BaseSolid.h"

using namespace std;
	// Using namespace from mat_base

namespace fe_base{		// Put classes into namespace fe_base
    
  
//! Solid185 	==> 8-Node-3D-Structural-Brick
  /*! Implemented as described in:
    Methode der Finiten Elemente,
    H.R.Schwarz, Teubner Studienbuecher */
  
  class Solid185:public BaseSolid{
    
  public:
    // Constructors
    // ------------
    //! Standard constructor for Solid185 element
    Solid185()		: BaseSolid(NodeCount) {}
    //! Standard destructor for Solid185 element
    virtual ~Solid185()	{}
    //! Virtual clone function, needed in PtrVector
    virtual Solid185* Clone() const { return new Solid185(*this); }
    
    // Get functions
    // -------------
    //! Returns the number of nodes of an element.
    virtual unsigned  GetNodeCount() 	const { return NodeCount;}
    //! Returns the ID-number of the element.
    virtual unsigned  GetId()        	const { return Id; }
    //!Returns the name of the element.    
    virtual string    GetName()		const { return Name; }
    //! Returns a set which contains the degrees of freedom of the element.
    virtual StructDofSet    GetDofSet()       const { return ElementDofSet; }
    //! Returns the size of the element stiffnes matrix.
    virtual unsigned  GetEMSize()      const { return ElementDofSet.count()*NodeCount; }
    //! Returns the matrix of integration points
    virtual Dense_Matrix GetIntPoints() const { return IntPoints; };

    //Other functions
    // --------------------------------
    //! Compute element stiffness matrix (ESM)
    /*! CalcEM() is a member function of the derived class. This gives more 
      flexibility to FE - Models with more than one elementtype. Inside CalcEM()
      the one can choose different integration schemes (full integration and the 
      b-Bar method). The functions used inside are implemented in the base-class
      and explained in more detail there. The matrix of integration points and 
      weights are stored in the base-class. The number of integration points has to be adjuste
    */
    virtual Dense_Matrix CalcEM();   // Compute element stiffness matrix
 
    //! Evaluate the shape functions derived with respect to standard cooords xi, eta, zeta.
    /*! EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to 
      the standard coords xi, eta, zeta in a specific point.  Its rows correspond to a specific 
      shape function. The cols to the direction of derivation
      The function has no return value and is therefore  called with a matrix for the shape functions
      which is beeing changed in the function and the evaluation point given as follows.
      EvalPoints is a Matrix of points,  where the shape functions are to be evaluated. 
      IP determines the specific point, means row of Evalpoints.
    */
    void EvalDerivedShapeFunc(Dense_Matrix, const Dense_Matrix , int) const;
    void EvalShapeFunc(Dense_Matrix, const Dense_Matrix , int);

    Dense_Matrix GetlNodeCoords();

   // Static data members
   // -------------------
    const static int    NodeCount;
    const static int    Id;
    const static string Name;
    static StructDofSet ElementDofSet;
    static Dense_Matrix IntPoints;

 };  //of class Solid185
} // of namespace

#endif


