//-----------------------------------------------------------------------------
// Plane182.h
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
#ifndef Plane182_h
#define Plane182_h Plane182_h

#include "BaseSolid.h"

using namespace std;
	// Using namespace from mat_base

namespace fe_base{		// Put classes into namespace fe_base
  
  //! Plane182 	==> 4-Node-2D-Structural-Plane
  class Plane182:public BaseSolid{

  public:
    // Constructors
    // ------------
    Plane182()  	: BaseSolid(NodeCount) {}
    virtual ~Plane182()	{}
    //! Virtual clone function, needed in PtrVector
    virtual Plane182* Clone() const { return new Plane182(*this); }

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


    //! Evaluate the shape functions derived with respect to standard cooords xi and eta.
    /*! EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to 
      the standard coords xi and eta in a specific point.  Its rows correspond to a specific 
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
    const static int             NodeCount;
    const static int             Id;
    const static string          Name;  
    static StructDofSet          ElementDofSet;
    static Dense_Matrix IntPoints;

  }; //of Class Plane182
  
} // of namespace

#endif


