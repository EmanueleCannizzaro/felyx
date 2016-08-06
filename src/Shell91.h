//-----------------------------------------------------------------------------
// Shell91.h
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
   

#ifndef Shell91_h
#define Shell91_h Shell91_h

//#include "StructElement.h"

#include "LayeredShell.h"

using namespace std;

namespace fe_base{		// Put classes into namespace fe_base

  //! Shell91 	==> 8-Node-Structural-Shell
  /*! Implemented as described in:
    Analysis of thick and thin shell structures by curved finite elements,
    International Journal for Numerical Methods in Engineering, Vol.2, 419-451 (1970),
    Sohrabuddin Ahmad*/

  class Shell91:public LayeredShell{
    
  public:
    // Constructors
    // ------------
    //! Standard constructor for Shell91 element
    Shell91()             :LayeredShell(NodeCount){}
    //! Standard destructor for Shell91 element
    virtual ~Shell91()	{}
    //! Virtual clone function, needed in PtrVector
    virtual Shell91* Clone() const { return new Shell91(*this); } 
     
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
    //! Returns the matrix of integration points for a specific layer number
    virtual Dense_Matrix GetIntPoints(const unsigned layerNr) const;
 
    //Other functions
    // --------------------------------

    //! Compute element stiffness matrix (ESM), by calling a function from LayeredShell class
    virtual Dense_Matrix CalcEM()  { return CalcEMbyLayerSuperposition();  }

    //! Compute element mass matrix (EMM), by calling a function from LayeredShell class
    virtual Dense_Matrix CalcEmm() { return CalcEmmbyLayerSuperposition(); }
    

    //! Evaluate the shape functions derived with respect to standard coords xi, eta, zeta.
    /*! EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to 
      the standard coords xi, eta, zeta in a specific point.  Its rows correspond to a specific 
      shape function. The cols to the direction of derivation
      The function has no return value and is therefore  called with a matrix for the shape functions
      which is beeing changed in the function and the evaluation point given as follows.
      EvalPoints is a Matrix of points,  where the shape functions are to be evaluated. 
      IP determines the specific point, means row of Evalpoints.
    */
    virtual Dense_Matrix EvalDerivedShapeFunc( const Dense_Matrix::Row IntPoint ) const;
    
    
    Dense_Matrix GetlNodeCoords();

    //! Evaluate the derived deformations derived with respect to standard coords xi. eta, zeta.
    /*! EvalDerivedDeformations is a Matrix of the values of the derived deformations with respect to 
      the standard coords xi, eta, zeta in a specific point.  Its rows correspond to a specific 
      shape function. The cols to the direction of derivation
      The function has no return value and is therefore  called with a matrix for the shape functions
      which is beeing changed in the function and the evaluation point given as follows.
      EvalPoints is a Matrix of points,  where the shape functions are to be evaluated. 
      IP determines the specific point, means row of Evalpoints.
    */
   virtual void EvalDerivedDeformations( Dense_Matrix, const Dense_Matrix, const Dense_Matrix, const Dense_Matrix,unsigned );
   

   // Static data members
   // -------------------

    const static unsigned    NodeCount;
    const static unsigned    Id;
    const static string      Name;
    static StructDofSet      ElementDofSet;  
    static Dense_Matrix      IntPoints; 
  
    

  }; //of class Shell91

} // of namespace

#endif


