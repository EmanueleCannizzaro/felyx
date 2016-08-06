//-----------------------------------------------------------------------------
// Darcy3D8.h
//
// begin     : Mar 2 2005
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
#ifndef Darcy3D8_h
#define Darcy3D8_h Darcy3D8_h

#include "LcmElement.h"

using namespace std;
	// Using namespace from mat_base

namespace fe_base{		// Put classes into namespace fe_base
    
//! Darcy3D8 	==> 8-Node-3D-Thermal-Brick
  class Darcy3D8:public LcmElement{
    
  public:
    // Constructors
    // ------------
    //! Standard constructor for Darcy3D8 element
    Darcy3D8()		: LcmElement(NodeCount) { }
    //! Standard destructor for Darcy3D8 element
    virtual ~Darcy3D8()	{}
    //! Virtual clone function, needed in PtrVector
    virtual Darcy3D8* Clone() const { return new Darcy3D8(*this); }
    
    // Get functions
    // -------------
    //! Returns the number of nodes of an element.
    virtual unsigned  GetNodeCount() 	const { return NodeCount;}
    //! Returns the ID-number of the element.
    virtual unsigned  GetId()        	const { return Id; }
    //!Returns the name of the element.
    virtual string    GetName()		const { return Name; }
    //! Returns a set which contains the degrees of freedom of the element.
    virtual LcmDofSet    GetDofSet()       const { return ElementDofSet; }
    //! Returns the size of the element stiffnes matrix.
    virtual unsigned  GetEMSize()      const { return ElementDofSet.count()*NodeCount; }


    //Other functions
    // --------------------------------

   //! Evaluate the shape functions derived with respect to standard cooords xi, eta, zeta.
    /*! EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to 
      the standard coords xi, eta, zeta in a specific point.  Its rows correspond to a specific 
      shape function. The cols to the direction of derivation
      The function has no return value and is therefore  called with a matrix for the shape functions
      which is beeing changed in the function and the evaluation point given as follows.
      EvalPoints is a Matrix of points,  where the shape functions are to be evaluated. 
      IP determines the specific point, means row of Evalpoints. The shape functions
      are given in volume coordinates L1,L2,L3 and L4.
    */
    static Dense_Matrix GetIntegrationPoints();
    void EvalDerivedShapeFunc(Dense_Matrix&, Dense_Vector&);

    //   Dense_Matrix GetlNodeCoords();
    Dense_Matrix CalcEM();
    void EvalJacobi(Dense_Matrix&, Dense_Matrix&);
    void stCoord2globalCoord(const Dense_Matrix&,  Dense_Matrix& );
    void Local2Global(Dense_Matrix&);

    //! Compute element velocity
    virtual void CompEleVelocity();
    void CompPressGrad(Dense_Vector&, Dense_Vector&);

    //! Compute flux between nodes of element
    virtual void Flux(Dense_Matrix&, unsigned, unsigned, unsigned, unsigned);

    virtual void EvalCV();

    virtual void CompVolFlux();

    //! Overwritten function
    virtual void EvalCentroid();

    // Static data members
    // -------------------
    const static int    NodeCount;
    const static int    Id;
    const static string Name;
    static LcmDofSet ElementDofSet;
    static Dense_Matrix IP;

  private:

    Dense_Vector scs, scsd, scsab; // interface nodes a-b
    Dense_Vector sds, sdsb, sdsac; // interface nodes a-c
    Dense_Vector sbs, sbsc, sbsad; // interface nodes a-d
    Dense_Vector sas, sasd, sasbc; // interface nodes b-c
    Dense_Vector      scsa, scsbd; // interface nodes b-d
    Dense_Vector      sasb, sascd; // interface nodes c-d

    Dense_Vector      ssab, ssac, ssad, ssbc, ssbd, sscd; // interface nodes c-d
    Dense_Vector      ssa, ssb, ssc, ssd;

    Dense_Matrix K_Mat;            // Permeability Matrix
                                   // Kxx Kxy Kxz
                                   // Kyx Kyy Kyz
                                   // Kzx Kzy Kzz

  };  //of class Darcy3D8
  
} // of namespace

#endif


