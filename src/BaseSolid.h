//-----------------------------------------------------------------------------
// BaseSolid.h
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
   

#ifndef BaseSolid_h
#define BaseSolid_h BaseSolid_h

#include "NumIntElem.h"

using namespace std;


namespace fe_base{		// Put classes into namespace fe_base

  class BaseSolid:public NumIntElem{

  public:
    BaseSolid()            : NumIntElem(), MaterialPtr(NULL) {};
    BaseSolid(unsigned n_) : NumIntElem(n_),MaterialPtr(NULL) {};

    virtual ~BaseSolid() {};

    virtual BaseSolid* Clone() const = 0;

    //CLASS FUNCTIONS
    //--------------

    //! Set MaterialPtr
    virtual void SetMaterialPtr( Material* Ptr ) { MaterialPtr = Ptr; };
    virtual Material* GetMaterialPtr() const { return MaterialPtr; };    

    virtual double EvalMass() const;

    //VIRTUAL FUNCTIONS
    //-----------------

    //! Eval Element Stiffness Matrix (ESM)
    /*!
      Solid185 uses a special integration scheme (b-bar) and is 
      has also its own member function.
    */
    //-------------------------------------------------------
    virtual Dense_Matrix CalcEM();

    virtual Dense_Matrix CalcEmm();

    virtual Dense_Matrix GetIntPoints() const = 0;

    virtual Dense_Matrix GetlNodeCoords() = 0;
    
    //! Evaluates the Jacobimatrix in 2 or 3 dimensional space and for shells
    /*! Depending on the number translational dof in z direction the function
      returns a 2x2- or 3x3- Jacobimatrix for a specific point given by the IPth
      row of the matrix evalpoints. The shape function are beeing evaluated inside this
      function.
    */
    void EvalJacobi(Dense_Matrix, Dense_Matrix, Dense_Matrix, unsigned);

    //! Convertes standard coordinates to global coordinates
    /*! The Jacobian is transfered to the function. Its inversion is done within the function
    */
    virtual void stCoord2globalCoord(const Dense_Matrix , Dense_Matrix);

    virtual void SetBMatrix( const Dense_Matrix, Dense_Matrix, BMatrixTypes);
    virtual void SetNMatrix( const Dense_Matrix, Dense_Matrix);

    virtual void EvalStresses(std::string ="" );

    // Data members
    // -------------------
    Dense_Matrix IntPointStress;

  private:
    //DATA MEMBERS
    //------------
    //! Pointer to material  
    Material* MaterialPtr;
  };
} //end of namespace

#endif
