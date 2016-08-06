//-----------------------------------------------------------------------------
// NumIntElem.h
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
   

#ifndef NumIntElem_h
#define NumIntElem_h NumIntElem_h

#include <stdexcept>
#include "StructElement.h"

using namespace std;


namespace fe_base{		// Put classes into namespace fe_base
  

  class NumIntElem:public StructElement{
    
    
  public:
  
    enum IntegrationScheme{AreaCoords3Points, VolumeCoords4Points, GaussQuadratic3D, 
			 GaussLinear3D, GaussQuadratic2D, GaussQubic2D,Gauss2x2x3,Gauss2x2x1,
			 Gauss3x3, Gauss3x2};
  
    enum BMatrixTypes{FullIntegration, ShellIntegration, bBarMethod};
  
    NumIntElem()            : StructElement()   {};
    NumIntElem(unsigned n_) : StructElement(n_) {};

    virtual ~NumIntElem() {};

    virtual NumIntElem* Clone() const = 0;

    //CLASS FUNCTIONS
    //---------------

    //! GetIntegrationPoints returns a matrix with points in standard coords and their weights
    /*! Get IntegrationPoints takes an enumerationtype  as input. The first two to four columns
      of the returned matrix are the coordinates of the integrationpoints in standard coordinates.
      These can be rectangular, voluminal or areal coordinates. The last column gives the weight 
      of the specific point for the integration. Since the integration is done in one loop only
      the weights are given as volume or area. The single integrationschemes are described in
      more detail in the implementation.
      This function needs to be static, because the Dense_Matrix IntPoints is static
    */
    static Dense_Matrix GetIntegrationPoints(IntegrationScheme);

    //! Three different integration schemes
    /*! The difference lies in the different assumptions for the solution of the differential
      equation. The bBarmethod evaluates the bBar therms inside the function.
      ShellIntegrationis used for shells...
    */     
    virtual void SetBMatrix( const Dense_Matrix, Dense_Matrix, BMatrixTypes) = 0;


    //VIRTUAL FUNCTIONS
    //-----------------

    //! virtual function to be defined in the derived classes
    virtual void stCoord2globalCoord (const Dense_Matrix, Dense_Matrix) = 0;

    //! Evaluates the shapefunctions in a specific point in the element
    /*! Virtual function for 2 and 3 dimensional elements implemented and described
      in the derived classes
    */
    // OLD FUNCTION, still needed by some elements
    virtual void EvalDerivedShapeFunc(Dense_Matrix, const Dense_Matrix, int) const{
      // The following exception should always be thrown, therefore no macro is used!
      throw std::logic_error("NumIntElem::EvalDerivedShapeFunc (OLD version): not implemented"); 
    }
    
    // NEW VERSION
    virtual Dense_Matrix EvalDerivedShapeFunc( const Dense_Matrix::Row IntPoint ) const{
      // The following exception should always be thrown, therefore no macro is used!
      throw std::logic_error("NumIntElem::EvalDerivedShapeFunc (NEW version): not implemented"); 
      return Dense_Matrix();
    }
    
    // Who needs this function ? - ok May 2004 
    virtual void EvalShapeFunc(Dense_Matrix, const Dense_Matrix, int);

    //!Returns the integration points stored in the derived elements
    virtual Dense_Matrix GetIntPoints() const = 0 ;

  };
} //end of namespace

#endif
