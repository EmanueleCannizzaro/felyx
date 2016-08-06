/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// StructElement.h
//
// begin     : Jan 3 2003
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
   

#ifndef StructElement_h
#define StructElement_h StructElement_h

#include <iostream>
#include <vector>

#include "Element.h"
#include "StructNode.h"		
#include "Material.h"
#include "Properties.h"
#include "StructDofSet.h"

//mtl-includes
#include "mtl/mtl_felyx_utils.h"

#include <boost/timer.hpp>


using namespace std;
using namespace mtl;

namespace fe_base{

  // Abstract base class for all elements
  // No instances of "StructElement" can be created
  class StructElement:public Element<StructNode, StructDofSet>{
  
  public:
    //! Enum used to identify layers for evaluation in layered elements
    enum  LayerEnum {BottomLayer, MiddleLayer, TopLayer};
  
    //! Matrices to hold stress and strain values:
    /*! First Dimension : 
        - Stress/strain vectors at Gauss points or nodes
        Second Dimension: 
        - Either components of stress/strain vector: sx, sy, sz, sxy, syz, sxz
        - Or some other values: Example Beam4: s_ax, s_b_y_max, s_b_z_max, s_xphi_max 
	    (see documentation of EvalStress() member function of the different elements)
    */
    //Dense_Matrix Strains;
    Dense_Matrix Stresses;
    
    // Constructors
    // ------------
    StructElement()			: Element<StructNode, StructDofSet>() { }
    StructElement(unsigned n_)		: Element<StructNode, StructDofSet>(n_){ }

    virtual ~StructElement(){ }
    //! Virtual clone function, needed in PtrVector
    virtual StructElement* Clone() const =0;

    //! Sets the pointer to material or gives a warning
    virtual void SetMaterialPtr( Material* );
    //! Sets the pointer to properties or gives a warning
    virtual void SetPropertiesPtr( PropertySet* );
    //! Sets the pointer to elecoordsys or gives a warning
    virtual void SetEleCoordSysPtr( CoordSys* );
    //! Set NodeDofSets based on the attached elements
    void SetNodeDofSet() const;  
    //! Sets the pointer to Laminate or gives a warning
    virtual void SetLaminatePtr( Laminate* );

    virtual unsigned      GetNodeCount()  const = 0;
    virtual unsigned      GetId()         const = 0;
    virtual string        GetName()       const = 0;
    virtual StructDofSet  GetDofSet()     const = 0;
    virtual unsigned      GetEMSize()     const = 0;


    /*!These functions return a pointer to NULL if the 
      derived Element doesnt possess such a pointer 
    */
    virtual Material*         GetMaterialPtr() const;
    virtual PropertySet*      GetPropertiesPtr() const;
    virtual CoordSys*         GetEleCoordSysPtr() const;
    virtual Laminate*         GetLaminatePtr() const;

    // Some other functions
    //---------------------
    
    //!Gives back the sum of the first three DOF's
    /*!This gives a feedback on the dimension of the element*/
    unsigned GetDimension();

    //! Eval Element Mass Matrix (EMM)
    virtual Dense_Matrix CalcEmm() = 0;

    //! Eval volume of an element
    virtual double EvalVolume() const;

    //! Eval mass of an element
    virtual double EvalMass() const = 0;
    
    //! Eval Euler buckling
    virtual double EvalEulerBuckling( double length_factor );

    typedef mtl::matrix< float_type , symmetric<lower>, envelope<>, row_major >::type EnvelopeMatrix;
    void AssembleElement2Gmm( EnvelopeMatrix& );

    ////
    // POSTPROCESSING
    ////
    
    // STRAINS
    // -------
    //! Eval strains at int points
    virtual Dense_Matrix evalStrains();
    
    //! Eval strains at int points of layer LayerEnum
    virtual Dense_Matrix evalStrains(const LayerEnum Layer );
    
    //! Eval strains at int points
    virtual Dense_Matrix EvalStrainsAtIntPoints();
    
    //! Eval strains at int points of layer Layer
    virtual Dense_Matrix EvalStrainsAtIntPoints( const unsigned Layer );
    
    //! Eval strains at int points in element coordiante system
    virtual Dense_Matrix evalElementDirectionStrains();
    
    //! Eval strains at int points of layer Layer in element coordinate directions.
    virtual Dense_Matrix evalElementDirectionStrains( const unsigned Layer );
    
    //! Eval strains at int points in material direction coordinate system
    virtual Dense_Matrix evalMaterialDirectionStrains();
    
    //! Eval strains at int points of layer Layer in main material directions.
    virtual Dense_Matrix evalMaterialDirectionStrains( const unsigned Layer );
    
    //! Eval max vonMises strain in element, based on calculations at int points
    virtual float_type evalMaxVonMisesStrain();
    
    //! Eval max failure criteria of type ("MaximumStress", "Tsai-Wu", "Tsai-Hill", and "Hashin".
    virtual float_type evalMaxFailureCriteria( std::string type );
    
    // STRESSES
    // --------
    //! Eval stresses
    virtual void EvalStresses( std::string StressType ="" );
    
    //! Eval stresses at int points in global coordinates
    virtual Dense_Matrix EvalStressesAtIntPoints();
    
    //! Eval stresses at int points of layer Layer at integration points
    virtual Dense_Matrix EvalStressesAtIntPoints( const unsigned Layer );
    
    //! Eval stresses at int points in element coordinate system
    virtual Dense_Matrix evalElementDirectionStresses();
    
    //! Eval stresses at int points of layer Layer in element coordinate directions.
    virtual Dense_Matrix evalElementDirectionStresses( const unsigned Layer );
    
    //! Eval stresses at int points in material coordinate system
    virtual Dense_Matrix evalMaterialDirectionStresses();
    
    //! Eval stresses at int points of layer Layer in main material directions
    virtual Dense_Matrix evalMaterialDirectionStresses ( const unsigned Layer );
    

    
    // I/O
    // ---
    // Provide << operator for all derived classes
    // --> operator must be public
    // Make use of the "Virtual Friend Function Idiom"
    friend ostream& operator<<(ostream&, const StructElement&);  
    
  };//of class StructElement
  
} // of namespace

#endif //StructElement_h

