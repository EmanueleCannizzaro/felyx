//-----------------------------------------------------------------------------
// LayeredShell.h
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

#ifndef LayeredShell_h
#define LayeredShell_h LayeredShell_h

#include "BaseShell.h"

namespace fe_base{

  class LayeredShell:public BaseShell {

  public:
    LayeredShell()            : BaseShell(), LaminatePtr(NULL) {};
    LayeredShell(unsigned n_) : BaseShell(n_), LaminatePtr(NULL) {};

    virtual ~LayeredShell() {};

    virtual LayeredShell* Clone() const = 0;

    //CLASS FUNCTIONS
    //--------------

    //! Set MaterialPtr
    virtual void SetLaminatePtr( Laminate* Ptr ) { LaminatePtr = Ptr; };
    void AppendMaterialOrientation( MaterialOrientation a_ ) { LocalMaterialOrientationsVector.push_back(a_); };

    //!returns the pointer to the material of the specified layer
    Material* GetLayerMaterialPtr( unsigned nl ) const { return LaminatePtr->GetLayerMaterialPtr(nl ); };
    double GetLayerThickness( unsigned nl )      const { return LaminatePtr->GetLayerThickness( nl );  };
//    double GetLayerAngle( unsigned nl )          const { return LaminatePtr->GetLayerAngle( nl );      };

      public:
    /*!Function returning the orientation of the layer specified by n
       It checks wether the LocalMaterialOrientationsVector is filled
       if it is filled it returns the orientation stored there
       else it returns the angle beeing a member of class:layer
    */
    MaterialOrientation GetMaterialOrientationOfLayer(unsigned n);

    //! Mapping LayerEnum (bottom, middle, top) to layer number
    unsigned GetLayerNr(const LayerEnum layer) const;

  protected:
    //! Compute element stiffness matrix by superposition of the layer stiffness matrices.
    Dense_Matrix CalcEMbyLayerSuperposition();
    //! Compute element mass matrix by superposition of the layer stiffness matrices.
    Dense_Matrix CalcEmmbyLayerSuperposition();

    //VIRTUAL FUNCTIONS
    //-----------------
    virtual double GetShellThickness() const { return LaminatePtr->GetLaminateThickness(); };
    virtual Laminate* GetLaminatePtr() const { return LaminatePtr; };
    virtual unsigned GetLayerCount() const { return LaminatePtr->GetLayerNumber(); };

    //! Returns the matrix of integration points for a specific layer number
    /*! Must be implemented by all layered shell elements */
    virtual Dense_Matrix GetIntPoints(const unsigned layerNr ) const =0;

    // POSTPROCESSING
    // --------------

    //! Evaluate mass of layered shell
    virtual double EvalMass() const;

    //! Eval strains at int points of layer LayerEnum
    virtual Dense_Matrix evalStrains( const LayerEnum Layer );
    
    //! Eval strains at int points of layer Layer
    virtual Dense_Matrix EvalStrainsAtIntPoints( const unsigned Layer );
    
    //! Eval stresses at int points of layer Layer at integration points
    virtual Dense_Matrix EvalStressesAtIntPoints( const unsigned Layer );

    //! Eval max vonMises strain in element, based on calculations at int points of bottom and top layer
    virtual float_type evalMaxVonMisesStrain();
    
    //! Eval strains at int points of layer Layer in element coordinate directions.
    virtual Dense_Matrix evalElementDirectionStrains( const unsigned Layer );
    
    //! Eval strains at int points of layer Layer in main material directions.
    virtual Dense_Matrix evalMaterialDirectionStrains( const unsigned Layer );
    
    //! Eval stresses at int points of layer Layer in element coordinate directions.
    virtual Dense_Matrix evalElementDirectionStresses( const unsigned Layer );
    
    //! Eval stresses at int points of layer Layer in main material directions
    virtual Dense_Matrix evalMaterialDirectionStresses ( const unsigned Layer );
    
    //! Eval maximum stress criteria of layer Layer at integration points
    virtual Dense_Matrix evalMaxStressCriteria ( const unsigned layer );
    
    //! Eval Tsai-Hill criteria of layer Layer at integration points
    virtual Dense_Matrix evalTsaiHillCriteria ( const unsigned layer );
    
    //! Eval Tsai-Wu criteria of layer Layer at integration points
    virtual Dense_Matrix evalTsaiWuCriteria ( const unsigned layer );
    
    //! Eval Hashin criteria of layer Layer at integration points
    virtual Dense_Matrix evalHashinCriteria ( const unsigned layer );

    //! Eval max failure criteria of type "MaximumStress", "Tsai-Wu", "Tsai-Hill", and "Hashin"
    virtual float_type evalMaxFailureCriteria( std::string type );
    
    //DATA MEMBERS
    //------------
    //! Pointer to Laminate
    Laminate* LaminatePtr;
    std::vector<MaterialOrientation> LocalMaterialOrientationsVector;

  };
} //end of namespace

#endif
