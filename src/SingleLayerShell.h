//-----------------------------------------------------------------------------
// SingleLayerShell.h
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
   

#ifndef SingleLayerShell_h
#define SingleLayerShell_h SingleLayerShell_h

#include "BaseShell.h"

using namespace std;


namespace fe_base{		// Put classes into namespace fe_base

  class SingleLayerShell:public BaseShell{

  public:
    SingleLayerShell()            : BaseShell(), MaterialPtr(NULL), PropertiesPtr(NULL) {};
    SingleLayerShell(unsigned n_) : BaseShell(n_), MaterialPtr(NULL), PropertiesPtr(NULL) {};

    virtual ~SingleLayerShell() {};

    virtual SingleLayerShell* Clone() const = 0;

    //CLASS FUNCTIONS
    //--------------
    //! Set MaterialPtr
    virtual void SetMaterialPtr( Material* Ptr ) { MaterialPtr = Ptr; };    

    /*! Set PropertiesPtr, only if element allows propertysets, otherwise
      PropertiesPtr is set to NULL --> therefore a virtual function */
    virtual void SetPropertiesPtr( PropertySet* Ptr ) {PropertiesPtr = Ptr;} ;

    virtual double EvalMass() const;

    virtual double GetShellThickness() const {return PropertiesPtr->GetDouble("Thickness"); };
    virtual double GetLayerAngle( unsigned = 0 );
    virtual unsigned GetLayerCount() const { return 1; };
    virtual PropertySet* GetPropertiesPtr() const { return PropertiesPtr; };


    Material* GetMaterialPtr() const { return MaterialPtr; };

    // POSTPROCESSING
    // --------------
    
    //! Eval strains at int points
    virtual Dense_Matrix evalStrains();
    virtual Dense_Matrix EvalStrainsAtIntPoints();
    
    //! Eval strains at int points in element coordiante system
    virtual Dense_Matrix evalElementDirectionStrains();
    
    //! Eval strains at int points in material direction coordinate system
    virtual Dense_Matrix evalMaterialDirectionStrains();
    
    //! Eval max vonMises strain in element, based on calculations at int points
    virtual float_type evalMaxVonMisesStrain();
    
    //! Eval max failure criteria of type  "MaximumStress", "Tsai-Wu", "Tsai-Hill", and "Hashin"
    virtual float_type evalMaxFailureCriteria( std::string type );
    
    //! Eval stresses at int points in global coordinates
    virtual Dense_Matrix EvalStressesAtIntPoints();
    
    //! Eval stresses at int points in element coordinate system
    virtual Dense_Matrix evalElementDirectionStresses();
    
    //! Eval stresses at int points in material coordinate system
    virtual Dense_Matrix evalMaterialDirectionStresses();
    
    
                       
  private:
    //DATA MEMBERS
    //------------
    //! Pointer to material  
    Material* MaterialPtr;
    //! Pointer to properties
    PropertySet* PropertiesPtr;

  };
} //end of namespace

#endif
