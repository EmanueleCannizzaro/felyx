//-----------------------------------------------------------------------------
// Layer.h
//
// begin     : July 2004
// copyright : (c) 2004 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
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
#ifndef FE_BASELAYER_H
#define FE_BASELAYER_H

#include "Material.h"

namespace fe_base {

  /**
  Class representing a single layer in a laminate for layered composite structures

  @author Oliver Koenig, Boris Meier, Marc Wintermantel, Nino Zehnder
  */
  class Layer{

  public:

    //Constructors
    Layer():MatPtr(NULL), thickness(0.0) { SetAngle(0.); }
    Layer( Material* ptr, double a, double t ) 
      : MatPtr(ptr), thickness(t)
    { SetAngle(a); }

    //Destructor
    ~Layer() {};

    //Get functions
    //-------------
    Material* GetMaterialPtr() const { return MatPtr; }
    double GetAngle() const; //returns the angle in degrees
    MaterialOrientation GetMaterialOrientation() const { return angle; };
    double GetThickness() const;
    std::string ClassName() const { return "Laminate"; }
    
    //Set functions
    //-------------
    void SetMaterial( Material* mat_ ) { MatPtr = mat_; }
    void SetAngle ( const double a_ );
    void SetThickness ( const double t_ ) { thickness = t_; }
  
    //! Print-function
   void printOn(ostream&) const;
   
   //!other functions
   bool operator==(const Layer&);
   bool operator!=(const Layer&);
    
  private:
    Material* MatPtr;
    MaterialOrientation angle;    //entered in degrees!!!
    double thickness;  
  };

};

#endif
