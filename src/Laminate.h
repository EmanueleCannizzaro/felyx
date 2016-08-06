//-----------------------------------------------------------------------------
// Laminate.h
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
#ifndef FE_BASELAMINATE_H
#define FE_BASELAMINATE_H

#include <vector>
#include "Layer.h"

#include "mtl/mtl_felyx_utils.h"

namespace fe_base {

/**
Derived vector class to hold a laminate consisting of layers.
The class was named LayeredMaterial in an earlier life!
@author Oliver Koenig, Boris Meier, Marc Wintermantel, Nino Zehnder
*/
  class Laminate : public std::vector<Layer*>
  {
  public:
    Laminate()
      : std::vector<Layer*>()
    {}

    Laminate( const unsigned n )
      : std::vector<Layer*>(n, NULL)
    {}

    ~Laminate(){}

    std::string ClassName() const { return "Laminate"; }

    //! Add top layer to the laminate
    inline void AddTopLayer( Layer* aLayer ){ push_back( aLayer ); }

    inline Layer* GetLayerPtr( const unsigned l ) const { return operator[](l); }
    inline double GetLayerThickness( const unsigned l ) const{ return operator[](l)->GetThickness(); }
    inline MaterialOrientation GetLayerMaterialOrientation( const unsigned l ) const { return operator[](l)->GetMaterialOrientation(); };
    inline Material* GetLayerMaterialPtr( const unsigned l ) const { return operator[](l)->GetMaterialPtr(); }
    inline unsigned GetLayerNumber() const { return size(); };

    double GetLaminateThickness() const;

    mtl::Dense_Matrix GetAMatrix( std::vector<MaterialOrientation> ) const;
    mtl::Dense_Matrix GetABDMatrix( std::vector<MaterialOrientation> ) const;

    //! Returns the thickness coordinate of the top interface of layer number nl
    double GetTopInterface( const unsigned ln ) const ;

    //! Returns the thickness coordinate of the bottom interface of layer number nl
    double GetBottomInterface( const unsigned ln ) const;

    void printOn(std::ostream&) const;

    //output functions
    friend std::ostream& operator<<(std::ostream&, const Laminate&); // overloaded << operator
    void print(std::ostream&) const;
  };

};

#endif
