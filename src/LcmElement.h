/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// LcmElement.h
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
   

#ifndef LcmElement_h
#define LcmElement_h LcmElement_h

#include <iostream>
#include <vector>

#include "Element.h"
#include "LcmNode.h"		
#include "Material.h"
#include "Properties.h"
#include "LcmDofSet.h"

//mtl-includes
#include "mtl/mtl_felyx_utils.h"

#include <boost/timer.hpp>


using namespace std;
using namespace mtl;

namespace fe_base{		// Put classes into namespace fe_base

  // Abstract base class for all elements
  // No instances of "LcmElement" can be created
  class LcmElement:public Element<LcmNode, LcmDofSet>{
  public:


    // Constructors
    // ------------
    LcmElement()			: Element<LcmNode, LcmDofSet>(), MaterialPtr(NULL), sx(0), sy(0) {  }
    LcmElement(unsigned n_)		: Element<LcmNode, LcmDofSet>(n_), MaterialPtr(NULL), sx(0), sy(0) {  }

    virtual ~LcmElement(){ }
    //! Virtual clone function, needed in PtrVector
    virtual LcmElement* Clone() const =0;

  
    virtual unsigned      GetNodeCount()  const = 0;
    virtual unsigned      GetId()         const = 0;
    virtual string        GetName()       const = 0;
    virtual LcmDofSet     GetDofSet()     const = 0;
    virtual unsigned      GetEMSize()     const = 0;

    double Get(string) const;

    virtual PropertySet*      GetPropertiesPtr() const;
    virtual CoordSys*         GetEleCoordSysPtr() const;
    virtual Laminate*  GetLaminatePtr() const;

    void SetNodeDofSet() const;

    //! Check if element is active.
    bool GetStatus() const {
      bool b = false;
      for (unsigned i=0; i<GetNodeCount(); ++i) {
	if (!GetNodeIter(i)->ExistBoundCon()) b = true;
	else if (GetNodeIter(i)->ExistBoundCon()) if (GetNodeIter(i)->BoundConPtr->getValue(bc::P)!=0 || GetNodeIter(i)->BoundConPtr->getValue(bc::V)!=0) b = true;
      }
      return b; 
    }	
    
    virtual void EvalCentroid();

    unsigned GetDimension();

    virtual void SetMaterialPtr( Material* Ptr ) { MaterialPtr = Ptr; };
    virtual Material* GetMaterialPtr() const { return MaterialPtr; };    

    virtual void SetEleCoordSysPtr( CoordSys* Ptr ) { EleCoordSysPtr = Ptr; };
    virtual void CompEleVelocity() {};
    virtual void CompVolFlux() {};
    virtual void EvalCV() {};

    Material* MaterialPtr;
    CoordSys* EleCoordSysPtr;

    double sx, sy, sz;
    double vx, vy, vz;
    Dense_Vector a, b, c, s;
    Dense_Matrix T;
    
    void setWarp(double w) { warp = w; }
    void setWeft(double w) { weft = w; }

    double getWarp() { return warp; }    
    double getWeft() { return weft; }    
  
  private:

    //Two permeabilities ...
    double K1, K2, K3;
    double warp, weft;


    // Some other functions
    //---------------------
    
  };//of class LcmElement
  
} // of namespace

#endif //LcmElement_h

