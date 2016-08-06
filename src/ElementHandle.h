/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// Element.h
//
// begin     : Jun 17 2004
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
   

#ifndef ElementHandle_h
#define ElementHandle_h ElementHandle_h

#include <iostream>
#include <vector>

#include "Handle.h"
//#include "ElementHeaders.h"
//#include "LcmElement.h"
#include "Darcy2D3.h"
#include "Darcy3D3.h"
//#include "Link1.h"
#include "LcmNode.h"
#include "LcmDofSet.h"
//#include "FELyX.h"


using namespace std;
using namespace mtl;

namespace fe_base{		// Put classes into namespace fe_base

//  template<class node_type, class dof_type> class Element;
//  class LcmElement;
 // class Darcy2D3;
  
  class ElementHandle{
    
    public:
      ElementHandle() {}

      void create(unsigned id) {
         switch (id)
           {
             case 55 : 
              ep = new Darcy2D3; break;
             case 70 : 
              ep = new Darcy3D3; break; 
             default : 
              cout << "Unknown element type." << endl; break;
           }
         }
    
   //   unsigned GetNodeCount() {return ep->GetNodeCount();}

    private:
      Handle<Element<LcmNode, LcmDofSet> > ep;
      //Handle<> ep;
 //Handle<double>     ep;
  };

  
} // of namespace

#endif

