//-----------------------------------------------------------------------------
// ElementFactory.h
//
// begin     : Feb 20 2002
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

#ifndef ElementFactory_h
#define ElementFactory_h ElementFactory_h

#include <iostream>

// Include factory and singleton from loki library
#include "utils/LokiFactory.h"     
#include "utils/LokiSingleton.h"   

// Include base class for the factory
#include "PtrVector.h"
#include "StructElement.h"
#include "LcmElement.h"

namespace fe_base{

  //! Define custom error object for the ElementFactory
  template <typename IdentifierType, class AbstractProduct>
  struct ElementFactoryError
  {
    static AbstractProduct* OnUnknownType(IdentifierType){
      cerr << endl << "ERROR in class ElementFactoryError: " << endl;
      cerr <<         "Unknown Type requested from ElementFactory" << endl << endl;
      exit(1);
    }
  }; 

  //! Define the Singleton Factory
  /*  Singleton Factory is used as described in
      "Modern C++ Design from Andrei Alexandrescu, p197ff"
      and implemented in "Loki".
  */
  typedef Loki::SingletonHolder
  <
  Loki::Factory  < StructElement, int, StructElement* (*)(), ElementFactoryError >
  > 
  StructElementFactory;

  typedef Loki::SingletonHolder
  <
  Loki::Factory  < LcmElement, int, LcmElement* (*)(), ElementFactoryError >
  > 
  LcmElementFactory;

};
 

#endif
