//-----------------------------------------------------------------------------
// StructAnalysisType.h
//
// begin     : Jan 12 2003
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
   

#ifndef StructAnalysisType_h
#define StructAnalysisType_h StructAnalysisType_h

#include "StructElement.h"
#include "StructBoundCon.h"
#include "StructNode.h"


namespace fe_base{	
    
  class StructAnalysisType{

  public:
    typedef StructElement element_type;
    typedef StructBoundCon bc_type;
    typedef StructNode node_type;
    typedef bc::StructBCtype bla;

  };

} // of namespace

#endif
