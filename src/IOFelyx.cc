//-----------------------------------------------------------------------------
// IOFelyx.cc
//
// begin     : Dec 1 2004
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

#include "IOFelyx.h"

void fe_base::create_element(PtrVector<StructElement*>::iterator eleit, int id_){
	*eleit = StructElementFactory::Instance().CreateObject( id_ ) ;
}


void fe_base::create_element(PtrVector<LcmElement*>::iterator eleit, int id_){
	*eleit = LcmElementFactory::Instance().CreateObject( id_ ) ; 
}

void fe_base::create_element(StructElement*& eleit, int id_){
	eleit = StructElementFactory::Instance().CreateObject( id_ ) ; 
}

void fe_base::create_element(LcmElement*& eleit, int id_){
	eleit = LcmElementFactory::Instance().CreateObject( id_ ) ; 
}
