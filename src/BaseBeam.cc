//-----------------------------------------------------------------------------
// BaseBeam.cc
//
// begin     : Nov 20 2001
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


#include "BaseBeam.h"


using namespace fe_base;

////
// Eval mass of an element, based on volume calculation
////
double BaseBeam::EvalMass() const{
  FELYX_RUNTIME_ASSERT( MaterialPtr != NULL,"BaseBeam::EvalMass()");
  return EvalVolume() * MaterialPtr->Get("rho");
}

double BaseBeam::EvalEulerBuckling( double length_factor ){
  FELYX_RUNTIME_THROW("BaseBeam::EvalEulerBuckling: not defined yet");
}
