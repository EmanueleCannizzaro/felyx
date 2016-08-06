//-----------------------------------------------------------------------------
// Layer.cc
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
#include <iomanip>

#include "Layer.h"

namespace fe_base {

    
void Layer::printOn(ostream& stream) const {
  stream << "Layer: MaterialPtr =  " << std::setw(10) << " Angle = " << std::setw(10) << GetAngle() << " Thickness = " << std::setw(10) << GetThickness();
}
bool Layer::operator==(const Layer& a)
{
	return ( thickness == a.thickness && angle == a.angle && MatPtr == a.MatPtr);
}

bool Layer::operator!=(const Layer& a)
{
	return !operator==(a);
}

double Layer::GetAngle() const
{
	double PI = 4*atan(1.);
	return angle.GetAngleMaterial1Direction()*180./PI;
}

void Layer::SetAngle ( const double a_ )
{
	double PI = 4*atan(1.);
	angle.SetAngleMaterial1Direction(a_*PI/180.);	
}

double Layer::GetThickness() const
{
    double h = thickness;
    
    return h;
}


};

