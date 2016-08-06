//-----------------------------------------------------------------------------
// SingleLayerShell.cc
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
   

#include "SingleLayerShell.h"
#include "PostProcessing.h"

using namespace fe_base;

////
// Eval mass of an element
////
double SingleLayerShell::EvalMass() const{
  FELYX_RUNTIME_ASSERT(MaterialPtr != NULL,"SingleLayerShell::EvalMass()");
  return EvalVolume() * MaterialPtr->Get("rho");
}


double SingleLayerShell::GetLayerAngle( unsigned nl_ ){

  return PropertiesPtr->GetDouble("Theta");

}


Dense_Matrix SingleLayerShell::evalStrains(){
  return BaseShell::evalStrains( GetIntPoints() );
}

Dense_Matrix SingleLayerShell::EvalStrainsAtIntPoints(){
  return BaseShell::EvalStrainsAtIntPoints( GetIntPoints() );
}
   
Dense_Matrix SingleLayerShell::evalElementDirectionStrains(){
  return BaseShell::evalElementDirectionStrains( GetIntPoints() );
}

Dense_Matrix SingleLayerShell::evalMaterialDirectionStrains(){
  // get the angle in rads
  double PI = 4*atan(1.);
  float_type angle = GetLayerAngle() * PI / 180;
  return BaseShell::evalMaterialDirectionStrains( GetIntPoints(), angle );
}
 
float_type SingleLayerShell::evalMaxVonMisesStrain(){
  Dense_Matrix Strains = evalStrains();
  float_type e=0.0, maxval=0.0;
  for (unsigned i=0; i < Strains.nrows(); ++i){
    e = vonMises3D( mtl::rows( Strains )[i] );
    if (e > maxval)
      maxval = e;
  }
  return maxval;
}

float_type SingleLayerShell::evalMaxFailureCriteria( std::string type ){
  
  // get the int points
  Dense_Matrix IntPoints = GetIntPoints();  
  
  // get the angle in rads
  double PI = 4*atan(1.);
  float_type angle = GetLayerAngle() * PI / 180;
  
  // get the material
  Material* MatPtr = GetMaterialPtr();
  
  int id;
  if ( type == "MaximumStress")   id = 1;
  else if ( type == "Tsai-Hill")  id = 2;
  else if ( type == "Tsai-Wu")    id = 3;
  else if ( type == "Hashin")     id = 4;
  else std::cerr << "SingleLayerShell: Failure criteria type " << type << " not known!" << endl;
  
  switch(id) {
      case 1: {
        Dense_Matrix maxstress = evalMaxStressCriteria ( IntPoints, angle, MatPtr );
        double max(0.);
        // find the maximum value
        for (unsigned i = 0; i < maxstress.nrows(); ++i){
          for (unsigned j = 0; j < maxstress.ncols(); ++j){
            if (maxstress(i,j) > max) max = maxstress(i,j);
          }
        }
        return max;
      }
      
      case 2: {
        Dense_Matrix maxstress = evalTsaiHillCriteria ( IntPoints, angle, MatPtr );
        double max(0.);
        // find the maximum value
        for (unsigned i = 0; i < maxstress.nrows(); ++i){
          for (unsigned j = 0; j < maxstress.ncols(); ++j){
            if (maxstress(i,j) > max) max = maxstress(i,j);
          }
        }
        return max;
      }
      
      case 3: {
        Dense_Matrix maxstress = evalTsaiWuCriteria ( IntPoints, angle, MatPtr );
        double max(0.);
        // find the maximum value
        for (unsigned i = 0; i < maxstress.nrows(); ++i){
          for (unsigned j = 0; j < maxstress.ncols(); ++j){
            if (maxstress(i,j) > max) max = maxstress(i,j);
          }
        }
        return max;
      }
      
      case 4: {
        Dense_Matrix maxstress = evalHashinCriteria ( IntPoints, angle, MatPtr );
        double max(0.);
        // find the maximum value
        for (unsigned i = 0; i < maxstress.nrows(); ++i){
          for (unsigned j = 0; j < maxstress.ncols(); ++j){
            if (maxstress(i,j) > max) max = maxstress(i,j);
          }
        }
        return max;
      }   
  }
  return 0;
}

Dense_Matrix SingleLayerShell::EvalStressesAtIntPoints(){
  // get the angle in rads
  double PI = 4*atan(1.);
  float_type angle = GetLayerAngle() * PI / 180;
  return BaseShell::EvalStressesAtIntPoints( GetIntPoints(), GetMaterialPtr(), angle );
}

Dense_Matrix SingleLayerShell::evalElementDirectionStresses(){
  // get the angle in rads
  double PI = 4*atan(1.);
  float_type angle = GetLayerAngle() * PI / 180;
  return BaseShell::evalElementDirectionStresses( GetIntPoints(), GetMaterialPtr(), angle );
}

Dense_Matrix SingleLayerShell::evalMaterialDirectionStresses(){
  // get the angle in rads
  double PI = 4*atan(1.);
  float_type angle = GetLayerAngle() * PI / 180;
  return BaseShell::evalMaterialDirectionStresses( GetIntPoints(), GetMaterialPtr(), angle );
}


