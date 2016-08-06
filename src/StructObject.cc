//-----------------------------------------------------------------------------
// StructObject.cc
//
// begin     : Jan 2 2003
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

#include "StructObject.h"

#include "PostProcessing.h"

using namespace felyx;


////
// Postprocess operations
////
void StructObject::EvalStresses(string StressType ){

  // Eval element stresses
  if (noise>0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Postprocess: Evaluate stresses" << endl;
  EvalElementStresses( Elements, StressType );
  
//   for (vector<StructNode>::iterator it=Nodes.begin(); it !=Nodes.end(); ++it )
//     cout << "sp " <<      (*it).GetStresses()[0] << ", " << (*it).GetStresses()[1] << ", " << (*it).GetStresses()[2] << ", " << (*it).GetStresses()[3] << ", " << (*it).GetStresses()[4] << ", " << (*it).GetStresses()[5] << endl;
 }


/////
//  Eval total mass of active elements
////
float_type StructObject::EvalMass() const{
   
  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval total mass of active Finite Elements in the model. M = ";
  float_type mass = fe_base::EvalMass( Elements );
  if (noise > 0) OutStream << mass << endl;
  return mass;
}

////
// Eval compliance of structure based on nodal loads and deformations
////
float_type StructObject::EvalCompliance() const{

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval compliance of structure based on nodal loads and deformations. C = ";
  float_type compliance = fe_base::EvalCompliance( Nodes );
  if (noise > 0) OutStream << compliance << endl;
  return compliance;
}

////
// Eval maximum vonMises strain
////
float_type StructObject::evalMaxVonMisesStrain() const{

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval maximum von Mises strain . e_vonMises_max = ";
  float_type e_max = fe_base::evalMaxVonMisesStrain(Elements);
  if (noise > 0) OutStream << e_max << endl;
  return e_max;
}
  
////
// Eval maximum deformation
//// 
float_type StructObject::EvalMaxDeformation() const{

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval maximum deformation . maxdef = ";
  float_type d_max = fe_base::GetMaximumDeformation(Nodes);
  if (noise > 0) OutStream << d_max <<endl;
  return d_max;
}

////
// Eval maximum failure criteria
////
float_type StructObject::EvalMaxFailureCriteria( std::string type ) const{
  
  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval maximum failure criteria of type " << type << " = ";
  float_type s_max = fe_base::EvalMaxFailureCriteria(type, Elements);
  if (noise > 0) OutStream << s_max << endl;
  return s_max;
}

////
// Print global status of actual FEM evaluation
////
void StructObject::PrintGlobalStatus(){
  OutStream << "# Print Global Status: " << endl;
  OutStream << "\t - Model loaded from        : " << InterfacePtr->GetLoadPath() << endl; 
  OutStream << "\t - # Elements               : " << Elements.size() << endl;
  OutStream << "\t - # Nodes                  : " << Nodes.size() << endl;
  OutStream << "\t - # Materials              : " << Materials.size() << endl;
  OutStream << "\t - # Active DOF's           : " << DofCount << endl;
  OutStream << "\t - # Memory needs GSM [MB]  : " 
      << (int)(Profile * sizeof(float_type) / (1024*1024) ) << endl;
#ifdef USE_FLOATS
  OutStream << "\t - Precision of val type    : float" << endl;
#else
  OutStream << "\t - Precision of val type    : double" << endl;
#endif
  OutStream << "\t - Maximum deformation      : " << GetMaximumDeformation(Nodes) << endl;
  OutStream << "\t - Maximum vonMises stress (nodal solution)  : " << EvalMaxVonMisesStress(Nodes) << endl;
}
