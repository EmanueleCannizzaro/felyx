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

#include "EfObject.h"

#include "PreProcessing.h"
#include "Lanczos.h"

using namespace felyx;



////
// Print global status of actual FEM evaluation
////
void EfObject::PrintGlobalStatus(){
  OutStream << "# Print Global Status: " << endl;
  OutStream << "\t - Model loaded from        : " << InterfacePtr->GetLoadPath() << endl; 
  OutStream << "\t - # Elements               : " << Elements.size() << endl;
  OutStream << "\t - # Nodes                  : " << Nodes.size() << endl;
  OutStream << "\t - # Materials              : " << Materials.size() << endl;
  OutStream << "\t - # Active DOF's           : " << DofCount << endl;
  OutStream << "\t - # Memory needs GSM [MB]  : " 
      << (int)(Profile * sizeof(float_type) / (1024*1024) ) << endl;
  OutStream << "\t - # Memory needs GMM [MB]  : " 
      << (int)(Profile * sizeof(float_type) / (1024*1024) ) << endl;
#ifdef USE_FLOATS
  OutStream << "\t - Precision of val type    : float" << endl;
#else
  OutStream << "\t - Precision of val type    : double" << endl;
#endif
}


////
// Function which calculates the Eigenfrequencies
////
int EfObject::CalcEf(unsigned nEf, bool Ev, std::vector<double> &EfRes){
  
  for (PtrVector<Material*>::iterator it=Materials.begin(); it!=Materials.end(); ++it) {
    if ((*it)->Get("rho")<=0) {
      cerr << "ERROR in CalcEf " << endl;
      cerr << "There are materials with no density specified. " << endl;
      exit(1);
    }
  }

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Link node DOF's to GSM/GMM index" << endl; 
  DofCount = LinkNodes2Gsm(Nodes, Elements);

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval envelope and profile of GSM : ";

  dense1D<unsigned> Envelope(DofCount);
  Profile = EvalEnvelope(Envelope, Elements);

  if (noise > 0)
    OutStream << Profile
	      << " -> Memory needs of GSM: "
	      << (int)( Profile * sizeof(float_type) / (1024*1024)) << " MB " << endl;

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval envelope and profile of GMM : ";

  Profile = EvalEnvelope(Envelope, Elements);

  if (noise > 0)
    OutStream << Profile
	      << " -> Memory needs of GMM: "
	      << (int)( Profile * sizeof(float_type) / (1024*1024)) << " MB " << endl;

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Initialize GSM of size : " << DofCount << endl;
  typedef mtl::matrix< float_type , symmetric<lower>, envelope<>, row_major >::type EnvelopeMatrix;
  EnvelopeMatrix GSM(Envelope, Profile );

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Initialize GMM of size : " << DofCount << endl;
  EnvelopeMatrix GMM(Envelope, Profile );

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Evaluate ESM's and assemble them to GSM" << endl;
  AssembleGM(Elements, GSM);

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Evaluate EMM's and assemble them to GMM" << endl;
  AssembleGmm(Elements, GMM);

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Compute Eigenfrequencies" << endl;
  int solved = mtl::lanczos(GSM,GMM,Envelope,nEf,Ev,EfRes);

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status : " << solved << endl;

  return solved;

}

/////
//  Eval total mass of active elements
////
float_type EfObject::EvalMass() const{
   
  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval total mass of active Finite Elements in the model. M = ";
  float_type mass = fe_base::EvalMass( Elements );
  if (noise > 0) OutStream << mass << endl;
  return mass;
}

