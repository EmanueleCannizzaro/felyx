//-----------------------------------------------------------------------------
// StructElement.cc
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

#include "StructElement.h"

using namespace fe_base;



////
// Eval Volume of an element - dummy function, cause function is not implemented in
// every element yet!
////
double StructElement::EvalVolume() const{
  std::cerr << GetName() << "::EvalVolume() not implemented! Doing nothing, returning 0.0."<< std::endl;
  return 0.0;
}

////
// Eval Euler Buckling of an element - dummy function, cause function is not implemented in
// every element yet!
////
double StructElement::EvalEulerBuckling( double l_fac){
  std::cerr << GetName() << "::EvalEulerBuckling() not implemented! Doing nothing, returning 0.0."<< std::endl;
  return 0.0;
}

Dense_Matrix StructElement::evalStrains(){
  std::cerr << GetName() << "::evalStrains() not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}

Dense_Matrix StructElement::evalStrains( const LayerEnum Layer ){
  std::cerr << GetName() << "::evalStrains(LayerEnum Layer) not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}


Dense_Matrix StructElement::EvalStrainsAtIntPoints(){
  std::cerr << GetName() << "::EvalStrainsAtIntPoints() not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    

Dense_Matrix StructElement::EvalStrainsAtIntPoints( const unsigned Layer ){
  std::cerr << GetName() << "::EvalStrainsAtIntPoints(Layer) not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    

Dense_Matrix StructElement::evalElementDirectionStrains(){
  std::cerr << GetName() << "::evalElementDirectionStrains() not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    

Dense_Matrix StructElement::evalElementDirectionStrains( const unsigned Layer ){
  std::cerr << GetName() << "::evalElementDirectionStrains(Layer) not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    

Dense_Matrix StructElement::evalMaterialDirectionStrains(){
  std::cerr << GetName() << "::evalMaterialDirectionStrains() not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    

Dense_Matrix StructElement::evalMaterialDirectionStrains( const unsigned Layer ){
  std::cerr << GetName() << "::evalMaterialDirectionStrains(Layer) not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}

float_type StructElement::evalMaxVonMisesStrain(){
  std::cerr << GetName() << "::evalMaxVonMisesStrain() not implemented! Doing nothing."<< std::endl;
  return 0.0;
}

float_type StructElement::evalMaxFailureCriteria( std::string type ){
  std::cerr << GetName() << "::evalMaxFailureCriteria(string) not implemented! Doing nothing."<< std::endl;
  return 0.0;
}
////
// Eval stresses of an element not implemented for every element yet!
////
void StructElement::EvalStresses( std::string StressType){
  std::cerr << GetName() << "::EvalStresses() not implemented! Doing nothing."<< std::endl;
}

Dense_Matrix StructElement::EvalStressesAtIntPoints(){
  std::cerr << GetName() << "::EvalStressesAtIntPoints() not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    
Dense_Matrix StructElement::EvalStressesAtIntPoints( const unsigned Layer ){
  std::cerr << GetName() << "::EvalStressesAtIntPoints(Layer) not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    
Dense_Matrix StructElement::evalElementDirectionStresses(){
  std::cerr << GetName() << "::evalElementDirectionStresses() not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    
Dense_Matrix StructElement::evalElementDirectionStresses( const unsigned Layer ){
  std::cerr << GetName() << "::evalElementDirectionStresses(Layer) not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    
Dense_Matrix StructElement::evalMaterialDirectionStresses(){
  std::cerr << GetName() << "::evalMaterialDirectionStresses() not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}
    
Dense_Matrix StructElement::evalMaterialDirectionStresses ( const unsigned Layer ){
  std::cerr << GetName() << "::evalMaterialDirectionStresses(Layer) not implemented! Doing nothing."<< std::endl;
  return Dense_Matrix();
}

////
// Set the material pointer of an element, only for the elements containing one...
////
void StructElement::SetMaterialPtr( Material* dummy_ ){

}

////
// Set the properties pointer of an element, only for the elements containing one...
////
void StructElement::SetPropertiesPtr( PropertySet* dummy_ ){
  std::string error = "Element type " + GetName() +" has no PropertiesPtr. SetPropertiesPtr() fails";
  FELYX_RUNTIME_THROW( error.c_str() );
}

////
// Set the EleCoordSys pointer of an element, only for the elements containing one...
////
void StructElement::SetEleCoordSysPtr( CoordSys* dummy_ ){
  std::string error = "Element type " + GetName() +" has no EleCoordSysPtr. SetEleCoordSysPtr() fails";
  FELYX_RUNTIME_THROW( error.c_str() );
}


////
// Set the laminate pointer of an element, only for the elements containing one...
////
void StructElement::SetLaminatePtr( Laminate* dummy_ ){
  std::string error = "Element type " + GetName() +" has no LaminatePtr. SetLaminatePtr() fails";
  FELYX_RUNTIME_THROW( error.c_str() );
}

unsigned StructElement::GetDimension(){
  StructDofSet myDofSet_ = GetDofSet();
  unsigned Dim = myDofSet_[5]+myDofSet_[4]+myDofSet_[3];
  return Dim;
}


Material* StructElement::GetMaterialPtr() const {
  return NULL;
}


PropertySet* StructElement::GetPropertiesPtr() const {
  return NULL;
}


CoordSys* StructElement::GetEleCoordSysPtr() const {
  return NULL;
}

Laminate* StructElement::GetLaminatePtr() const {
  return NULL;
}


// /*! Compare the existing NodeDofVec of every node with the DofVec
//   write the more free one back to NodeDofVec */
void StructElement::SetNodeDofSet() const{
  for (unsigned i = 0 ; i < NodeVec.size() ; ++i){
    NodeVec[i]->NodeDofSet |= GetDofSet();     //bitwise or
  }
}



////
// Implementation of Assembly of one element into global structure
////
void StructElement::AssembleElement2Gmm( EnvelopeMatrix& globMat ){

  // Evaluate mass matrix of element
  Dense_Matrix M = CalcEmm();

 // Transform EM, if there are any rotated nodal coordinate systems
  TransformEM(M);

  unsigned i, j;
  int globMatRow, globMatCol, nodenr, dofnr, dofpos;
  bool swapping = 0;

  for ( i=0; i < M.nrows(); ++i ){                      // loop over lower triangle of symmetric M

    //    cout << endl << " i = " << i << endl;
    nodenr = i / GetDofSet().count();                   // nodenr in this element
    //    cout << " nodenr = " << nodenr << endl;
    dofnr  = i % GetDofSet().count() + 1;               // dofnr-th DOF is looked at right now
    //    cout << " dofnr = " << dofnr << endl;
    dofpos = GetDofSet().GetDofPos( dofnr );            // this corresponds to dofpos  (x,y,z,rx,ry,rz)
    //    cout << " dofpos = " << dofpos << endl;
    globMatRow = NodeVec[nodenr]->GetGMIndex(dofpos);   // Get appropriate globMat row
    //    cout << " globMatRow = " << globMatRow << endl;

    if (globMatRow >= 0 ){                                      // If considered DOF is active (else globMatRow = -1)

      for ( j = 0; j <= i ; ++j ){                      // loop over lower triangle of symmetric K

        nodenr = j / GetDofSet().count();               // nodenr in this element
        dofnr  = j % GetDofSet().count() + 1;           // dofnr-th DOF is looked at right now
        dofpos = GetDofSet().GetDofPos( dofnr );        // this corresponds to dofpos  (x,y,z,rx,ry,rz)
        globMatCol = NodeVec[nodenr]->GetGMIndex(dofpos);       // Get appropriate globMat col

        if (globMatCol >= 0 ){                          // If considered DOF is active (else globMatCol = -1)

          // if considered dof is located in upper half of globMat, swap indices
          if ( globMatCol > globMatRow ){
            swapping = 1;
            swap( globMatRow, globMatCol );
          }

          globMat( globMatRow, globMatCol ) += M( i,j );        // Add M(i,j) to globMat

          // Swap indices back for further use
          if (swapping){
            swap( globMatRow, globMatCol );
            swapping = 0;
          }
        }
      }

    }
  }
}



////
//// I/O
////

// Overloading << operator for all derived classes
// using the "Virtual Friend Function Idiom"
ostream& fe_base::operator<< (ostream& stream, const StructElement& e)
{
  // Print element and all its associated entities
  stream.setf(ios::right, ios::adjustfield);

  stream << "Id          : " << e.GetId();
  stream << " ( Status : " << e.GetStatus() << " ) " << endl;
  if ( !e.NodeVec.empty()){
    stream << "Nodes       : " << endl;
    for (unsigned i=0; i < e.NodeVec.size(); ++i){
      stream.width(10); stream << i << " -- ";
      stream << *(e.NodeVec[i]) << endl;
    }
  }
  return stream;
}
