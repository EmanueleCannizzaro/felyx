//-----------------------------------------------------------------------------
// LayeredShell.cc
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


#include "LayeredShell.h"
#include "PostProcessing.h"

using namespace fe_base;

//! Compute element stiffness matrix by superposition of the layer stiffness matrices.
Dense_Matrix LayeredShell::CalcEMbyLayerSuperposition(){

  //The size of the square ESM is given as follows
  //Actually every node has only 5 DoF, which are translated to six global DoF
  unsigned  nSize = GetNodeCount() * GetDofSet().count(), ln;

  //Initialization of some Matrices
  Dense_Matrix FinalK(nSize, nSize);;
  mtl::set_value(FinalK,0.0);

  //This is set before the ESM calcualtion so that also for layered shells
  //only one set of Transformation Matrices is generated
//  SetMaterialTransformationVector();

  //Loop over all layers
  for ( ln = 0 ; ln < LaminatePtr->GetLayerNumber() ; ++ln ) {
    if ( LaminatePtr->GetLayerThickness(ln) != 0.0 ){
      mtl::add(SingleLayerESM( GetLayerMaterialPtr(ln), GetIntPoints(ln) , GetMaterialOrientationOfLayer(ln), GetLayerThickness(ln) ), FinalK );
    }
  }
  //****************************************
 //   std::cout << "Matrix: FinalK: " << std::endl;
 //   print_all_matrix(FinalK);
    //****************************************
  return FinalK;
}

//! Compute element mass matrix by superposition of the layer stiffness matrices.
Dense_Matrix LayeredShell::CalcEmmbyLayerSuperposition(){

  //The size of the square ESM is given as follows
  //Actually every node has only 5 DoF, which are translated to six global DoF
  unsigned  nSize = GetNodeCount() * GetDofSet().count(), ln;

  //Initialization of some Matrices
  Dense_Matrix FinalM(nSize, nSize);;
  mtl::set_value(FinalM,0.0);

  //Loop over all layers
  for ( ln = 0 ; ln < LaminatePtr->GetLayerNumber() ; ++ln ) {
    if ( LaminatePtr->GetLayerThickness(ln) != 0.0 ){
      mtl::add( SingleLayerEMM( GetLayerMaterialPtr(ln), GetIntPoints(ln) , GetMaterialOrientationOfLayer(ln) ), FinalM );
    }

  }

  return FinalM;
}


////
// Eval mass of an element
////
double LayeredShell::EvalMass() const{

  FELYX_RUNTIME_ASSERT(LaminatePtr != NULL,"LayeredShell::EvalMass()");
  double area_density=0.0;

  // Sum up over all layers
  for (unsigned l=0; l < GetLayerCount(); ++l)
    area_density += GetLayerThickness(l) * GetLayerMaterialPtr(l)->Get("rho");

  return area_density * EvalArea();
}

///
// Get layer number
///
unsigned LayeredShell::GetLayerNr(const LayerEnum layer) const{

  if (layer == MiddleLayer){
    return (unsigned) ( GetLayerCount()-1 / 2) ;
  }
  else if (layer == TopLayer){
    return GetLayerCount()-1;
  }

  // if BottomLayer or anything else
  return 0;
}

Dense_Matrix LayeredShell::evalStrains( const LayerEnum layer ){

  unsigned layerNr = GetLayerNr( layer );

  // Eval strains for these integration points
  return BaseShell::evalStrains( GetIntPoints( layerNr  ) );

}

Dense_Matrix LayeredShell::EvalStrainsAtIntPoints( const unsigned layer ){
  return BaseShell::EvalStrainsAtIntPoints( GetIntPoints( layer ) ); 
}

Dense_Matrix LayeredShell::EvalStressesAtIntPoints( const unsigned layer ){
  MaterialOrientation orient = GetMaterialOrientationOfLayer ( layer );
  Material* mat = GetLayerMaterialPtr( layer );
  return BaseShell::EvalStressesAtIntPoints( GetIntPoints( layer ), mat, orient.GetAngleMaterial1Direction() ); 
}

float_type LayeredShell::evalMaxVonMisesStrain(){

  float_type e=0.0, maxval=0.0;

  // Find max val of bottom layer
  Dense_Matrix Strains = evalStrains( BottomLayer );
  for (unsigned i=0; i < Strains.nrows(); ++i){
    e = vonMises3D( mtl::rows( Strains )[i] );
    if (e > maxval)
      maxval = e;
  }

  // Find max val of top layer
  Strains = evalStrains( TopLayer );
  for (unsigned i=0; i < Strains.nrows(); ++i){
    e = vonMises3D( mtl::rows( Strains )[i] );
    if (e > maxval)
      maxval = e;
  }

  return maxval;

}

Dense_Matrix LayeredShell::evalElementDirectionStrains( const unsigned layer ){
  return BaseShell::evalElementDirectionStrains( GetIntPoints( layer ) );
}

Dense_Matrix LayeredShell::evalMaterialDirectionStrains( const unsigned layer ){
  
  //unsigned layerNr = GetLayerNr( layer );
  
  MaterialOrientation orient = GetMaterialOrientationOfLayer( layer );
  
  return BaseShell::evalMaterialDirectionStrains( GetIntPoints( layer ), orient.GetAngleMaterial1Direction() );
}

Dense_Matrix LayeredShell::evalElementDirectionStresses( const unsigned layer ){
  
  MaterialOrientation orient = GetMaterialOrientationOfLayer ( layer );
  Material* mat = GetLayerMaterialPtr( layer );
  
  return BaseShell::evalElementDirectionStresses( GetIntPoints( layer ), mat, orient.GetAngleMaterial1Direction());
}

Dense_Matrix LayeredShell::evalMaterialDirectionStresses( const unsigned layer ){
  
  MaterialOrientation orient = GetMaterialOrientationOfLayer ( layer );
  Material* mat = GetLayerMaterialPtr( layer );
  
  return BaseShell::evalMaterialDirectionStresses( GetIntPoints( layer ), mat, orient.GetAngleMaterial1Direction());
}

Dense_Matrix LayeredShell::evalMaxStressCriteria( const unsigned layer ){
  
  MaterialOrientation orient = GetMaterialOrientationOfLayer ( layer );
  Material* mat = GetLayerMaterialPtr( layer );
  
  return BaseShell::evalMaxStressCriteria( GetIntPoints( layer ), orient.GetAngleMaterial1Direction(), mat);
}

Dense_Matrix LayeredShell::evalTsaiHillCriteria( const unsigned layer ){
  
  MaterialOrientation orient = GetMaterialOrientationOfLayer ( layer );
  Material* mat = GetLayerMaterialPtr( layer );
  
  return BaseShell::evalTsaiHillCriteria( GetIntPoints( layer ), orient.GetAngleMaterial1Direction(), mat);
}

Dense_Matrix LayeredShell::evalTsaiWuCriteria( const unsigned layer ){
  
  MaterialOrientation orient = GetMaterialOrientationOfLayer ( layer );
  Material* mat = GetLayerMaterialPtr( layer );
  
  return BaseShell::evalTsaiWuCriteria( GetIntPoints( layer ), orient.GetAngleMaterial1Direction(), mat);
}

Dense_Matrix LayeredShell::evalHashinCriteria( const unsigned layer ){
  
  MaterialOrientation orient = GetMaterialOrientationOfLayer ( layer );
  Material* mat = GetLayerMaterialPtr( layer );
  
  return BaseShell::evalHashinCriteria( GetIntPoints( layer ), orient.GetAngleMaterial1Direction(), mat);
}

float_type LayeredShell::evalMaxFailureCriteria(std::string type){
  
  int id;
  if ( type == "MaximumStress")   id = 1;
  else if ( type == "Tsai-Hill")  id = 2;
  else if ( type == "Tsai-Wu")    id = 3;
  else if ( type == "Hashin")     id = 4;
  else std::cerr << "SingleLayerShell: Failure criteria type " << type << " not known!" << endl;
  
  switch(id) {
      case 1: {
        //iterate over all layers
        double max(0.);
        for (unsigned l=0; l < GetLayerCount(); ++l){
          Dense_Matrix maxstress = evalMaxStressCriteria ( l );
          // find the maximum value 
          for (unsigned i = 0; i < maxstress.nrows(); ++i){
            for (unsigned j = 0; j < maxstress.ncols(); ++j){
              if (maxstress(i,j) > max) max = maxstress(i,j);
            }
          }
        }
        return max;
      }
      
      case 2: {
        //iterate over all layers
        double max(0.);
        for (unsigned l=0; l < GetLayerCount(); ++l){
          Dense_Matrix maxstress = evalTsaiHillCriteria ( l );
          // find the maximum value
          for (unsigned i = 0; i < maxstress.nrows(); ++i){
            for (unsigned j = 0; j < maxstress.ncols(); ++j){
              if (maxstress(i,j) > max) max = maxstress(i,j);
            }
          }
        }
        return max;
      }
      
      case 3: {
        //iterate over all layers
        double max(0.);
        for (unsigned l=0; l < GetLayerCount(); ++l){
          Dense_Matrix maxstress = evalTsaiWuCriteria ( l );
          // find the maximum value
          for (unsigned i = 0; i < maxstress.nrows(); ++i){
            for (unsigned j = 0; j < maxstress.ncols(); ++j){
              if (maxstress(i,j) > max) max = maxstress(i,j);
            }
          }
        }
        return max;
      }
      
      case 4: {
        //iterate over all layers
        double max(0.);
        for (unsigned l=0; l < GetLayerCount(); ++l){
          Dense_Matrix maxstress = evalHashinCriteria ( l );
          // find the maximum value
          for (unsigned i = 0; i < maxstress.nrows(); ++i){
            for (unsigned j = 0; j < maxstress.ncols(); ++j){
              if (maxstress(i,j) > max) max = maxstress(i,j);
            }
          }
        }
        return max;
      }   
  }
  return 0;
}

MaterialOrientation LayeredShell::GetMaterialOrientationOfLayer(unsigned n)
{
        if ( LocalMaterialOrientationsVector.size() == 0 )
                return GetLaminatePtr()->GetLayerMaterialOrientation(n);
        else
                return LocalMaterialOrientationsVector[n];
}


