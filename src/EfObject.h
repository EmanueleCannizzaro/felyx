/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */ 
//-----------------------------------------------------------------------------
// StructObject.h
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

#ifndef EfObject_h
#define EfObject_h EfObject_h


#include "FelyxObject.h"
#include "StructAnalysisType.h"


//#include <boost/timer.hpp>

using namespace std;
using namespace fe_base;
using namespace mtl;

namespace felyx{

  // ! Eigenfrequency Object class 
  class EfObject:public FelyxObject<StructAnalysisType>{

  public:
    
    // Constructors
    // ------------    
    //! Default constructor with optional noise argument
    EfObject( unsigned noise_ = 0 ) :
      FelyxObject<StructAnalysisType>( noise_ ){}
    
    //! File constructor, reading in FE model 
    //! Args: filenam, datadir, noiselevel (optional)
    EfObject( string filenam_, string datadir_, unsigned noise_ = 0 ):
      FelyxObject<StructAnalysisType>( filenam_, datadir_, noise_ ){}

    void PrintGlobalStatus();
    
    int CalcEf(unsigned nEf, bool Ev, std::vector<double> &EfRes);

    //! Eval total mass of active elements
    float_type EvalMass() const;  
    
    ////
    // Assembly: Evaluate element mass and fill GMM.
    ////
    template<class Matrix>
    void AssembleGmm(PtrVector<StructElement*>& Elements, Matrix& globMat) {
      
      PtrVector<StructElement*>::iterator eleit;
      for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ){
	if ( (*eleit)->GetStatus() )   	// check if element is active
	  (*eleit) -> AssembleElement2Gmm(globMat);	//Add the EMM to the GMM
      }
    }
    
  };

}

#endif
