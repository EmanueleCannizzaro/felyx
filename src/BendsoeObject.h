/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */ 
//-----------------------------------------------------------------------------
// StructObject.h
//
// begin     : Jan 2 2003
// copyright : (c) 2001 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {okoenig, wintermantel, nzehnder}@imes.mavt.ethz.ch
// www       : www.imes.ethz.ch/st
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

#ifndef BendsoeObject_h
#define BendsoeObject_h BendsoeObject_h

#include "FelyxObject.h"
#include "StructAnalysisType.h"

using namespace std;
using namespace fe_base;
using namespace mtl;

namespace felyx{

 typedef matrix<float_type, symmetric<lower>, envelope<>, row_major >::type EnvelopeMatrix;
 typedef matrix<float_type, rectangle<>, dense<>, column_major>::type Dense_Matrix_1;
 typedef dense1D<float_type> Vector;


  // ! Bendsoe Object class 
  class BendsoeObject:public FelyxObject<StructAnalysisType>{

  public:
    
    // Constructors
    // ------------    
    //! Default constructor with optional noise argument
    BendsoeObject( unsigned _noise = 0 ) :
      FelyxObject<StructAnalysisType>( _noise ){}
    
    //! File constructor, reading in FE model 
    //! Args: filenam, datadir, noiselevel (optional)
    BendsoeObject( string _filenam, string _datadir, unsigned _noise = 0 ):
    FelyxObject<StructAnalysisType>( _filenam, _datadir, _noise ){}

    void PrintGlobalStatus();
    
    int CalcBend();

    int FE(Vector,dense1D<unsigned>,EnvelopeMatrix&,Vector,Vector,Vector);
    
    int calc_xnew(Vector&,Vector,unsigned,vector<unsigned>,double);

    int check(Vector,Vector&,vector<vector<unsigned> >&,vector<unsigned>);

    void PrintFile(Vector,unsigned);

    int topopt(EnvelopeMatrix&, dense1D<unsigned>);

    // Variables
    // ------------   
    int penal;      // penalty-factor
    int loop;       // iterations counter
    unsigned dim;   // dimension of optimization problem
    unsigned Esize; // number of elements
    double rmin;    // radius of the filter

  };

}

#endif
