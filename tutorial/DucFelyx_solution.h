//-----------------------------------------------------------------------------
// DucFelyx.h
//
// begin     : Jan 2005
// copyright : (c) 2005 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {koenig,wintermantel,zehnder}@even-ag.ch
// www       : even-ag.ch
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

#ifndef DucFelyx_h
#define DucFelyx_h DucFelyx_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

// Header for standard mathematical functions
#include <cmath>

// include FELyX object
#include "StructObject.h"

using namespace std;
using namespace felyx;

extern const double PI;

//! Definition of class DucFelyxObject
/*! Class is derived from StructObject, implementing all necessary FEM functionality
    for the ducati frame optimization. 
*/
class DucFelyxObject : public StructObject {

  public:

    // Constructors
    // ------------
    //! Empty constructor, doing nothing
    DucFelyxObject() : StructObject() {}

    //! File constructor, reading in FE model 
    //! Args: filenam, datadir, noiselevel (optional)
    DucFelyxObject( std::string filenam_, std::string datadir_, unsigned noise_ = 0 )
     : StructObject( filenam_, datadir_, noise_){}

    // Functions to apply different loadcases
    // --------------------------------------
    //! Function that deletes all forces and moments and applies the torsion loadcase
    void ApplyTorsionLoadcase();
    //! Function that deletes all forces and moments and applies the braking loadcase
    void ApplyBrakingLoadcase();

    // Eval functions
    // --------------

    //! Eval mass of structure in a range of element reference numbers in kg
    double EvalMass() const;

    //! Run FE-analysis for torsion loadcase
    //! Eval torsion stiffness of frame, including engine in Nm/degree
    double EvalTorsionStiffness() const;

    //! Eval maximum stress of structure, looking at elements in range of reference numbers
    double EvalMaxStress();

    //! Eval and set Beam properties for a certain Ri and t
    /*! @param tube	index of tube
        @param Ri		inner radius
        @param t		thickness
        @param final	bool to tell function if propery set is used for later ANSYS evaluations (=true)
    	or just for further optimizations within EO/FELyX (false)
    */
    void SetTubeProperties( unsigned tube, double Ri, double t, bool final );

    //! Return number of tubes that can be varied
    unsigned GetVarTubesCount() { return nVarTubes; }

  private:
    //! Number of tubes to be optimized - defined through parser in constructor
    unsigned nVarTubes;

};



#endif
