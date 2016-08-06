//-----------------------------------------------------------------------------
// StructObjectTest.cpp
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// FELyX Main header
#include "DucFelyx_solution.h"

int main( int ac, char* av[] ) {

  try {

    std::cout << "FE analysis of Ducati steel trellis frame using FELyX\n";


    DucFelyxObject FEM("DucatiFrameSolution.ansys", "../../tutorial", 2 );

    // Eval mass of frame
    double mass = FEM.EvalMass();

    FEM.ApplyBrakingLoadcase();

    // Solve FE problem
    FEM.SparseSolver();

    double MaxStress= FEM.EvalMaxStress() ;
    
    // Apply torsion loadcase
    FEM.ApplyTorsionLoadcase();

    // Solve FE problem
    FEM.SparseSolver();

    double Stiffness= FEM.EvalTorsionStiffness();


  }

  // Catch exceptions
  catch ( std::exception & e ) {
    cerr << e.what() << endl;
    return 1;
  }

  std::cout << "FELyX analysis done\n";

  return 0;
}
