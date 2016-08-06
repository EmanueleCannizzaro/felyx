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
#include "StructObject.h"

// Header for parser library
#include "utils/GetPot"

using namespace std;
using namespace felyx;

void print_help(const char*);

int main(int argc, char** argv){
  
  // Command line parsing using GetPot
  // --> Naming all I/O files
  // ---------------------------------
  GetPot cl(argc, argv);
  
  // help options
  if( cl.search(3, "--help", "-help", "-h") || cl.size() == 1 ) {
    print_help(argv[0]);
    exit(0);
  }
  
  // parsing all kinda variables
  string datadir                = cl("-datadir" , "../../fedata");
  string fname                  = cl("-fname"   , "StructModel-3Dmix-1131e.ansys");
  
  std::string bwalgoDefault, solverDefault;
#ifdef HAVE_PARDISO  
  bwalgoDefault = "mmd";
  solverDefault = "pardiso";
#else
  bwalgoDefault = "sloan";
  solverDefault = "skyline";
#endif

  string bandwidthAlgorithm     = cl("-bwalgo"  , bwalgoDefault.c_str() );
  string SolverType             = cl("-solver"  , solverDefault.c_str() );
  int noiselevel                = cl("-noise"   , 1);
  string precond                = cl("-precond" , "rict");
  float_type droptol                = cl("-tol", 1e-6);
  int lfil                      = cl("-lfil", 10);    
  float_type shift                  = cl("-shift", 0);

  // list datadir
  if( cl.search(3, "--list", "-ls", "-l") ) {
    cout << endl << "Listing files in datadir: " << datadir << endl;
    string list = "ls " + datadir;
    system( list.c_str() );
    cout << endl;
    exit(0);
  }
  
  cout << endl;
  cout << "START StructObjectTest " << endl;
  cout << "====================================================" << endl;

  try{
    StructObject FEM( fname, datadir, noiselevel );
    
    FEM.SaveAnsysModel();
    
    // Commented, since the solid elements do not provide this functionality yet!
    //FEM.EvalMass();

    FEM.NodesReordering( bandwidthAlgorithm );
    
    if ( SolverType == "skyline" )
      FEM.DirectSolver();
    else if ( SolverType == "pcg" )
      FEM.IterativeSolver(1e-8, precond, lfil, droptol, shift);
    else if (SolverType == "pardiso" )
      FEM.SparseSolver();
    else
      cerr << "WARNING: No valid solver type specified, evaluating nothing! " << endl;
    

//    if ( StressType == "nodal" || StressType == "element" )
//      FEM.EvalStresses(StressType);
    
    
    FEM.EvalCompliance();
    
    FEM.SaveResults();
    
    cout << endl;
    FEM.PrintGlobalStatus();
    cout << endl;
   
 
  }
  
  // Catch exceptions
  catch(std::exception& e) {
    cerr << e.what() << endl;
  }
  
  cout << "====================================================" << endl;
  cout << "END of StructObjectTest" << endl << endl;

  return 0;
}

void
print_help(const char* Application)
{
  cout << endl;
  cout << "StructFELyX - Structural static analysis with FELyX" << endl << endl;
  cout << "USAGE:" << endl;
  cout << "--help, -help, -h," << endl;
  cout << "  get some help about this program." << endl << endl;
  cout << "--list, -ls, -l" << endl;
  cout << "  list all files in the data directory " << endl << endl;
  cout << "Type a command line like the following:" << endl;
  cout << Application << " -fname=StructModel-3Dmix-1131e.ansys -datadir=../../fedata" << endl << endl;
  cout << "The variables that can be specified are:" << endl;
  cout << " -fname   filename of FE-model in ANSYS-archive format to be read " << endl;
  cout << "         (default: felyx-3Dmix-2-1131.ansys)" << endl;
  cout << " -datadir  relative path where this file's located" << endl;
  cout << "         (default: ../fedata) " << endl;
  cout << " -bwalgo  Bandwidth minimization algorithm; choose between: " << endl;
  cout << "      - none " << endl;
  cout << "      - sloan " << endl;
  cout << "      - reversed_sloan " << endl;
  cout << "      - cuthill_mckee " << endl;
  cout << "      - reversed_cuthill_mckee "<< endl;
  cout << " -solver  Choose between: " << endl;
  cout << "      - skyline " << endl;
  cout << "        Direct skyline solver, based on MTL envelope storage format, using ATLAS dot-product if available" << endl;
  cout << "      - pardiso " << endl;
  cout << "        Direct solver PARDISO from the university of Basel, Switzerland." << endl;
  cout << "        Taken from Intel's MKL v7.0, which therefore must be available at compile time!" << endl;
  cout << "      - pcg " << endl;
  cout << "        Preconditioned Gradient Solver, calling ITL solver " << endl;
  cout << "      -noise   Choose amount of information printed to console (from 0=nothing to 2=everything)" << endl;
  cout << "      -stress  Type of Stress computation (none, nodal, element)" << endl;
  cout << endl << endl;
}


