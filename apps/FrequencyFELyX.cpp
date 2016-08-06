//-----------------------------------------------------------------------------
//
// begin     : Dec 6 2001
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

// FELyX Main header
#include "EfObject.h"

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
  string bandwidthAlgorithm     = cl("-bwalgo"  , "sloan");
  int noiselevel                = cl("-noise"   , 1);
  unsigned nEf                  = cl("-nEf"     , 3);
  bool Ev                       = cl("-Ev"      , false);
  
  // list datadir
  if( cl.search(3, "--list", "-ls", "-l") ) {
    cout << endl << "Listing files in datadir: " << datadir << endl;
    string list = "ls " + datadir;
    system( list.c_str() );
    cout << endl;
    exit(0);
  }
  
  cout << endl;
  cout << "START EfObjectTest " << endl;
  cout << "====================================================" << endl;
  
  EfObject FEM( fname, datadir, noiselevel );
  FEM.SaveAnsysModel();

  FEM.NodesReordering( bandwidthAlgorithm );

  std::vector<double> EfRes(nEf, 0.0);
  FEM.CalcEf(nEf,Ev,EfRes);

  cout << endl;
  FEM.PrintGlobalStatus();
  cout << endl;

  cout << "====================================================" << endl;
  cout << "END of EfObjectTest" << endl << endl;

  return 0;
}

void
print_help(const char* Application)
{
  cout << endl;
  cout << "FrequencyFELyX - Frequency Analysis of structures with FELyX" << endl << endl;
  cout << "USAGE:" << endl;
  cout << "--help, -help, -h," << endl;
  cout << "     get some help about this program." << endl << endl;
  cout << "--list, -ls, -l" << endl;
  cout << "     list all files in the data directory " << endl << endl;
  cout << "Type a command line like the following:" << endl;
  cout << Application << " -fname=StructModel-3Dmix-1131e.ansys -datadir=../fedata" << endl << endl;
  cout << "The variables that can be specified are:" << endl;
  cout << "     -fname   filename of FE-model in ANSYS-archive format to be read " << endl;
  cout << "             (default: felyx-3Dmix-2-1131.ansys)" << endl;
  cout << "     -datadir  relative path where this file's located" << endl;
  cout << "             (defualt: ../fedata) " << endl;
  cout << "     -bwalgo  Bandwidth minimization algorithm; choose between: " << endl;
  cout << "              - none " << endl;
  cout << "              - sloan " << endl;
  cout << "              - reversed_sloan " << endl;
  cout << "              - cuthill_mckee " << endl;
  cout << "              - reversed_cuthill_mckee "<< endl;
  cout << "     -noise   Choose amount of information printed to console (from 0=nothing to 2=everything)" << endl;
  cout << "     -nEf     Number of Eigenfrequencies to be computed " << endl;
  cout << "     -Ev      Turn on/off Eigenvector computation. Eigenvectors will be written into eigenvectors.txt in matlab-format. " << endl;
  cout << endl << endl;
}


