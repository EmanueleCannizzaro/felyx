//-----------------------------------------------------------------------------
// IOLcmBc.cc
//
// begin     : Nov 24 2003
// copyright : (c) 2003 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
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
   
#include "IOLcmBc.h"

////
//// Constructors
////
IOLcmBc::IOLcmBc( vector<LcmNode>&		Nodes_,
		  vector<CoordSys>& 	        NodeCoordSysList_,
		  PtrVector<LcmElement*>&        Elements_,
		  PtrVector<Material*>&	        Materials_,
		  vector<PropertySet>&	        Properties_,
		  vector<LcmBoundCon>& 		BoundaryConditions_,
		  vector<Layer>&                    Layers_,
          vector<Laminate>&                 Laminates_,
		  string 			        Fname_, 
		  string 			        DataDir_ ) 
  : IOFelyx<LcmAnalysisType>( Nodes_, NodeCoordSysList_, Elements_, Materials_, Properties_, BoundaryConditions_, Layers_, Laminates_ ) {
  SetFileNames( Fname_, DataDir_, ".ansys", ".ansysnew", ".felyxres" );
}

////
//// Member functions
////

string IOLcmBc::LoadModel( string Fname_, string DataDir_, vector<TdepBC>& FileBC) { 

  ifstream 		IN;			// File stream
  string line;
  string word;

  // Modify File name settings, if necessary
  SetFileNames( Fname_, DataDir_ );
  
  // Create path of file to be openend
  string path = GetLoadPath();

  // Open File to read
  IN.open( path.c_str(), ios::in );
  if( !IN.is_open() ){
    cout << "ERROR in IOLcmBct::LoadModel()" << endl;
    cout <<"Couldn't open file: "<< path << endl;
    exit(1);
  }

  // Go to begin of file
  IN.clear();
  IN.seekg( 0, ios::beg );
  
  unsigned i=0;
  vector<unsigned> position, n_nodes;
  // First loop through file
  while(!IN.eof()){
    if(getline(IN,line)) {
      position.push_back(IN.tellg());
    }
    if(getline(IN,line)) {
      position.push_back(IN.tellg());
      bool go=true;
      int pos=0, n=0;
      while(go) {
	pos=line.find(" ",pos+1);
	if (pos==-1) go=false;
	++n;
      }
      n_nodes.push_back(n);
    }
    ++i;
  }

  // Go to begin of file
  IN.clear();
  IN.seekg( 0, ios::beg );
  
  // Second loop through file
  unsigned node_nr;
  for (unsigned k=0; k<i-1; ++k){
    TdepBC bc;
    bc.NodeNr.resize(0);
    if(IN >> bc.p_start >> bc.p_end >> bc.t_start >> bc.t_end) {}
    IN.ignore(10000,'\n');
    for (unsigned j=0; j<n_nodes[k]; ++j) {
      if(IN >> node_nr) {}
      (bc.NodeNr).push_back(node_nr);
    }
    FileBC.push_back(bc);
  }
//   while(!IN.eof()){
//     LcmBoundCon* bcPtr;
//     bcPtr = new LcmBoundCon;
//     if(IN >> bcPtr->p_start >> bcPtr->p_end >> bcPtr->t_start >> bcPtr->t_end) {}
//     IN.ignore(10000,'\n');
//     for (unsigned j=0; j<n_nodes[i]; ++j) {
//       if(IN >> node_nr) {}
//       Nodes[node_nr].set(bcPtr);
//     }
//     ++i;
//   }
  return path;
}
