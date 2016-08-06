//-----------------------------------------------------------------------------
// IOPatranFormatModel.inl
//
// begin     : Mar 18 2001
// copyright : (c) 2001 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {okoenig, wintermantel, nzehnder}@imes.mavt.ethz.ch
// www       : www.imes.ethz.ch/st
//
// $Id: IOPatranLoadModel.inl,v 1.1.1.1 2004/12/13 17:55:16 okoenigl Exp $
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
   
#include "IOPatranFormat.h"

////
//// Member functions
////


template<class analysis_type>
void IOPatranFormat<analysis_type>::readOutFile(string path) {
  
  ifstream IN;
  
  // Open File to read
  IN.open( path.c_str(), ios::in );
  if( !IN.is_open() ){
    cout << "ERROR in IOPatranFormat::LoadModel()" << endl;
    cout <<"Couldn't open file: "<< path << endl;
    exit(1);
  }

  map<int, int> node2CID;
  map<int, int> elem2PID;
  map<int, int> elem2coord;
  map<int, int> bc2node;
  map<int, int> CID2coord;
  
  map<int, int> elemTypeMap;
  map<string, typename analysis_type::bla> bcTypeMap;
  
  elemTypeMap[3] = Darcy3D3::Id;
  elemTypeMap[5] = Darcy3D4::Id;
  elemTypeMap[8] = Darcy3D8::Id;
  
  bcTypeMap["HEAT"] = (typename analysis_type::bla)V;
  bcTypeMap["TEMP"] = (typename analysis_type::bla)P;
  
  // index, tmp vars
  int i;
  string line;
  
  CoordStore coord(&NodeCoordSysList);

  // packet header in .out file
  int it, id, iv, kc, n1, n2, n3, n4, n5;

  // summary (packet 26)
  int numNode, numElem, numMat, numProp, numCoord;

  IN.clear();
  getline(IN, line);
  while(!IN.eof()) {

    // read packet header
    LineReader reader(line);
    it = reader.Int(2);
    id = reader.Int(8);
    iv = reader.Int(8);
    kc = reader.Int(8);
    n1 = reader.Int(8);
    n2 = reader.Int(8);
    n3 = reader.Int(8);
    n4 = reader.Int(8);
    n5 = reader.Int(8);
    
    // first id = 1 -> 0
    id--;

    switch(it) {
      case 1: /* node */ {
        double x, y, z;
        int cid;
        
        assert(kc == 2); // # of data cards
  
        getline(IN, line);
        reader = LineReader(line);
        x = reader.Double(16);
        y = reader.Double(16);
        z = reader.Double(16);
        
        getline(IN, line);
        reader = LineReader(line);
        reader.skip(1+1+8+8);
        cid = reader.Int(8);  
        cid--;
        if(cid > 0) node2CID[id] = cid;
        
        Nodes[id].set(x, y, z);
               
      } break;
        
      case 2: /* element */ {
        int nodes, pid;
        int node;
        
        double theta1, theta2, theta3;
  
        getline(IN, line);
        reader = LineReader(line);
        nodes = reader.Int(8);
        reader.skip(8);
        pid = reader.Int(8);
        reader.skip(8);
        theta1 = reader.Double(16);
        theta2 = reader.Double(16);
        theta3 = reader.Double(16);
     
        elem2coord[id] = coord.add(theta1, theta2, theta3);
        
        kc--;
        
        elem2PID[id] = pid; // pid++/-- !!
        
        typename PtrVector<typename analysis_type::element_type*>::iterator eleit = Elements.begin();
        eleit += id;
        create_element(eleit, elemTypeMap[iv]);
        
        while(nodes > 0) {
          getline(IN, line);
          reader = LineReader(line);
          kc--;
          for(i=0;(i < 10) && (nodes > 0);i++,nodes--) {
            node = reader.Int(8);
            node--;
            Elements[id]->SetNodeIter(i, Nodes.begin()+node);
          }
        }
        
        while(kc-- > 0) {
          getline(IN, line);
        }
  
      } break;
  
      case 3: /* material */ {
        Material *mat = NULL;
        double rho = 0, t = 0, e11 = 0, e22 = 0, e33 = 0;
        
        i = 0;
        while(kc-- > 0) {
          getline(IN, line);
          switch(i) {
            case 0:
              reader = LineReader(line);
              reader.skip(16);
              rho = reader.Double(16);
              t = reader.Double(16);
              break;
            case 5:
              reader = LineReader(line);
              reader.skip(16);
              e11 = reader.Double(16);
              e22 = reader.Double(16);
              e33 = reader.Double(16);
              break;
            default:
              break;
          }
          i++;
        }

        switch(iv) {
          case 1: /* isotropic */
            mat = new IsotropicMaterial(0 ,0 , 0, e11, rho);
            break;
          
          case 3: /* 2d/3d orthotropic */
            mat = new OrthotropicMaterial(0, 0, 0, 
	    				  0, 0, 0, 
					  0, 0, 0, 
					  0,
					  e11, e22, e33,
					  rho);
            break;

          case 2: /* 2d anisotropic */
          case 6: /* 3d anisotropic */
          default:
              cout << "  Unknown Material Type: " << iv << endl;
              exit(1);
            break;
        }
        
        if(mat != NULL) {
          Materials[id] = mat;
        } else {
          cout << "  Could not Create Material" << endl;
          exit(1);
        }
        
      } break;
                assert(kc == 1); // # of data cards

      case 4: /* property */ {
        double value;
  
        getline(IN, line);
        reader = LineReader(line);
        reader.skip(16);
        value = reader.Double(16);
        kc--;
  
        iv--;
        Properties[id] = PropertySet(Materials[iv], id);
        Properties[id].Set("Thickness", value);

        while(kc-- > 0) {
          getline(IN, line);
        }
      } break;
      
      case 5: /* coordinate frame */ {
      
        double r11, r21, r31, r12, r22, r32, r13, r23, r33;

        assert(kc == 4); // # of data cards
        
        if(iv != 1) {        
          cout << "Spherical/Cylindrical CoordSys are not supported." << endl;
          exit(1);
        }
                
        getline(IN, line);
        getline(IN, line);
        reader = LineReader(line);
        reader.skip(4*16);
        r11 = reader.Double(16);
        getline(IN, line);
        reader = LineReader(line);
        r21 = reader.Double(16);
        r31 = reader.Double(16);
        r12 = reader.Double(16);
        r22 = reader.Double(16);
        r32 = reader.Double(16);
        getline(IN, line);
        reader = LineReader(line);
        r13 = reader.Double(16);
        r23 = reader.Double(16);
        r33 = reader.Double(16);
        
        Dense_Matrix rotmat(3,3);
        
        rotmat(0, 0) = r11;
        rotmat(1, 0) = r21;
        rotmat(2, 0) = r31;
        rotmat(0, 1) = r12;
        rotmat(1, 1) = r22;
        rotmat(2, 1) = r32;
        rotmat(0, 2) = r13;
        rotmat(1, 2) = r23;
        rotmat(2, 2) = r33;
        
        NodeCoordSysList.push_back(CoordSys(CoordSys::cartesian, CoordSys::rad, CoordSys::Rotvec, CoordSys::Mat2Rotvec(rotmat)));
        CID2coord[id] = NodeCoordSysList.size()-1;
        
      } break;
  
      case /* boundary cond */ 10: {
        double value;
        
        assert(kc == 1); // # of data cards
  
        getline(IN, line);
        reader = LineReader(line);
        value = reader.Double(16);
        
        typename analysis_type::bc_type bc;
        bc.set(bcTypeMap["TEMP"], value);
        BoundaryConditions.push_back(bc);
        bc2node[BoundaryConditions.size()-1] = id;
        
      } break;
        
      case 25: /* title */ {
        assert(kc == 1); // # of data cards
        while(kc-- > 0) {
          getline(IN, line);
        }
      } break;
          
      case 26: /* summary */ {
        assert(kc == 1); // # of data cards
        
        numNode =  n1;
        numElem =  n2;
        numMat =   n3;
        numProp =  n4;
        numCoord = n5;
  
        // init
        Nodes.resize(numNode);
        Elements.resize(numElem);
        Materials.resize(numMat);
        Properties.resize(numProp);
        
        NodeCoordSysList.clear();
        BoundaryConditions.clear();
  
        while(kc-- > 0) {
          getline(IN, line);
        }
      } break;
          
      case 99: /*  of file */ {
        while(kc-- > 0) {
          getline(IN, line);
        }
        
      } break;
  
      // unknown
      default: {
        cout << "  Ignoring unknown Packet of Type " << it << endl;
        while(kc-- > 0) {
          getline(IN, line);
        }
      } break;
    }

    getline(IN, line);
  }

  IN.close();
  
  // insert missing pointers
  
  map<int, int>::iterator iter;
  
  iter = node2CID.begin();
  while(iter != node2CID.end()) {
    Nodes[(*iter).first].set(&(NodeCoordSysList[CID2coord[(*iter).second]]));
    iter++;
  }
    
  iter = elem2coord.begin();
  while(iter != elem2coord.end()) {
    Elements[(*iter).first]->SetEleCoordSysPtr(&(NodeCoordSysList[(*iter).second]));
    iter++;
  }

  iter = elem2PID.begin();
  while(iter != elem2PID.end()) {
    int pid = (*iter).second;
    if(pid < 0) {
      pid = 1 - pid;
      Elements[(*iter).first]->SetMaterialPtr(Materials[pid]);
    } else {
      pid--;
      Elements[(*iter).first]->SetPropertiesPtr(&(Properties[pid]));
      Elements[(*iter).first]->SetMaterialPtr(Properties[pid].MaterialPtr);
    }
    iter++;
  }
  
  iter = bc2node.begin();
  while(iter != bc2node.end()) {
    Nodes[(*iter).second].set(&BoundaryConditions[(*iter).first]);
    iter++;
  }
    
}

template<class analysis_type>
void IOPatranFormat<analysis_type>::readLayupFile(string path) {
  
  ifstream IN;
    
  // Open File to read
  IN.open( path.c_str(), ios::in );
  if( !IN.is_open() ){
    IN.close();
    cout << "  No Layup file found." << endl;
    return;
  } else {
    cout << "  Load Data from " << path << endl;
  }
    
  string line;
  istringstream iss;
  
  int code, elem, nil;
  double warp, weft, double_nil;
  
  IN.clear();
  getline(IN, line);
  while(!IN.eof()) {
    iss.str(line);
    if(line.size() > 3) {
      if(line.substr(0,3) == "430") {
        iss >> code >> nil >> elem >> nil >> double_nil >> warp >> weft;
        elem--;
        if( dynamic_cast<LcmElement*>(Elements[elem]) != NULL ) {
          dynamic_cast<LcmElement*>(Elements[elem])->setWarp(warp);
          dynamic_cast<LcmElement*>(Elements[elem])->setWeft(weft);
        }
      }
    }
    getline(IN, line);
  }  
 
  IN.close();
}

////
//// Implementation of the Member function LoadModel of class IOPatranFormat
////
template<class analysis_type>
string IOPatranFormat<analysis_type>::LoadModel ( string Fname_, string DataDir_ ) {

  string path;
  string FnameLayup;
  
  FnameLayup = Fname_;
  FnameLayup.replace(FnameLayup.rfind(".out"), 6, ".Layup");
  
  // .out file
  SetFileNames( Fname_, DataDir_ );
  path = GetLoadPath();
  
  cout << path << endl;
  readOutFile(path);
  cout << "  Finished loading " << path << endl;
  
  // .Layup file
  SetFileNames( FnameLayup, DataDir_ );
  path = GetLoadPath();
  readLayupFile(path);

  cout << "  Finished loading ";
  return path;
}


