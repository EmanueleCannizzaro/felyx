//-----------------------------------------------------------------------------
// Entrapment.h
//
// begin     : Dec 6 2001
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
   

#ifndef Entrapment_h
#define Entrapment_h Entrapment_h

using namespace std;
using namespace mtl;


namespace fe_base{	//begin of namespace fe_base
  
  class Entrapment{
    
  public:
    // Data members
    // ------------
    vector<LcmNode*> MemberNodes;
    double Volume;		
    bool used;	
    LcmBoundCon *BoundConPtr;	
    
    // Constructors
    // ------------
    Entrapment(){
      used = false;
      BoundConPtr = new LcmBoundCon;
      BoundConPtr->set(P,0);
    };
    ~Entrapment(){
      // delete BoundConPtr;
    };
    
    // Sets and gets
    // -------------
    
    void CompVolumes(){
      vector<LcmNode*>::iterator nodeit;
      Volume=0;
      for ( nodeit = MemberNodes.begin(); nodeit != MemberNodes.end(); ++nodeit ){
	Volume += (*nodeit)->Volume * (1.0 - (*nodeit)->Fill);
      }
    };

    void AddMemberNode(LcmNode* node){
      MemberNodes.push_back(node);
    };
    
    void SetNodePressure(){
      vector<LcmNode*>::iterator nodeit;
      for ( nodeit = MemberNodes.begin(); nodeit != MemberNodes.end(); ++nodeit ){
	(*nodeit)->Pressure = BoundConPtr->getValue(P);
	if ((*nodeit)->Fill < 0.8) {
	  (*nodeit)->setTempBCptr(BoundConPtr);

	}
      } 
    };

    void dissolve() {
      BoundConPtr = NULL;
      vector<LcmNode*>::iterator nodeit;
      for ( nodeit = MemberNodes.begin(); nodeit != MemberNodes.end(); ++nodeit ){
	(*nodeit)->setTempBCptr(NULL);
      }
    };

    void PrintNodes(){
      vector<LcmNode*>::iterator nodeit;
      cout << "Print all Nodes" << endl;
      for ( nodeit = MemberNodes.begin(); nodeit != MemberNodes.end(); ++nodeit ){
	cout << ", " << (*nodeit)->number;
      }
      cout << endl;
    };
    
  };
  
  
}	// end of namespace fe_base

#endif
