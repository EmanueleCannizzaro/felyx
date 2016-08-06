//-----------------------------------------------------------------------------
// LcmObject.cc
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

#include "LcmObject.h"

using namespace felyx;

vector<unsigned> LcmObject::Perm(0);

////
//  Loading time depending bc file; optional args: filename, datadir
////
string LcmObject::LoadLcmBcFile( string Fname_, string DataDir_ )
{ 
  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Loading LCM Boundary Conditions from : ";
  string path = LcmBcFileInterface.LoadModel(Fname_, DataDir_, FileBC); 
  if (noise > 0) OutStream << path << endl;
  return path;
}

void LcmObject::PreProcessing() { 

  it = 0;
  solved = 0;
  time =0;
  conjunctionlength=0;
  
  neighbours.clear();
  vents.clear();
  Entraps.clear();
  
  unsigned jj=0;
  for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
    ++jj;
    nodeit->number = jj;
    nodeit->Volume=0;
  }

  // Create list of neighbours
  neighbours.resize(Nodes.size()+1);    
  EvalNeighbours();

  // Detect vents
  DetectVents();

  // Evaluate centroid of elements
  for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ){
    (*eleit)->EvalCentroid();
    (*eleit)->EvalCV();
  }

  // Avoid error
  if (output>Nodes.size()) output=Nodes.size();

  // Remove ancient output files
  if (output>0) system("rm *.dat");
}

void LcmObject::DirectSolver() {
	
  // Set new boundary conditions
  SetNewBC();
  
  if (noise > 0) cout << "----------- Beginning of Timesteps" << endl;

  // Start of time steps
  for (it=0; it<1.5*Nodes.size() && solved==0; ++it){

    if (noise > 1) cout << it << " -------------------------------- time " << time << endl;

    // Set temporary BCs if there are
    if (FileBC.size()>0) SetTdepBC();

    // Remove boundary condition of completely filled nodes
    for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
      if (nodeit->ExistBoundCon() && nodeit->Fill>=1 && !nodeit->ExistOrigBoundCon()){
	nodeit->BoundConPtr = NULL;
      }
    }
    
    // Assemble stiffness matrix and apply boundary conditions
    DofCount = LinkNodes2Gsm(Nodes, Elements);
    dense1D<unsigned> Envelope(DofCount);
    Profile = EvalEnvelope(Envelope, Elements);
    typedef mtl::matrix< float_type , symmetric<lower>, envelope<>, row_major >::type EnvelopeMatrix;
    EnvelopeMatrix GSM(Envelope, Profile );
    AssembleGM(Elements, GSM);
    DofSolution.resize( DofCount );
    ApplyLoads(Nodes, GSM, DofSolution, Envelope);
    
    // Solve set of equations K * p = v in order to obtain pressure distribution
    solved = skyline_solve(GSM, DofSolution, Envelope); 
    if (solved==1) cerr << "Error in Skyline-Solver " << endl;
    
    // Store pressures
    for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
      int GMindex = nodeit->GetGMIndex(0);
      if ( GMindex >=0 )
      nodeit->Pressure = DofSolution[GMindex];
      else
      nodeit->Pressure = 0;
    }

    // Compute element velocities
    for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit)
      if ((*eleit)->GetStatus()) (*eleit)->CompEleVelocity(); 
      else { (*eleit)->vx=0; (*eleit)->vy=0; (*eleit)->vz=0;}

    // Compute volumetric fluxes    
    for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit )
      nodeit->q = 0;

    for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit )
      (*eleit)->CompVolFlux();

    // Check for air entrapments
    if (entrapment>0) UpdateEntrapments(); 

    // Check if there are vents to be closed
    CloseVent();
    if (vents.size()<=0) {solved=1; cout << "vents closed" << endl; }

    // Compute timestep dt
    dt = Comp_dt();
    time += dt;

    // Print current state to file
    if (output > 0) PrintFile();

  } // End of time step

  if (noise > 0) cout << "----------- End of Timesteps" << endl;

  FillingTime = time;
  OutStream << "--< " << my_timer.elapsed() << " >-- " << "LCM simulation done in " << it << " timesteps, time at end of simulation [s]: " << time << endl;
}

void LcmObject::IterativeSolver() {

  typedef mtl::matrix< float_type , symmetric<lower>, envelope<>, row_major >::type EnvelopeMatrix;

  // Assemble stiffness matrix
  DofCount = LinkNodes2Gsm(Nodes, Elements);
  dense1D<unsigned> Envelope(DofCount);
  Profile = EvalEnvelope(Envelope, Elements);
  EnvelopeMatrix GSM_orig(Envelope, Profile );
  EnvelopeMatrix GSM(Envelope, Profile );
  AssembleGM(Elements, GSM);
  DofSolution.resize( DofCount );
  Dense_Vector Forces(DofCount,0);
  symmetric_tag detlef;
  twod_copy(GSM,GSM_orig,detlef);
  
  // Vector that shows which line of the gms corresponds to which node
  vector<unsigned> Nodes_index(DofCount);
  for ( unsigned i=0; i < Nodes.size(); ++i ){
    Nodes_index[Nodes[i].Idx2GM] = i;
  }

  SetNewBC();

  // Start of time steps
  for (it=0; it<1.5*Nodes.size() && solved==0; ++it){

    if (noise > 1) cout << it << " -------------------------------- time " << time << endl;

    if (FileBC.size()>0) SetTdepBC();

    // Remove boundary condition of completely filled nodes
    for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
      if (nodeit->ExistBoundCon() && nodeit->Fill>=1 && !nodeit->ExistOrigBoundCon()){
	nodeit->BoundConPtr = NULL;
      }
    }

    twod_copy(GSM_orig,GSM,detlef);
    ApplyLoads(Nodes, GSM, Forces, Envelope);

    // Update stiffness matrix
    EnvelopeMatrix::iterator matrixit;
    EnvelopeMatrix::OneD::iterator oneDit;
    unsigned i=0,j=0;
    for( matrixit = GSM.begin(); matrixit != GSM.end(); ++matrixit, ++j){
      // if (Nodes[Nodes_index[j]].Fill < 1) {
      if (Nodes[Nodes_index[j]].ExistBoundCon()) {
	for(oneDit = (*matrixit).begin(), i = oneDit.index();  oneDit != (*matrixit).end(); ++oneDit, ++i)
	  *oneDit = 0;
	--oneDit; *oneDit = 1;
      }
      else {
	for(oneDit = (*matrixit).begin(), i = oneDit.index();  oneDit != (*matrixit).end(); ++oneDit, ++i) {
	  // if (Nodes[Nodes_index[i]].Fill < 1){ 
	  if (Nodes[Nodes_index[i]].ExistBoundCon()) {
	    *oneDit = 0;}
	}
      }
    }

    // Solve set of equations K * p = v in order to obtain pressure distribution
    solved = mtl::myCG(GSM, Forces, DofSolution, 1e-6);

    // Store pressures
    for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
      if (nodeit->ExistBoundCon())
	nodeit->Pressure = nodeit->BoundConPtr->getValue(P);
      else 
	nodeit->Pressure = DofSolution[nodeit->Idx2GM];
    }

    // Compute element velocities
    for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit )
      (*eleit)->CompEleVelocity();

    // Compute volumetric flux
    for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit )
      nodeit->q = 0;
    for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit )
      (*eleit)->CompVolFlux();

    // Check for air entrapments
    if (entrapment>0) UpdateEntrapments(); 

    // Check if there are vents to be closed, if yes, build new stiffness matrix
    if (CloseVent()) {
      DofCount = LinkNodes2Gsm(Nodes, Elements);
      Envelope.resize(DofCount);
      Profile = EvalEnvelope(Envelope, Elements);
      EnvelopeMatrix tmp_orig(Envelope, Profile);
      EnvelopeMatrix tmp(Envelope, Profile );
      AssembleGM(Elements, tmp);
      DofSolution.resize( DofCount );
      Forces.resize(DofCount);
      mtl::set(Forces,0);
      symmetric_tag detlef;
      twod_copy(tmp,tmp_orig,detlef);
      GSM = tmp;
      GSM_orig = tmp_orig;

      Nodes_index.resize(DofCount);
      for ( unsigned i=0; i < Nodes.size(); ++i ){
	Nodes_index[Nodes[i].Idx2GM] = i;
      }
    }
    
    // Compute timestep dt
    dt = Comp_dt();
    time += dt;

    // Print current state to file
    if (output > 0) PrintFile();

  } // End of time steps

  FillingTime = time;
  OutStream << "--< " << my_timer.elapsed() << " >-- " << "LCM simulation done in " << it << " timesteps, time at end of simulation [s]: " << time << endl;
}

// Check for vents going through all nodes and write them in vector "vents"
void LcmObject::DetectVents(){

  vector<LcmNode*> temp, temp2;
  vector<LcmNode*>::iterator tempit, tempit2;
  vector<unsigned>::iterator tempit3;
  for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
    if (nodeit->ExistBoundCon()){
      if (nodeit->BoundConPtr->getValue(P)==0 && nodeit->BoundConPtr->getValue(V)==0) temp.push_back(&(*nodeit));
    }
  }
  unsigned i=0;
  while (temp.size() != 0) {
    temp2.clear();
    temp2.push_back(temp[0]);
    vents.push_back(temp2);
    tempit = temp.begin();
    temp.erase(tempit);
    for (unsigned j=0; j < vents[i].size(); ++j) {
      for (tempit2 = temp.begin(); tempit2 != temp.end();) {
	tempit3 = find(neighbours[vents[i][j]->number].begin(), neighbours[vents[i][j]->number].end(), (*tempit2)->number);
	if (tempit3 != neighbours[vents[i][j]->number].end()) {
	  vents[i].push_back(*tempit2);
	  temp.erase(tempit2);
	}
	else ++tempit2;
      }
    }
    ++i;
  }

  if (noise > 0) cout << vents.size() << " Vent(s) detected" << endl;
}

void LcmObject::SetNewBC() {

  for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
    if (nodeit->ExistBoundCon()){
      if (nodeit->BoundConPtr->getValue(P)==0 && nodeit->BoundConPtr->getValue(V)==0) 
	nodeit->setFill(0);
      else {nodeit->setFill(1); }
    }
    else {
      nodeit->setTempBCptr(voidBoundConPtr);
      nodeit->setFill(0);
    }
  }
}

void LcmObject::SetTdepBC() {
  // go through all time dependent boundary conditions
  for (unsigned i=0; i<FileBC.size(); ++i) {

    // if bc is active at the moment
    if (FileBC[i].t_start<=time && FileBC[i].t_end>=time && FileBC[i].t_start>=0) {   
      double pnow = FileBC[i].p_start + (FileBC[i].p_end - FileBC[i].p_start)/
	(FileBC[i].t_end - FileBC[i].t_start)*(time - FileBC[i].t_start);
      LcmBoundCon* lbc;
      lbc = new LcmBoundCon;
      lbc->set(P,pnow); 
      for (unsigned j=0; j<FileBC[i].NodeNr.size(); ++j) {
	// if there is a bc
	LcmNode* nodePtr = &Nodes.at(Perm[FileBC[i].NodeNr[j]-1]); 
	if (nodePtr->ExistOrigBoundCon()) {

	  nodePtr->OrigBoundConPtr->set(P,pnow);
	  if (nodePtr->ExistBoundCon()) {

	    nodePtr->BoundConPtr->set(P,pnow); }
	    }
	else {

	  nodePtr->set(lbc);
	  nodePtr->setFill(1);
	}
      }
    }
    else if (FileBC[i].t_start<0 && !Nodes[Perm[FileBC[i].NodeNr[0]-1]].ExistOrigBoundCon()) {  
      bool go=true;
      for (unsigned j=0; j<FileBC[i].NodeNr.size(); ++j) {
	go = (Nodes[Perm[FileBC[i].NodeNr[j]-1]].Fill>0.9 && go); 
      }
      if (go) {
	LcmBoundCon* lbc;
	lbc = new LcmBoundCon;
	lbc->set(P,FileBC[i].p_start);
	FileBC[i].t_start=time;
	for (unsigned j=0; j<FileBC[i].NodeNr.size(); ++j) {
	  LcmNode* nodePtr = &Nodes[Perm[FileBC[i].NodeNr[j]-1]];
	  nodePtr->set(lbc);
	  nodePtr->setFill(1);
	}
      }
    }
    // if bc is not active at the moment
    else { 
      for (unsigned j=0; j<FileBC[i].NodeNr.size(); ++j) {
	LcmNode* nodePtr = &Nodes[Perm[FileBC[i].NodeNr[j]-1]];
	nodePtr->OrigBoundConPtr=NULL;
	if (nodePtr->Fill==1) nodePtr->BoundConPtr=NULL;
      }
    }
  }
}

void LcmObject::EvalNeighbours() {

  vector<unsigned>::iterator iter;
  for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ){
    unsigned NodeCount = (*eleit)->GetNodeCount();
    for (unsigned i = 0; i<NodeCount; ++i) {
      unsigned temp = (*eleit)->GetNodeIter(i)->number;
      
      for (unsigned j = 1; j<NodeCount; ++j) {
	
	iter = find(neighbours[temp].begin(), neighbours[temp].end(),(*eleit)->GetNodeIter((i+j)%NodeCount)->number);
	if (iter == neighbours[temp].end())
	  { neighbours[temp].push_back((*eleit)->GetNodeIter((i+j)%NodeCount)->number); }
      }    
    }
  }
}

double LcmObject::Comp_dt() {

  double ts=-1.0;
  double temp=1e9;
  vector< std::vector<LcmNode>::iterator > nodeitervec;
  //vector<LcmNode>::iterator tempptr;

  //  if (entrapment>0) solved=1;
  for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
    // if node is about to get filled
    if (nodeit->q>0 && nodeit->Fill<1.0) {
      //      if (entrapment>0) {if (nodeit->entrapment==0) solved=0;}
      ts = (1.0-nodeit->Fill)*nodeit->Volume/nodeit->q;
      if (ts<temp){temp = ts;  }
    }
    // if node is about to get void
    if ( nodeit->q<0 && nodeit->Fill>0 && !nodeit->ExistOrigBoundCon() ) {
      ts = - nodeit->Fill*nodeit->Volume/nodeit->q;
      if (ts<temp){temp = ts;  }
    }
  }
  if (temp==1e9) { temp = 0; solved = 1;}

  if (overfill > 1.0){
    temp *= overfill; //--it;
  }
  // Update fill status
  unsigned a=0;
  for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
    if (nodeit->q>0 && nodeit->Fill<1.0) {
      nodeit->Fill += nodeit->q*temp/nodeit->Volume;
      if (nodeit->Fill >= 0.99999999) { if (nodeit->Fill > 1.0 && a) {++it; ++a; cout << "jetzt hauts " << endl;} nodeit->time=time; nodeit->Fill = 1.0; nodeitervec.push_back(nodeit);}
    }
    bool isentrapneighbour = false;
    for (unsigned u=0; u<neighbours[nodeit->number].size(); ++u)
      {
	if (Nodes[neighbours[nodeit->number][u]-1].entrapment>0 && Nodes[neighbours[nodeit->number][u]-1].Fill<0.1) isentrapneighbour = true; 
      }
    if (nodeit->entrapment>0) {isentrapneighbour = true;}
    if (nodeit->q<0 && ( !nodeit->ExistOrigBoundCon() || nodeit->OrigBoundConPtr->getValue(P)==0  ) && isentrapneighbour) {
      nodeit->Fill += nodeit->q*temp/nodeit->Volume;
      if (nodeit->Fill < 0.0) {nodeit->Fill = 0.0;}
    }
  }

  if (temp > 50*time/it && it>Nodes.size()/5.0) {solved=1; temp=0;}

  if (solved == 0 && !nodeitervec.empty()) {
    if (nodeitervec[0]->q > 0)  { nodeitervec[0]->Fill = 1.0;}
    else { nodeitervec[0]->Fill = 0.0; } }

  // start conjunction
  bool isconjunction=false;
  double u;
  for ( unsigned k=0; k<nodeitervec.size(); ++k ) {
    vector<vector<double> > V;
    for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ) {
      for (unsigned i=0; i<(*eleit)->GetNodeCount(); ++i) {
	if ((*eleit)->GetNodeIter(i) == nodeitervec[k]) {
	  if ((*eleit)->vx!=0 || (*eleit)->vy!=0 || (*eleit)->vz!=0) {
	    vector<double> a(3);
	    a[0]=(*eleit)->vx;
	    a[1]=(*eleit)->vy;
	    a[2]=(*eleit)->vz;
	    V.push_back(a);
	  }
	}
      }
    }
    isconjunction=false;
    if (V.size()>2) {
      for (unsigned i=0; i<V.size()-1; ++i) {
	for (unsigned j=i+1; j<V.size(); ++j) {
	  u = (V[i][0]*V[j][0] + V[i][1]*V[j][1] + V[i][2]*V[j][2]);
	  u /= sqrt(V[i][0]*V[i][0]+V[i][1]*V[i][1] + V[i][2]*V[i][2]) * sqrt(V[j][0]*V[j][0] + V[j][1]*V[j][1] + V[j][2]*V[j][2]);
	  if (u<-0.75) isconjunction=true;
	}
      }
    }
    if (V.size()==2) {
      u = (V[0][0]*V[1][0] + V[0][1]*V[1][1] + V[0][2]*V[1][2]);
      u /= sqrt(V[0][0]*V[0][0]+V[0][1]*V[0][1] + V[0][2]*V[0][2]) * sqrt(V[1][0]*V[1][0] + V[1][1]*V[1][1] + V[1][2]*V[1][2]);
      if (u<-0.75) isconjunction=true;
    }
    if (isconjunction) {conjunctionlength += sqrt(nodeitervec[k]->Volume); 
    // cout << "vol " << nodeptrvec[k]->Volume << ", " << conjunctionlength << endl; 
    //   nodeptrvec[k]->time = 100*u;
    }
  }
  // end counjunction
  return temp;
}

void LcmObject::PrintFile() {

  bool doPrint=false;
  if(output==double(int(output))) { 
    if(it%(Nodes.size()/int(output))==0 )
      doPrint=true; }
  else { 
    if(int(time/output)!=int((time-dt)/output))
      doPrint=true;
  }

  if (solved==1 || it==0) doPrint=true;

  if(doPrint)
    {     
      char buffer[80];
  
      cout << "Write output file ";
      sprintf(buffer, "%04d.dat", it);
      puts(buffer);
      
      string fname = buffer;
      fstream FS(fname.c_str(), ios::out | ios::trunc);
    
      FS << "title = \"LCM-Simulation, time = " << time << " s\" " << endl;
      FS << "variables = \"x\", \"y\", \"z\", \"Fill\", \"Pressure\", \"Time\", \"Flux\", \"Entrapment\", \"bc\", \"temp_bc\", \"vx\", \"vy\", \"vz\" " << endl;
      FS << "zone n=" << Nodes.size() << ", e=" << Elements.size() << ", f=fepoint, et=";
      if ((*Elements.begin())->GetNodeCount()==3) FS << "triangle " << endl;
      else if ((*Elements.begin())->GetNodeCount()==4) FS << "tetrahedron " << endl;
      else if ((*Elements.begin())->GetNodeCount()==8) FS << "brick " << endl;
      FS << endl;
      
      for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
	FS << nodeit->Cx << "\t" << nodeit->Cy << "\t" << nodeit->Cz << "\t" << nodeit->Fill << "\t" << nodeit->Pressure << "\t" << nodeit->time << "\t" << nodeit->q << "\t" << nodeit->entrapment;
	if (nodeit->OrigBoundConPtr == NULL) FS <<  "\t" << "-100000";
	else FS << "\t" << nodeit->OrigBoundConPtr->getValue(P);
	if (nodeit->BoundConPtr == NULL) FS <<  "\t" << "-100000";
	else FS << "\t" << nodeit->BoundConPtr->getValue(P);

	bool go=true;
	for ( eleit = Elements.begin(); eleit != Elements.end() && go; ++eleit ) {
	  if ((*eleit)->GetNodeIter(0) == nodeit ) {
	    FS << "\t" << (*eleit)->vx << "\t" << (*eleit)->vy << "\t" << (*eleit)->vz << endl;
	    go = false;
	  }
	}
	if (go)
	  FS << "\t" << "0" << "\t" << "0" << "\t" << "0" << endl;
      }

      FS << endl;
      
      for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ){
	for (unsigned i=0; i<(*eleit)->GetNodeCount(); ++i){
	  FS << (*eleit)->GetNodeIter(i)->number << "\t";
	}
	FS << endl;
      }
      FS << endl;
      FS.close();     
    }
}

void LcmObject::UpdateEntrapments(){
  vector<LcmNode*>::iterator nodeptrit, nodeptrit2, nodeptrit3;
  int n=EvalEntrapment();
  vector<Entrapment> EntrapsTmp(n);
  // call constructor for each vector entry
  for (int i=0; i<n; ++i) {
    Entrapment a;
    EntrapsTmp[i]=a;
  }
  // if there are or were any entrapments
  if (n>0 || Entraps.size()>0) {
    vector<Entrapment*> Vec1, Vec2;
    bool go=true;
    // evaluate vector EntrapTmp
    for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
      if (nodeit->entrapment>0) { EntrapsTmp[nodeit->entrapment-1].AddMemberNode(&(*nodeit)); }
    }


    // tracking of entrapments
    if (noise > 0) cout << "Total Entrapments " << EntrapsTmp.size() << endl;
    for (unsigned i=0; i<EntrapsTmp.size(); ++i) {
      if (!EntrapsTmp[i].used) {
	Vec1.push_back(&EntrapsTmp[i]); EntrapsTmp[i].used=true;
	for ( nodeptrit = EntrapsTmp[i].MemberNodes.begin(); nodeptrit != EntrapsTmp[i].MemberNodes.end(); ++nodeptrit) {
	  for (unsigned j=0; j<Entraps.size() && go; ++j) {
	    if (!Entraps[j].used) {
	      nodeptrit2 = find(Entraps[j].MemberNodes.begin(), Entraps[j].MemberNodes.end(), *nodeptrit);
	      // if node belongs to both (actual and ancient) entrapments
	      if (nodeptrit2 != Entraps[j].MemberNodes.end()) {
		Vec2.push_back(&Entraps[j]); Entraps[j].used=true; go=false;
		// do check whether entrapment has split up
		for ( nodeptrit2 = Entraps[j].MemberNodes.begin(); nodeptrit2 != Entraps[j].MemberNodes.end(); ++nodeptrit2) {
		  for (unsigned k=i+1; k<EntrapsTmp.size(); ++k) {
		    if (!EntrapsTmp[k].used) {
		      nodeptrit3 = find(EntrapsTmp[k].MemberNodes.begin(), EntrapsTmp[k].MemberNodes.end(), *nodeptrit2);
		      // if entrapment has split up
		      if (nodeptrit3 != EntrapsTmp[k].MemberNodes.end()) {
			Vec1.push_back(&EntrapsTmp[k]); EntrapsTmp[k].used=true;
		      }
		    }
		  }
		}
		
	      }
	    }
	  }
	  go=true;
	}
	CompareEntrapments(Vec1, Vec2);
	Vec1.clear(); Vec2.clear();
      }
    }
    for (unsigned i=0; i<Entraps.size(); ++i)
      if (!Entraps[i].used) Entraps[i].dissolve(); 
    // set used flags to false
    for (unsigned i=0; i<EntrapsTmp.size(); ++i) EntrapsTmp[i].used=false;
    for (unsigned i=0; i<Entraps.size(); ++i) Entraps[i].used=false;
  }   
  
  Entraps.clear();
  Entraps = EntrapsTmp;
}

void LcmObject::CompareEntrapments(vector<Entrapment*>& Vec1_, vector<Entrapment*>& Vec2_){
  vector<LcmNode*>::iterator nodeptrit;

  if (Vec1_.size() == Vec2_.size()) {
    if (noise > 1) cout << "Entrapment update" << endl;
    Vec1_[0]->CompVolumes();
    Vec1_[0]->BoundConPtr->set(P, (Vec2_[0]->BoundConPtr->getValue(P) + entrapment) * Vec2_[0]->Volume / Vec1_[0]->Volume - entrapment);
  if (noise > 1) cout << "Volume = " << Vec1_[0]->Volume << ", Pressure = " << Vec1_[0]->BoundConPtr->getValue(P) << endl;
  Vec1_[0]->SetNodePressure();
  }

  if (Vec1_.size() > Vec2_.size()) {
    if (Vec2_.size()==0) {
      if (noise > 0) cout << "New entrapment has been detected" << endl;
      // Vec1_[0]->Pressure=0;
      Vec1_[0]->CompVolumes();
      Vec1_[0]->SetNodePressure();
    }
    if (Vec2_.size()>0) {
      if (noise > 0) cout << "Entrapment has splitted up" << endl;
      double TotalVolume = 0;
      for (unsigned i=0; i<Vec1_.size(); ++i) {Vec1_[i]->CompVolumes(); TotalVolume += Vec1_[i]->Volume; }
      double NewPressure = (Vec2_[0]->BoundConPtr->getValue(P) + entrapment) * Vec2_[0]->Volume / TotalVolume - entrapment;
      for (unsigned i=0; i<Vec1_.size(); ++i) {Vec1_[i]->BoundConPtr->set(P, NewPressure);
      Vec1_[i]->SetNodePressure();}
    }
  }
  
  if (Vec1_.size() < Vec2_.size()) {
    if (noise > 0) cout << "Join of entrapments" << endl;
    double TotalVolume = 0, TotalVP = 0;
    for (unsigned i=0; i<Vec2_.size(); ++i) {TotalVolume += Vec2_[i]->Volume; TotalVP += Vec2_[i]->Volume*Vec2_[i]->BoundConPtr->getValue(P);}
    Vec1_[0]->CompVolumes();
    double OldPressure = TotalVP / TotalVolume;
    Vec1_[0]->BoundConPtr->set(P, (OldPressure + entrapment) * TotalVolume / Vec1_[0]->Volume - entrapment);
    Vec1_[0]->SetNodePressure();
  }
}


// sets Nodes[i].entrapment to
// 0 if there is a free way from the node to a vent
// 1, 2, ... if node belongs to an air entrapment
int LcmObject::EvalEntrapment() {
  boost::timer timer;
  timer.restart();

  for (unsigned ii=0; ii<Nodes.size(); ++ii){ 
    Nodes[ii].entrapment=-1;}
  int entrap_n = 0;

  // set all nodes with free way to 0
  vector<vector<LcmNode*> >::iterator nodevecptr;
  vector<LcmNode*>::iterator nodeptrit;  
  for (nodevecptr=vents.begin(); nodevecptr!=vents.end(); ++nodevecptr){
    // cout << " -------------_++++++++++++++++============== " << endl;
    for (nodeptrit=(*nodevecptr).begin(); nodeptrit!=(*nodevecptr).end(); ++nodeptrit){ 
      // cout << "einmal " << nodevecptr->size() << endl;
      if ((*nodeptrit)->Fill < 1.0) SetFreeWay(*nodeptrit);
      // cout << "einmalfertig " << endl;
    }
  }

  for (unsigned ii=0; ii<Nodes.size(); ++ii) { 
    if (Nodes[ii].Fill < 0.85 && Nodes[ii].entrapment!=0) { 
      //  cout << Nodes[ii].Fill << " ---- " << endl;
      for (unsigned jj=0; jj<Nodes.size(); ++jj){ 
	if (Nodes[jj].entrapment==-2) Nodes[jj].entrapment=-1;}
      int entrap = EntrapCheck(&Nodes[ii]);
      if (entrap==-3)
	{ ++entrap_n;
	Nodes[ii].entrapment=entrap_n;}
      else
	Nodes[ii].entrapment=entrap;
      // cout << "-----soulution for node " << ii+1 << ", " << Nodes[ii].entrapment << endl;
    }
    // if there is a free way from the node to a vent, set pressure and bc = 0
    if (Nodes[ii].entrapment==0) {
      Nodes[ii].Pressure=0; 
      Nodes[ii].setTempBCptr(voidBoundConPtr);
    }
  }


  return entrap_n;
}

int LcmObject::SetFreeWay(LcmNode* node){
  if (node->entrapment != 0 && node->Fill < 0.85) {
    node->entrapment=0;
    for (unsigned i=0; i<neighbours[node->number].size(); ++i) {
      if (Nodes[neighbours[node->number][i]-1].Fill < 0.85) SetFreeWay(&Nodes[neighbours[(*node).number][i]-1]);
    }
  }
  return 0;
}

int LcmObject::EntrapCheck(LcmNode* node){

  (*node).entrapment=-2;
  int temp=0;
  if ((*node).ExistOrigBoundCon()){
    if ((*node).OrigBoundConPtr->getValue(P) == 0)
      { return 0;}}
  for (unsigned i=0; i<neighbours[(*node).number].size(); ++i) {
    // Ventil ???
    if (Nodes[neighbours[(*node).number][i]-1].ExistOrigBoundCon()){
      if (Nodes[neighbours[(*node).number][i]-1].OrigBoundConPtr->getValue(P) == 0 && Nodes[neighbours[(*node).number][i]-1].Fill < 1)
	{ return 0;}}
    // Einschluss ???
    if (Nodes[neighbours[(*node).number][i]-1].entrapment >= 0)
      {(*node).entrapment = Nodes[neighbours[(*node).number][i]-1].entrapment;
      return Nodes[neighbours[(*node).number][i]-1].entrapment;}
    // einen Knoten weiter ???
    if ( Nodes[neighbours[(*node).number][i]-1].entrapment != -2  && Nodes[neighbours[(*node).number][i]-1].Fill < 1)
      { temp = EntrapCheck(&Nodes[neighbours[(*node).number][i]-1]);
      if (temp >= 0)
	return temp;  }
  }
  return -3;
}

bool LcmObject::CloseVent() {
  vector<vector<LcmNode*> >::iterator nodevecptr;
  vector<LcmNode*>::iterator nodeptrit;  
  bool anyventclosed=false;

  for (nodevecptr=vents.begin(); nodevecptr!=vents.end();){
    bool close=true;
    for (nodeptrit=(*nodevecptr).begin(); nodeptrit!=(*nodevecptr).end(); ++nodeptrit){ 
      if ((*nodeptrit)->Fill < 1.0) close=false;
    }
    if (close) {
      for (nodeptrit=(*nodevecptr).begin(); nodeptrit!=(*nodevecptr).end(); ++nodeptrit){ 
	(*nodeptrit)->BoundConPtr=NULL;
	(*nodeptrit)->OrigBoundConPtr=NULL;
      }
      if (noise > 0) cout << "Vent closed at time " << time << endl;
      vents.erase(nodevecptr);
      anyventclosed=true;
    }
    else ++nodevecptr;
  }
  return anyventclosed;
}

void LcmObject::PostProcessing() {
  TotalVolume=0; FilledVolume=0;
  for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit ){
    TotalVolume += nodeit->Volume;
    FilledVolume += nodeit->Volume * nodeit->Fill;
  }
}

////
// Print global status of actual FEM evaluation
////
void LcmObject::PrintGlobalStatus(){
  PostProcessing();
  OutStream << "# Print Global Status: " << endl;
  OutStream << "\t - Model loaded from        : " << InterfacePtr->GetLoadPath() << endl; 
  OutStream << "\t - # Elements               : " << Elements.size() << endl;
  OutStream << "\t - # Nodes                  : " << Nodes.size() << endl;
  OutStream << "\t - # Materials              : " << Materials.size() << endl;
  OutStream << "\t - # Total Volume           : " << TotalVolume << endl;
  OutStream << "\t - # Filled Volume fraction : " << FilledVolume/TotalVolume << endl;
#ifdef USE_FLOATS
  OutStream << "\t - Precision of val type    : float" << endl;
#else
  OutStream << "\t - Precision of val type    : double" << endl;
#endif
}




