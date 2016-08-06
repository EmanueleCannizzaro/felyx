//-----------------------------------------------------------------------------
// IOAnsysSaveModel.cc
//
// begin     : Mar 18 2001
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



////
//// Implementation of the Member function SaveModel of class IOAnsysFormat
////
template<class analysis_type>
string IOAnsysFormat<analysis_type>::SaveModel ( string Fname_, string DataDir_ ) {

  //// -------------
  //// File handling
  //// -------------

  // Modify File name settings, if necessary
  SetFileNames( Fname_, DataDir_ );

  // Create path of file to be openend
  string path = GetSavePath();

  ofstream OUT(path.c_str());

  //// -----------------
  //// Print File Header
  //// -----------------

  OUT<<"/COM, FELyX FE Model, exported to ANSYS archive format"<<endl;
  OUT<<"/COM, ------------------------------------------------"<<endl;
  OUT<<"/COM,"<<endl<<endl;

  OUT<<"FINISH"<<endl;
  OUT<<"/CLEAR,START"<<endl;

  OUT<<"/PREP7"<<endl;
  OUT<<"/NOPR"<<endl;
  OUT<<"/TITLE, Model exported by FELyX Toolbox - File: " << path << endl;


  //// -----------------------------
  //// Print Ansys Scalar Parameters
  //// -----------------------------

  OUT << endl << "! SCALAR PARAMETERS " << endl;
  OUT <<         "! ----------------- " << endl;
  map<string,double>::const_iterator it = AnsysParameters.begin();
  while (it != AnsysParameters.end() ){
    OUT <<"*SET, " << it->first << " , " << it->second << endl;
    it++;
  }


  //// -------------------
  //// Print Element Types
  //// -------------------

  OUT << endl << "! ELEMENT TYPES " << endl;
  OUT <<         "! ------------- " << endl;

  typename PtrVector<typename analysis_type::element_type*>::const_iterator eleit_begin, eleit;
  vector<unsigned> refnumbers;
  unsigned refnumber, Id;

  for( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ){
    if ( (*eleit)->GetStatus() ){       // check if element is active
      refnumber = (*eleit)->GetRefNumber();
      Id = (*eleit)->GetId();

      // Store each element reference number only once
      if ( find(refnumbers.begin(), refnumbers.end(), refnumber) == refnumbers.end() ) {
        refnumbers.push_back(refnumber);

        if ( Id == 393 ) Id = 93;
        if ( Id == 391 ) Id = 91;

        OUT <<"ET,  " << refnumber << " , " << Id << endl;
      }
    }
  }
  OUT << endl;


  //// -------------------
  //// Print Key Options
  //// -------------------

  OUT << endl << "! Key Options " << endl;
  OUT <<         "! ----------- " << endl;

  unsigned  ln = 0;
  map<unsigned,unsigned> refnumber2maxLn;

  for( eleit = Elements.begin(); eleit != Elements.end(); ++eleit ){
    if ( (*eleit)->GetLaminatePtr() != NULL ){
      ln = (*eleit)->GetLaminatePtr()->GetLayerNumber();
      refnumber = (*eleit)->GetRefNumber();
      if ( refnumber2maxLn.find( refnumber ) == refnumber2maxLn.end() )
        refnumber2maxLn[refnumber] = ln;
      else {
        if ( ln > refnumber2maxLn[refnumber] )
          refnumber2maxLn[refnumber] = ln;
      }
    }
  }

  map<unsigned,unsigned>::iterator myIter = refnumber2maxLn.begin();

  while ( myIter != refnumber2maxLn.end() ){
    if (myIter->second < 16) myIter->second = 16;
    OUT << "KEYOP, " << myIter->first << ", 1, " << myIter->second;
    OUT << endl << "KEYOP, " << myIter->first << ", 2, 1";
    ++myIter;
  }

//   if ( maxLn != 0 ){
//     if (maxLn < 16) maxLn = 16;
//     OUT << "KEYOP, " << refnumber << ", 1, " << maxLn;
//     OUT << endl << "KEYOP, " << refnumber << ", 2, 1" ;
//   }


  //// -------------------
  //// Print RealConstants
  //// -------------------
  OUT << endl << endl << "! REAL CONSTANTS" << endl;
  OUT <<         "! --------------" << endl;

  //Find the number of properties used in this model
  const int PropCount = Properties.size();

  map<unsigned, string>::reverse_iterator rIter;
  map<unsigned, string>::iterator Iter;

  unsigned RCEntryCount = 0, MaxRCEntryCount = 0;

  //find MaxRCEntryCount
  for ( int i = 0 ; i < PropCount ; ++i ){
    if ( Properties[i].GetId() == 0 )
      cerr << endl << "ERROR writing PropertySet: The Id of the Property is 0!! " << endl;
    else {
      rIter = ElementNr2RCSetPosition[Properties[i].GetId()].rbegin();
      RCEntryCount = rIter->first;
      if ( RCEntryCount >= MaxRCEntryCount )
        MaxRCEntryCount = RCEntryCount+1;
    }
  }

  unsigned TotalPropCount = 0;
  //find MaxRCEntryCount in the Laminates
  for (unsigned i = 0 ; i < Laminates.size() ; ++i ){
    ++TotalPropCount;
    unsigned ln = Laminates[i].GetLayerNumber();
    RCEntryCount = 12+(6*ln);
    if ( RCEntryCount > MaxRCEntryCount )
      MaxRCEntryCount = RCEntryCount;
  }

  TotalPropCount += PropCount;

  //write the first lines of the RLBLOCK command
  OUT << "RLBLOCK," << setw(7) << TotalPropCount << "," <<
    setw(7) << TotalPropCount << "," << setw(7) <<
    MaxRCEntryCount << "," << setw(7) << "7\n" <<
    "(2i8,6g16.9)\n" << "(7g16.9)";

  //loop over all Properties
  for ( int i = 0 ; i < PropCount ; i++ ){

    rIter = ElementNr2RCSetPosition[Properties[i].GetId()].rbegin();
    RCEntryCount = rIter->first;
    if ( RCEntryCount > MaxRCEntryCount )
      MaxRCEntryCount = RCEntryCount;

    OUT.setf(ios::right, ios::adjustfield);

    Iter = ElementNr2RCSetPosition[Properties[i].GetId()].begin();

    //Every RC Set has at least 6 entrys
    OUT << endl;
    OUT  << setw(8) <<  i+1; //first in the row is the number of the RCSet
    if (RCEntryCount > 6)
      OUT << setw(8) << RCEntryCount+1 << " "; //second is the number of RC entrys
    else
      OUT << setw(8) << "6" << " "; //second is the number of RC entrys

    OUT.precision(9);
    OUT.setf(ios::left, ios::adjustfield);
    OUT.setf(ios::showpoint);

    for ( unsigned k = 0 ; k < 6 ; k++ ){
      if ( Iter->first == k ){
        OUT << setw(16) << Properties[i].GetDoubleMap()[Iter->second];
        Iter++;
      }
      else {
        OUT << setw(16) << 0.0 ;
      }
    }
    //... some have more
    if ( RCEntryCount > 6 ){
      OUT << '\n' << " ";
      for ( unsigned k = 6 ; k <= RCEntryCount ; k++ ){
        if ( Iter->first == k ){
          OUT << setw(16) << Properties[i].GetDoubleMap()[Iter->second];
          Iter++;
        }
        else {
          OUT << setw(16) << 0.0 ;
        }
      }
    }
  }

  //Find the number of layered materials, which have to be written back
  //to propertysets...
  unsigned  LaminateCount = PropCount, LayerCount = 0;
  RCEntryCount = 0;
  for (unsigned i = 0 ; i < Laminates.size(); ++i ){
      Laminate* myLaminatePtr = &Laminates[i];
      LayerCount = myLaminatePtr->GetLayerNumber();
      ++LaminateCount;

      //Calculate the number of values
      //The last 3 are all 0 and must not be written
      RCEntryCount = 12 + (6*LayerCount)-3;
      //Writing the data
      OUT.setf(ios::right, ios::adjustfield);
      OUT << endl;
      OUT << setw(8) << LaminateCount;
      OUT << setw(8) << RCEntryCount << " ";

      OUT.precision(9);
      OUT.setf(ios::left, ios::adjustfield);
      OUT.setf(ios::showpoint);

      //write the first 12 values
      OUT << setw(16) << double(LayerCount);
      //There are some parameters provided by ansys which are not
      //used within FELyX
      for (int k = 0 ; k < 5 ; ++k )
        OUT << setw(16) << 0.0;
      OUT << endl;
      for (int k = 0 ; k < 6 ; ++k )
        OUT << setw(16) << 0.0;

      //write the layerwise data
      //LineCount organizes the endl statements
      unsigned LineCount = 7;
      for (unsigned k = 0 ; k < LayerCount ; ++k ){

        //Find the materialnumber for the layer
        Material* LayerMaterialPtr = myLaminatePtr->GetLayerMaterialPtr(k);
        vector<Material*>::iterator myItr = Materials.begin();
        double MatNo = 1.0 ;
        while((*myItr) != LayerMaterialPtr && myItr != Materials.end() ){
          ++MatNo;
          ++myItr;
        }

        OUT << setw(16) << MatNo;
        //      ++LineCount;
        if ( LineCount++ % 7  == 0 )
          OUT << endl;
        OUT << setw(16) << myLaminatePtr->GetLayerPtr(k)->GetAngle();
        if ( LineCount++ % 7  == 0 )
          OUT << endl;
        OUT << setw(16) << myLaminatePtr->GetLayerPtr(k)->GetThickness();

        //dont write this for the last layer
        if ( k != LayerCount-1 ){
          if ( LineCount++ % 7  == 0 )
            OUT << endl;
          for (unsigned l = 0 ; l < 3 ; ++l){
            OUT << setw(16) << 0.0;
          if ( LineCount++ % 7  == 0 )
            OUT << endl;
          }
        }
      }
    }

  ////------------------------
  ////Print Coordinate Systems
  ////------------------------
  OUT << endl << endl <<"! Coordinate Systems " ;
  OUT << endl << "! ----------------" << endl;

   vector<CoordSys*> written_cs; //to keep track of what systems have been written

  if ( NodeCoordSysList.size() != 0 ){
    typename vector<typename analysis_type::element_type*>::iterator iter; //since nodal coordsys are stored with all others
    vector<CoordSys*>::iterator csit;
    iter = Elements.begin();

    CoordSys GlobalCoordinateSystem(CoordSys::cartesian, CoordSys::deg, CoordSys::Euler312,0.,0.,0.,0.,0.,0.);

    int CSnumber;

    while ( iter != Elements.end() ){

      if ( (*iter) -> GetEleCoordSysPtr() != NULL && find(written_cs.begin(), written_cs.end(),(*iter) -> GetEleCoordSysPtr()) == written_cs.end() && *((*iter) -> GetEleCoordSysPtr()) != GlobalCoordinateSystem ){

        written_cs.push_back((*iter) -> GetEleCoordSysPtr());
        //get CoordSysNumber
        csit = find(written_cs.begin(), written_cs.end(), (*iter) -> GetEleCoordSysPtr() );
        CSnumber = 11 + distance( written_cs.begin(), csit); //local coordsys have to have numbers larger than 10


        //write the first line
        OUT << "LOCAL,R5.0,LOC,";
        OUT.setf(ios::right, ios::adjustfield);
        OUT.width(10); OUT << CSnumber << ",";
        OUT.width(3); OUT << ((*iter)->GetEleCoordSysPtr() -> GetCSType()) - 1 << ",";

        OUT.width(16); OUT << (*iter)->GetEleCoordSysPtr()->Get("x0") << ",";
        OUT.width(16); OUT << (*iter)->GetEleCoordSysPtr()->Get("y0") << ",";
        OUT.width(16); OUT << (*iter)->GetEleCoordSysPtr()->Get("z0") << endl;

        //write the second line
        OUT << "LOCAL,R5.0,ANG,";
        OUT.width(10); OUT << CSnumber << ",";
        OUT.width(3); OUT << ((*iter)->GetEleCoordSysPtr() -> GetCSType()) - 1 << ",";

        OUT.width(16); OUT << (*iter)->GetEleCoordSysPtr()->Get("Thxy") << ",";
        OUT.width(16); OUT << (*iter)->GetEleCoordSysPtr()->Get("Thyz") << ",";
        OUT.width(16); OUT << (*iter)->GetEleCoordSysPtr()->Get("Thzx") << endl;

        //write the third line
        OUT << "LOCAL,R5.0,PRM,";
        OUT.width(10); OUT << CSnumber << ",";
        OUT.width(3); OUT << ((*iter)->GetEleCoordSysPtr() -> GetCSType()) - 1 << ",";

        OUT.width(16); OUT <<  1.0  << ",";
        OUT.width(16); OUT <<  1.0  << endl;

        //write the fourth line
        OUT << "CSCIR,";
        OUT.width(10); OUT << CSnumber << ",";
        OUT.width(10); OUT << 0 << ",";
        OUT.width(10); OUT << 0 << endl;
      }
      ++iter;
    }
    OUT << "CSYS,     0" << endl;
  }
  else {
    OUT << endl << "! No special coordinate systems defined" << endl;
  }

  //// -------------------
  //// Print Nodes
  //// --------------------

  OUT << endl << endl << "! NODES" << endl;
  OUT <<         "! -----" << endl;
  OUT << "NBLOCK,6,SOLID" << endl;
  OUT << "(3i8,6e16.9)" << endl;

  // The following format specifier are necessary
  // to exactly match ANSYS conventions...
  //  mtl::dense1D<double> coords(3,0);
  Dense_Vector  coordsys(3,0), coords(3,0);
  for (unsigned i=0; i < Nodes.size(); i++){



    // 3i8
    OUT.setf(ios::right, ios::adjustfield);
    OUT << setw(8) << i+1;
    OUT << setw(8) << "0" << setw(8) << "0" << " ";

    // 6e16.9
    coords = Nodes[i].GetCoords();

    OUT.precision(9);
    OUT.setf(ios::left, ios::adjustfield);
    OUT.setf(ios::showpoint);
    OUT.width(16); OUT << coords[0];
    OUT.width(16); OUT << coords[1];
    OUT.width(16); OUT << coords[2];


    // Print nodal coordinate systems, if any
    if ( Nodes[i].GetNodeCoordSysEuler312(coordsys) ){

      OUT.width(16); OUT << coordsys[0] <<" ";
      OUT.width(16); OUT << coordsys[1] <<" ";
      OUT.width(16); OUT << coordsys[2] <<" ";
    }
    OUT << endl;
  }
  OUT << "N,R5.3,LOC, -1," << endl;
  OUT << endl;


  //// -------------------
  //// Print Elements
  //// --------------------

  // Ansys Format:
  //    Field 1 - The material number.
  //    Field 2 - The element type number.
  //    Field 3 - The real constant number.
  //    Field 4 - The section ID attribute (beam section) number.
  //              See elements BEAM188 and BEAM189 for more information.
  //    Field 5 - The element coordinate system number.
  //    Field 6 - The birth/death flag.
  //    Field 7- The solid model reference number.
  //    Field 8 - The element shape flag.
  //    Field 9 - The number of nodes per element
  //    Field 10 - The exclude key (p-elements).
  //    Field 11 - The element number.
  //    Field 12-19 - The node numbers.
  //    The next line will have the additional node numbers if there are more than eight.

  // Some variables
  int matIndex, propIndex = 0, CSnumber;
  int nIndex;
  PtrVector<Material*>::iterator matit;
  vector<PropertySet>::iterator propit;
  typename vector<typename analysis_type::node_type>::iterator nit;
  vector<CoordSys*>::iterator csit;

  OUT << endl << "! ELEMENTS" << endl;
  OUT <<         "! --------" << endl;
  OUT << "EBLOCK,19,SOLID" << endl;
  OUT << "(19i7)" << endl;

  // initialize an iterator for "distance" commands
  eleit_begin = Elements.begin();
  eleit = eleit_begin;

  typename analysis_type::element_type* TmpElePtr = NULL;

  while (eleit != Elements.end() ){             //go through all elements

    if ( (*eleit)->GetStatus() ){       // check if element is active

      Id = (*eleit)->GetId();

      //Check if the element is a triangular shell and make a temporary
      //copy of such an object. In this manner the felyx model stays unchanged
      //checking if there are triangular elements. If true the node list will be expanded
      if ( Id == 393 || Id == 391 ){
        if ( Id == 393 ){
          //TmpElePtr = StructElementFactory::Instance().CreateObject( 93 ) ;
          create_element(TmpElePtr,93);
          TmpElePtr->SetMaterialPtr((*eleit)->GetMaterialPtr());
          TmpElePtr->SetPropertiesPtr((*eleit)->GetPropertiesPtr());
          TmpElePtr->SetRefNumber((*eleit)->GetRefNumber());
          TmpElePtr->SetEleCoordSysPtr((*eleit)->GetEleCoordSysPtr());
        }
        if ( Id == 391 ){
          //TmpElePtr = StructElementFactory::Instance().CreateObject( 91 ) ;
          create_element(TmpElePtr,91);
          TmpElePtr->SetLaminatePtr((*eleit)->GetLaminatePtr());
          TmpElePtr->SetRefNumber((*eleit)->GetRefNumber());
          TmpElePtr->SetEleCoordSysPtr((*eleit)->GetEleCoordSysPtr());
        }

        //copy the nodes to the new shell
        TmpElePtr->NodeVec[0] = (*eleit)->NodeVec[0];
        TmpElePtr->NodeVec[1] = (*eleit)->NodeVec[1];
        TmpElePtr->NodeVec[2] = (*eleit)->NodeVec[2];
        TmpElePtr->NodeVec[3] = (*eleit)->NodeVec[2];
        TmpElePtr->NodeVec[4] = (*eleit)->NodeVec[3];
        TmpElePtr->NodeVec[5] = (*eleit)->NodeVec[4];
        TmpElePtr->NodeVec[6] = (*eleit)->NodeVec[2];
        TmpElePtr->NodeVec[7] = (*eleit)->NodeVec[5];
      }
      else {
        TmpElePtr = (*eleit);
      }

      // Print material type : Field 1
      // Get index for materials
      if ( TmpElePtr->GetMaterialPtr() != NULL) {
        matit = find(Materials.begin(), Materials.end(), TmpElePtr->GetMaterialPtr() );
        matIndex = distance( Materials.begin(), matit);
      }
      else
        matIndex = 0;

      OUT << setw(7) << matIndex+1;

      // Print element refernce number  : Field 2
      OUT << setw(7)  << TmpElePtr->GetRefNumber();

      // Print propertyset  : Field 3
      if ( TmpElePtr->GetPropertiesPtr() != NULL)
      {
        propit=Properties.begin();
        while ( *propit != *(TmpElePtr->GetPropertiesPtr()) && propit != Properties.end() )
        {
                ++propit;
        }
        if ( propit != Properties.end() )
                propIndex = std::distance( Properties.begin(), propit);
      }
      else if ( TmpElePtr -> GetLaminatePtr() != NULL ) {

        //std::vector<Laminate>::iterator lamit = static_cast<std::vector<Laminate>::iterator > ( TmpElePtr->GetLaminatePtr() );

        std::vector<Laminate>::iterator lamit=Laminates.begin();
        while ( *lamit != *(TmpElePtr->GetLaminatePtr()) && lamit != Laminates.end() )
        {
                ++lamit;
        }
        if ( lamit != Laminates.end())
                propIndex = std::distance( Laminates.begin(), lamit );
        propIndex += Properties.size();
      }
      else{
        propIndex = 0;
      }
      OUT << setw(7) << propIndex+1;

      // Following field is  not implemented yet:
      OUT << setw(7) << "1"; //: Field 4

      //The element coordinate system : Field 5
      if (TmpElePtr->GetEleCoordSysPtr() != NULL){
        csit = find(written_cs.begin(), written_cs.end(), TmpElePtr -> GetEleCoordSysPtr() );
        CSnumber = 11 + distance( written_cs.begin(), csit); //local coordsys have to have numbers larger than 10
      }
      else CSnumber = 0; //set to the global coordsys
      OUT << setw(7) << CSnumber;


      // Following fields are not implemented yet:
      OUT << setw(7) << "0"; //: Field 6
      OUT << setw(7) << "0"; //: Field 7
      OUT << setw(7) << "0"; //: Field 8

      // Number of nodes in element : Field 9
      OUT << setw(7) << TmpElePtr->GetNodeCount();

      // Following field is not implemented : Field 10
      OUT << setw(7) << "0";

      // Print element index : Field 11
      OUT << setw(7) << distance(eleit_begin, eleit)+1;

      // Print node indices : Fields 12 - 19
      for (unsigned n=0; n < TmpElePtr->GetNodeCount(); n++){

        if ( TmpElePtr->NodeVec.empty() ) {
          cerr << "ERROR in SaveFEData-Ansys.cc" << endl;
          cerr << "There are elements without nodes " << endl;
          exit(1);
        }

        nIndex = distance( Nodes.begin(), TmpElePtr->NodeVec[n]);

        // If there are more than 8 nodes, go to next line      for (unsigned mn = 0 ; mn < Materials.size() ; ++mn ) {
        if ( n == 8 )
          OUT << endl;

        OUT<< setw(7) <<  nIndex+1;
      }

      OUT << endl;
    }
    ++eleit;

  }
  OUT << setw(7) << "-1" << endl;
  OUT << endl;


  //// -------------------
  //// Print Materials
  //// --------------------

  // Only implemented for isotropic materials so far !!!
  // Doesn't check if material is isotrop yet !!!

  OUT << endl << "! MATERIALS" << endl;
  OUT <<         "! ---------" << endl;



  for (unsigned i=0; i < Materials.size(); ++i){ //loop over all materials

    if ( Materials[i]->ClassName() == "IsotropicMaterial" ){

      OUT << "MPTEMP,R5.0, 1, 1, 0. ," << endl;
      OUT << "MPDATA,R5.0, 1,EX," << setw(7) << i+1 ;
      OUT << ", 1, " << Materials[i]->Get("E") << endl;

      OUT << "MPTEMP,R5.0, 1, 1, 0. ," << endl;
      OUT << "MPDATA,R5.0, 1,PRXY," << setw(7) << i+1 ;
      OUT << ", 1, " << Materials[i]->Get("nu") << endl;

      OUT << "MPTEMP,R5.0, 1, 1, 0. ," << endl;
      OUT << "MPDATA,R5.0, 1,DENS," << setw(7) << i+1 ;
      OUT << ", 1, " << Materials[i]->Get("rho") << endl;
    }
    else if ( Materials[i]->ClassName() == "TransverseIsotropicMaterial23" ){

      OUT << "MPTEMP,R5.0, 1, 1, 0.0 ," << endl;
      OUT << "MPDATA,R5.0, 1,EX," << setw(7) << i+1 ;
      OUT << ", 1, " << Materials[i]->Get("E1") << endl;

      OUT << "MPTEMP,R5.0, 1, 1, 0.0 ," << endl;
      OUT << "MPDATA,R5.0, 1,EY," << setw(7) << i+1 ;
      OUT << ", 1, " << Materials[i]->Get("E2") << endl;

      OUT << "MPTEMP,R5.0, 1, 1, 0.0 ," << endl;
      OUT << "MPDATA,R5.0, 1,GXY," << setw(7) << i+1 ;
      OUT << ", 1, " << Materials[i]->Get("G12") << endl;
      
      OUT << "MPTEMP,R5.0, 1, 1, 0.0 ," << endl;
      OUT << "MPDATA,R5.0, 1,GYZ," << setw(7) << i+1 ;
      OUT << ", 1, " << Materials[i]->Get("G23") << endl;

      OUT << "MPTEMP,R5.0, 1, 1, 0.0 ," << endl;
      OUT << "MPDATA,R5.0, 1,PRXY," << setw(7) << i+1 ;
      OUT << ", 1, " << Materials[i]->Get("nu12") << endl;

      OUT << "MPTEMP,R5.0, 1, 1, 0. ," << endl;
      OUT << "MPDATA,R5.0, 1,DENS," << setw(7) << i+1 ;
      OUT << ", 1, " << Materials[i]->Get("rho") << endl;
    }
    else {
      cerr << endl << "ERROR: In IOAnsysSaveModel.cc"
           << endl << "ERROR: material doesnt exist or is not implemented yet" << endl;
      exit(1);
    }
  }
  OUT << endl;
  
  //// -------------------
  //// Print Failure Criteria values
  //// --------------------

  // Implemented for transverse isotropic materials only!!!

  OUT << endl << "! FAILURE CRITERIA" << endl;
  OUT <<         "! ----------------" << endl;
  
  // iterating over all materials while checking whether a FC exists
  for (unsigned i=0; i < Materials.size(); ++i){ //loop over all materials
    if ( Materials[i]->ClassName() == "TransverseIsotropicMaterial23" ){
      if ( Materials[i]->Get("Xt") != 0. ) {  // the fc exists
        OUT << "FC, " << i+1 << ", TEMP, , " << endl;
        OUT << "FC, " << i+1 << ",S,XTEN," << setw(7) << Materials[i]->Get("Xt") << endl;
        OUT << "FC, " << i+1 << ",S,XCMP," << setw(7) << -Materials[i]->Get("Xc") << endl;
        OUT << "FC, " << i+1 << ",S,YTEN," << setw(7) << Materials[i]->Get("Yt") << endl;
        OUT << "FC, " << i+1 << ",S,YCMP," << setw(7) << -Materials[i]->Get("Yc") << endl;
        OUT << "FC, " << i+1 << ",S,ZTEN," << setw(7) << 0.0001 << endl;
        OUT << "FC, " << i+1 << ",S,ZCMP," << setw(7) << -0.0001 << endl;
        OUT << "FC, " << i+1 << ",S,XY,"   << setw(7) << Materials[i]->Get("S")  << endl;
        OUT << "FC, " << i+1 << ",S,YZ,"   << setw(7) << 0.0001  << endl;
        OUT << "FC, " << i+1 << ",S,XZ,"   << setw(7) << 0.0001  << endl;
        OUT << "FC, " << i+1 << ",S,XYCP," << -1.0 << "," << -1.0 << "," << -1.0 << "," << -1.0 << "," << -1.0 << "," << -1.0 << endl;
        OUT << "FC, " << i+1 << ",S,YZCP," << -1.0 << "," << -1.0 << "," << -1.0 << "," << -1.0 << "," << -1.0 << "," << -1.0 << endl;
        OUT << "FC, " << i+1 << ",S,XZCP," << -1.0 << "," << -1.0 << "," << -1.0 << "," << -1.0 << "," << -1.0 << "," << -1.0 << endl;
      }
    }
    else if ( Materials[i]->ClassName() == "IsotropicMaterial" ){
      // no failure criteria for isotropic material defined yet
    }
    else {
      cerr << endl << "ERROR: In IOAnsysSaveModel.cc"
           << endl << "ERROR: material doesnt exist or is not implemented yet" << endl;
      exit(1);
    }
  }
  OUT << endl;
  



  //// -------------------------
  //// Print Boundary Conditions
  //// -------------------------

  OUT << endl << "! BOUNDARY CONDTIONS" << endl;
  OUT <<         "! ------------------" << endl;

  // loop through nodes
  for (unsigned i=0; i < Nodes.size(); i++){
    if ( Nodes[i].ExistBoundCon() ){

    print_bc(Nodes[i].BoundConPtr, OUT, i);

    }
  }

  OUT<<endl;
  OUT<<"/GO"<<endl;
  OUT<<"FINISH"<<endl;

  //// ---------------------------------------------------------------
  //// Close File and finish Load Function
  //// ---------------------------------------------------------------
  OUT.close();

  return path;

}

template<class analysis_type>
void IOAnsysFormat<analysis_type>::print_bc(StructBoundCon* BCPtr, ofstream& OUT, unsigned i ){
  double value =0;
      // check for all BC's and print'em if they are set
      if ( BCPtr->get(Dx,value) )
        OUT << "D," << setw(7) << i+1 << ",UX  ," << setw(15) << value << endl;
      if ( BCPtr->get(Dy,value) )
        OUT << "D," << setw(7) << i+1 << ",UY  ," << setw(15) << value << endl;
      if ( BCPtr->get(Dz,value) )
        OUT << "D," << setw(7) << i+1 << ",UZ  ," << setw(15) << value << endl;
      if ( BCPtr->get(Rx,value) )
        OUT << "D," << setw(7) << i+1 << ",ROTX  ," << setw(15) << value << endl;
      if ( BCPtr->get(Ry,value) )
        OUT << "D," << setw(7) << i+1 << ",ROTY  ," << setw(15) << value << endl;
      if ( BCPtr->get(Rz,value) )
        OUT << "D," << setw(7) << i+1 << ",ROTZ  ," << setw(15) << value << endl;
      if ( BCPtr->get(Fx,value) )
        OUT << "F," << setw(7) << i+1 << ",FX  ," << setw(15) << value << endl;
      if ( BCPtr->get(Fy,value) )
        OUT << "F," << setw(7) << i+1 << ",FY  ," << setw(15) << value << endl;
      if ( BCPtr->get(Fz,value) )
        OUT << "F," << setw(7) << i+1 << ",FZ  ," << setw(15) << value << endl;
      if ( BCPtr->get(Mx,value) )
        OUT << "F," << setw(7) << i+1 << ",MX  ," << setw(15) << value << endl;
      if ( BCPtr->get(My,value) )
        OUT << "F," << setw(7) << i+1 << ",MY  ," << setw(15) << value << endl;
      if ( BCPtr->get(Mz,value) )
        OUT << "F," << setw(7) << i+1 << ",MZ  ," << setw(15) << value << endl;
}

template<class analysis_type>
void IOAnsysFormat<analysis_type>::print_bc(LcmBoundCon* BCPtr, ofstream& OUT, unsigned i ){
  double value =0;
      // check for all BC's and print'em if they are set
      if ( BCPtr->get(P,value) )
        OUT << "D," << setw(7) << i+1 << ",P  ," << setw(15) << value << endl;
}
