//-----------------------------------------------------------------------------
// IOAnsysLoadModel.cc
//
// begin     : Mar 18 2001
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
   


////
//// Implementation of the Member function LoadModel of class IOAnsysFormat
////
template<class analysis_type>
string IOAnsysFormat<analysis_type>::LoadModel ( string Fname_, string DataDir_ ) {

  //// ------------------------------------------------------------------------
  //// Temporary data constructs needed through this procedural member function
  //// ------------------------------------------------------------------------
  
  
  ifstream 		IN;			// File stream
  string 		word;			// Smallest entity, which is read: a word
  int 			pos=0;			// Actual position in file
  int       fcpos=0;
  
  map<int, int> 	ElementNr2Type;		// Map from Ansys number to element type number.
  vector<int>   	TempNodeNumbers;	// Temporarly store node numbers 
  vector<CoordSys> 	TempCoordSysList;	// Temporarly store nodecoordinate systems.
  vector<typename analysis_type::bc_type> 	TempBoundaryConditions;	// Temporarily store boundary conditions systems.
  std::vector<int> 		TempBCNumbers;		// Temporarly store BC numbers
  map<int, Material*> 	MaterialNr2Ptr;		// Map for material numbers --> <number, ptr>
  map<int, int> 	PropertiesRcNr2Etype;	// Map < Real Constant Nr , Element Type Nr>.
  map<int, PropertySet*> PropertyNr2Ptr;    	// Map from real constant number to property pointer.
  map<int, CoordSys*>   CoordSysNr2Ptr;         // Map from coordsys number to coordsys pointer
  map<string, double> 	MaterialValues;		// Map for material values <key, value>.
  map<string, double>   FCValues;         // Map for failure criteria values <key, value>
  map<int, Laminate*>   RCSetNr2LaminatePtr;       // Map to assign the RCSet to a Laminate


  // All kinda dummy variables
  string varname, matid, bclabel, fcid;
  int number, oldnumber, strsize, type, matnr, oldmatnr, zero, ii, fcnr, oldfcnr;
  unsigned typenr, count;
  int materialnr, realconstantnr, coordsysnr, nodenr, nodecount;
  int RCSetNumber, etype, linewidth;


  double value, fcvalue;
  double x, y, z, thxy, thyz, thzx;

  vector<int>::iterator nodenumberit;
  PropertySet tmpProperty;

  //// -------------
  //// File handling
  //// -------------
  
  // Modify File name settings, if necessary
  SetFileNames( Fname_, DataDir_ );

  // Create path of file to be openend
  string path = GetLoadPath();

  // Open File to read
  IN.open( path.c_str(), ios::in );
  if( !IN.is_open() ){
    cout << "ERROR in IOAnsysFormat::LoadModel()" << endl;
    cout <<"Couldn't open file: "<< path << endl;
    exit(1);
  }

  //// ------------------------------------------------------------------
  //// Stepping once through input file, collecting all kinda information
  //// ------------------------------------------------------------------

  number = type = matnr = oldmatnr = oldfcnr = 0;

  // Go to begin of file
  IN.clear();
  IN.seekg( 0, ios::beg );
  
  // Loop through file
  do{
    pos = IN.tellg();	// Get position of begin of line
    getline(IN,word);	// Read one line (until delimiter '\n' is found)

    // Get detail information about nodes
    if( word.substr( 0, Details["Nodes"].Id1.size() ) == Details["Nodes"].Id1 ){
      Details["Nodes"].Startloc = pos;
      IN.ignore(300, '\n');
      
      // Loop over all nodes until end keyword of nodeblock
      do{
	getline(IN,word);
	Details["Nodes"].Count++;
      }while( word.substr( 0, Details["Nodes"].Id2.size() ) != Details["Nodes"].Id2);
      Details["Nodes"].Count--;
    }
    
    // Get numbers for element types 
    if( word.substr( 0, Details["ElementTypes"].Id1.size() ) == Details["ElementTypes"].Id1 ){
      IN.seekg(pos,ios::beg);
      IN.ignore(100,',');
      IN >> number;				// Read number of element type
      IN.ignore(100,',');
      IN >> type;			  	// Read type
      IN.seekg( pos+word.size() , ios::beg );
      
      // Write the Ansys number found into the map "ElementNr2Type"
      ElementNr2Type[number] = type;
    }
    
    // Get numbers of elements 
    if( word.substr( 0, Details["Elements"].Id1.size() ) == Details["Elements"].Id1){
      Details["Elements"].Startloc = pos;
      
      int number =0, entry =0, nodecount =0, etype =0, rcnr =0;
      
      IN.ignore(300,'\n'); 			// Ignore one line
      IN>> number;				// Read first int on next line
      do{
	entry=1;				// Entry 1 on line is read
	do{
	  IN >> number;				// Read integer
	  ++entry;				// Increment entry

	  // Find element types for appropriate real constant sets
	  if(entry == 2 ){
	    etype = ElementNr2Type[number];

	    // if this element has a property set, then
	    if ( ElementNr2RCSetPosition.find( etype ) != ElementNr2RCSetPosition.end() ){
	      IN >> rcnr;			// read real constant set number
	      entry++;
	      // if this real constant set nr occurs first time
	      if ( PropertiesRcNr2Etype.find(rcnr) == PropertiesRcNr2Etype.end() ){  
		// create new entry in PropertiesRcNr2Etype
		PropertiesRcNr2Etype[rcnr] = etype;
	      }
	    }
	  }
	  
	  if(entry == 9)			// Entry 9 --> # nodes per element
	    nodecount = number;
	  
	}while(entry < 11+nodecount );	// each element got 11 numbers + nodenumbers
	Details["Elements"].Count++;	// increment total count
		
	IN >> number;			// read first integer on next line
      }while(number != -1);		// do this until '-1' is first integer
    }
    
    // Get numbers of materials 
    if( word.substr( 0, Details["Materials"].Id1.size() ) == Details["Materials"].Id1){
      if( Details["Materials"].Count == 0 )
	Details["Materials"].Startloc = pos;
      
      IN.seekg(pos,ios::beg);		// go to begin of line
      IN.ignore(100,',');		// ignore 4 entries
      IN.ignore(100,',');
      IN.ignore(100,',');
      IN.ignore(100,',');
      IN >> matnr;			// Read Material number
      
      if( matnr != oldmatnr ){		// If Material Nr changed, increment count
	Details["Materials"].Count++;
	oldmatnr = matnr;
      }
      getline(IN,word);			// Read until end of line
    }
    
    // Get the startloc and the number of the failure criteria values
    if( word.substr( 0, Details["FailureCriteria"].Id1.size() ) == Details["FailureCriteria"].Id1){
      if( Details["FailureCriteria"].Count == 0 )
      Details["FailureCriteria"].Startloc = pos;
      
      IN.seekg(pos,ios::beg);   // go to begin of line
      IN.ignore(100,',');   // ignore 1 entry
      IN >> fcnr;           // read the fc number
      
      if( fcnr != oldfcnr ){ // if FC number is changed
        Details["FailureCriteria"].Count++;
        oldfcnr = fcnr;
        }
        getline(IN,word);
      
    }
    
    
    // Get numbers of real constants
    if( word.substr( 0, Details["RealConstants"].Id1.size() ) == Details["RealConstants"].Id1){
      Details["RealConstants"].Startloc = pos;
      IN.seekg(pos, ios::beg);		// go to begin of line
      getline(IN,word,',');		
      IN >> Details["RealConstants"].Count;		// read numbers of Real Constant sets
      IN.ignore(300,'\n');		// Read until end of line
    }


    // Get numbers for boundary conditions
    if( word.substr( 0, Details["BoundaryConditions"].Id1.size() ) 
	==  Details["BoundaryConditions"].Id1
	|| 
	word.substr( 0, Details["BoundaryConditions"].Id2.size() ) 
	== Details["BoundaryConditions"].Id2 ) {
      
      if(Details["BoundaryConditions"].Count ==0)
	Details["BoundaryConditions"].Startloc = pos;
      Details["BoundaryConditions"].Count++;

      IN.seekg(pos, ios::beg);
      IN.ignore(100,',');
      IN >> number;
      IN.ignore(200,'\n');	// read until end of line
      
      if ( std::find( TempBCNumbers.begin(), TempBCNumbers.end(), number) == TempBCNumbers.end() )
	TempBCNumbers.push_back(number);
    }

    // Get numbers of coordinate systems defined
    if( word.substr( 0, Details["CoordSystems"].Id1.size() ) == Details["CoordSystems"].Id1){
      if( Details["CoordSystems"].Count == 0 )
	Details["CoordSystems"].Startloc = pos;
      
      Details["CoordSystems"].Count++;      
    }

    // Get all kinda scalar parameters defined in ANSYS 
    // check if numerical or character parameter by testing if there is a " ' " in this line
    // only if a numerical value is found, action takes place !
    if( word.substr( 0, Details["Parameters"].Id1.size() ) 
	== Details["Parameters"].Id1 
	&& word.find("'") > word.size() ){
      IN.seekg(pos, ios::beg);
      IN.ignore(100,',');
      getline(IN,varname,',');   // get variable name
      IN.seekg(-1,ios::cur);     // one character back
      IN.ignore(100,',');        // ignore til next ','
      IN >> value;               // read value
      IN.ignore(200,'\n');	// read until end of line
      
      // remove white spaces at end of string      
      varname = varname.substr( 0, varname.find(" ",0) );
      
      AnsysParameters[varname] = value;                        
    }
    
  }while( !IN.eof() );


  //// ---------------------------------------------------------------
  //// Read nodes and nodal coordinatesystems and store'em accordingly
  //// ---------------------------------------------------------------

  // Resize data structs
  Nodes.resize(            Details["Nodes"].Count );
  TempCoordSysList.resize( Details["Nodes"].Count );
  TempNodeNumbers.resize(  Details["Nodes"].Count );

  ////
  //// a. ) Read nodes into Nodes and TempNodeNumbers
  //// 	    Read all node coordinate systems into TempCoordSysList

  IN.clear();						// Clear all kinda flags of "ifstream"
  pos = Details["Nodes"].Startloc;			// Go to begin of node block
  IN.seekg(pos, ios::beg);
  IN.ignore(300,'\n');IN.ignore(300,'\n');		// Ignore the first two lines
  
  for (unsigned n=0; n < Nodes.size(); n++){		// Loop over all nodes
    
    pos = IN.tellg();
    IN.ignore( 300 , '\n' );
    linewidth = IN.tellg();
    linewidth -= pos;			// determine linewidth
    IN.seekg(pos,ios::beg);		// go back to begin of line
    
    // Set all variables  of a node to zero
    x = y = z = thxy = thyz = thzx = 0.0;
    
    IN >> number;			// read node number
    IN >> zero;
    IN >> zero;			// read next two ints
    IN >> x;				// read first coordinate
    
    if ( linewidth > 45 ){		// check if ANSYS wrote 2nd coord
      IN >> y;			// read second coordinate
      if ( linewidth > 65 ){		// check if ANSYS wrote 3rd coord
	IN >> z;			// read third coordinate
	
	if ( linewidth > 80 ){			// check if ANSYS nodal coord sys
	  IN >> thxy;				// read 1st nodal coord angle
	  if ( linewidth > 95 ){		// check if ANSYS nodal coord sys
	    IN >> thyz;			// read 2nd  nodal coord angle
	    if ( linewidth > 110 ){		// check if ANSYS nodal coord sys
	      IN >> thzx;			// read 3rd  nodal coord angle
	    }
	  }
	}
      }
    }
    // Write values to containers
    TempNodeNumbers[n] = number;
    Nodes[n].set(x,y,z);
    TempCoordSysList[n].SetType(CoordSys::cartesian);
    TempCoordSysList[n].SetAngleType(CoordSys::deg);
    TempCoordSysList[n].Set( CoordSys::Euler312, thxy, thyz, thzx );
    
    IN.seekg(-1,ios::cur);		// go to next line
    IN.ignore(300,'\n');
  }

  ////
  //// b.) Copy TempCoordSysList (for every node a coord sys is stored)
  ////     to NodeCoordSysList (every coord sys exists only once)

  CoordSys ZeroCoordSys(CoordSys::cartesian, CoordSys::deg, CoordSys::Euler312, 0.0, 0.0, 0.0);
  vector<CoordSys>::iterator cs_it = NodeCoordSysList.begin();
  
  map<int, int> CoordSysPos;			// Map: <pos in TempCoorSysList, pos in CoordSysList >
  
  // loop through temporary coord sys list
  for (unsigned i=0; i < TempCoordSysList.size(); ++i) {
    
    // if coord sys is not equal global system
    if ( TempCoordSysList[i] != ZeroCoordSys ){
      
      // check if TempCoordSysList[i] is already in CoordSysList
      cs_it = find(NodeCoordSysList.begin(), NodeCoordSysList.end(), TempCoordSysList[i] );
      
      // If coordinate system doesn't exist yet in CoordSysList, create it
      // Remember to newly evaluate the iterators value, cause of resizing the vector !
      if ( cs_it == NodeCoordSysList.end() ){
	NodeCoordSysList.push_back( TempCoordSysList[i] );
	cs_it = NodeCoordSysList.end();
	cs_it--;
	}
      
      // put positions to map
      CoordSysPos[i] = distance( NodeCoordSysList.begin() , cs_it );      

    }
  }


  //Assign the size of the coordinate systems vector (element and nodal cs!!!)
  //prior to assigning the pointers from nodes to 
  //--------------------------------------------------------------------------
  unsigned csindex = NodeCoordSysList.size(); //this value is used later...
  NodeCoordSysList.resize( csindex + Details["CoordSystems"].Count );

  ////
  //// c.) Create pointers to the node coordinate systems and assign the csptr
  //// 
  map<int, int>::iterator pos_it = CoordSysPos.begin();
  
  // loop through map
  while ( pos_it != CoordSysPos.end() ){
    // pos_it->first gives index in "Nodes"
    // pos_it->second gives index in "CoordSysList"
    Nodes[ pos_it->first ].set( &NodeCoordSysList[ pos_it->second ] ); //sets the pointer
    pos_it++;
  }  



  
  //// ---------------------------------------------------------------
  //// Read Materials
  //// ---------------------------------------------------------------

  number=oldnumber=strsize=fcnr=oldfcnr=0;
  value=0.0;
  
  // Reserve storage for the material vector
  Materials.reserve( Details["Materials"].Count );
  
  // Go to begin of element block in file
  IN.clear();						// Clear all kinda flags of "ifstream"
  pos = Details["Materials"].Startloc;		// Go to begin of element block
  IN.seekg(pos, ios::beg);
  

  IsotropicMaterial Isotropic;
  TransverseIsotropicMaterial23 TransverseIsotropic23;
  
  for (int m=0; m < Details["Materials"].Count; m++){
    MaterialValues.clear();		// Clear map
    IN.seekg(pos, ios::beg);
    // Read material values of one material
    while (1) {
      pos = IN.tellg();			// Remember start position of line
      
      getline(IN, word, ',');		// Get first word on line
      
      // If first word on line is not equal "MPDATA", exit
      if( word.substr( 0, Details["Materials"].Id1.size() ) != Details["Materials"].Id1){
      ++number;  // increment, because of the corresponding failure criteria
	break;
      }
      
      IN.ignore(100,',');
      IN.ignore(100,',');
      getline(IN,matid,',');		// Get matid for value, example "EX"
      
      IN >> number;			// Get material number
      if (MaterialValues.size() == 0)	// first time, do something
	oldnumber = number;
      if ( number != oldnumber )	{	// upps, next material, thus get out of loop
	IN.seekg(pos, ios::beg);	// get back to begin of this material set
	break;
      }      
      
      IN.ignore(100,',');		// Get next ","
      IN.ignore(100,',');
      IN >> value;			// Get value
      IN.ignore(100,'\n');		// Ignore rest of line
      IN.ignore(100,'\n');		// Ignore next line
            
      strsize = matid.find(" ",0);
      MaterialValues[ matid.substr(0,strsize) ] = value;	// fill map, removing spaces
    }
    //WriteTransverseIsotropicMaterial23
    if ( (MaterialValues.find( "EY" ) != MaterialValues.end() && 
	 MaterialValues.find( "EZ" ) == MaterialValues.end() )
	|| (MaterialValues.find( "KYY" ) != MaterialValues.end() && 
	 MaterialValues.find( "KZZ" ) == MaterialValues.end() )
				){ //if material is transv. isotropic

      if (MaterialValues.find("EX")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "E1", MaterialValues["EX"] );
      else
	TransverseIsotropic23.Set( "E1", 0.0 );

      if (MaterialValues.find("EY")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "E2", MaterialValues["EY"] );
      else
	TransverseIsotropic23.Set( "E2", 0.0 );

      if (MaterialValues.find("GXY")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "G12", MaterialValues["GXY"] );
      else
	TransverseIsotropic23.Set( "G12", 0.0 );
	
      if (MaterialValues.find("GYZ")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "G23", MaterialValues["GYZ"] );
      else
	TransverseIsotropic23.Set( "G23", 0.0 );

      if (MaterialValues.find("PRXY")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "nu12", MaterialValues["PRXY"] );
      else
	TransverseIsotropic23.Set( "nu12", 0.0 );      
      if (MaterialValues.find("DENS")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "rho", MaterialValues["DENS"] );
      else
	TransverseIsotropic23.Set( "rho", 0.0 );

      if (MaterialValues.find("KXX")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "K1", MaterialValues["KXX"] );
      else
	TransverseIsotropic23.Set( "K1", 0.0 );

      if (MaterialValues.find("KYY")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "K2", MaterialValues["KYY"] );
      else
	TransverseIsotropic23.Set( "K2", 0.0 );

      if (MaterialValues.find("VISC")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "visc", MaterialValues["VISC"] );
      else
	TransverseIsotropic23.Set( "visc", 0.0 );
  
  /////
  // parsing for failure criteria values
  /////
  
  // Go to begin of fc block in file
  if (fcpos == 0) // initialize the first time it is run
    fcpos = Details["FailureCriteria"].Startloc;    // Go to begin of fc block
  
  IN.seekg(fcpos, ios::beg);
  
  FCValues.clear();   // Clear map
  // find a failure criteria with the same id as the current material
    while (1) {
      fcpos = IN.tellg();     // Remember start position of line
 
      getline(IN, word, ',');   // Get first word on line
      
      // If first word on line is not equal "FC", exit
      if( word.substr( 0, Details["FailureCriteria"].Id1.size() ) != Details["FailureCriteria"].Id1)
  break;
      
      IN >> fcnr;
      IN.ignore(10,',');
      
      if (FCValues.size() == 0) // first time, do something
        oldfcnr = fcnr;   
      if ( fcnr != oldfcnr )  { // upps, next fc, thus get out of loop
        IN.seekg(fcpos, ios::beg);  // get back to begin of this fc set
  break;
      }      
      
      if( number-1 == fcnr ){
        // distinguish between stress (S) and strain (E)
        getline(IN,word,',');
        //IN.ignore(10,',');
        getline(IN,fcid,','); // IN >> fcid;     // Get fcid, for example "XTEN"
          
        IN >> fcvalue;  // get the value
        IN.ignore(100,'\n');    // Ignore rest of line
        
        strsize = fcid.find(" ",0);
        FCValues[ word + fcid.substr(0,strsize) ] = fcvalue;  // fill map, removing spaces
      }
  }
  
      if (FCValues.find("SXTEN")  != FCValues.end() )
  TransverseIsotropic23.Set( "Xt", FCValues["SXTEN"] );
      else
  TransverseIsotropic23.Set( "Xt", 0.0 );
  
      if (FCValues.find("SXCMP")  != FCValues.end() )
  TransverseIsotropic23.Set( "Xc", abs(FCValues["SXCMP"]) );
      else
  TransverseIsotropic23.Set( "Xc", 0.0 );
  
      if (FCValues.find("SYTEN")  != FCValues.end() )
  TransverseIsotropic23.Set( "Yt", FCValues["SYTEN"] );
      else
  TransverseIsotropic23.Set( "Yt", 0.0 );
  
      if (FCValues.find("SYCMP")  != FCValues.end() )
  TransverseIsotropic23.Set( "Yc", abs(FCValues["SYCMP"]) );
      else
  TransverseIsotropic23.Set( "Yc", 0.0 );
  
      if (FCValues.find("SXY")  != FCValues.end() )
  TransverseIsotropic23.Set( "S", abs(FCValues["SXY"]) );
      else
  TransverseIsotropic23.Set( "S", 0.0 );
     
      //  Add material to PtrVector
      Materials.push_back( &TransverseIsotropic23 );
    }
    else if (MaterialValues.find("EY") != MaterialValues.end() &&
	     MaterialValues.find("EZ") != MaterialValues.end() ){ //if material is orthotropic
      cerr << endl << "ERROR: In IOAnsysLoadModel.cc reading in the Material!"
	   << endl << "ERROR: Ortotropic materials are not yet supported     "
	   << endl << "ERROR: Use Transverse Isotropic materials ( EZ  = 0 ) instead " << endl << endl;
    }
    else if (MaterialValues.find("KYY") != MaterialValues.end() &&
	     MaterialValues.find("KZZ") != MaterialValues.end() ){ //if material is orthotropic

      if (MaterialValues.find("KXX")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "K1", MaterialValues["KXX"] );
      else
	TransverseIsotropic23.Set( "K1", 0.0 );
      
      if (MaterialValues.find("KYY")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "K2", MaterialValues["KYY"] );
      else
	TransverseIsotropic23.Set( "K2", 0.0 );
      
      if (MaterialValues.find("VISC")  != MaterialValues.end() )
	TransverseIsotropic23.Set( "visc", MaterialValues["VISC"] );
      else
	TransverseIsotropic23.Set( "visc", 0.0 );
      Materials.push_back( &TransverseIsotropic23 );
    }
    else{ //if material is isotropic
      // Isotropic material:
      // Check if somebody used "Minor Poisson ratio"

      if (MaterialValues.find("NUXY") != MaterialValues.end() )
	MaterialValues["PRXY"] = MaterialValues["NUXY"];
      Isotropic.Set( "E", MaterialValues["EX"] );
      Isotropic.Set( "nu", MaterialValues["PRXY"] );

      // Check if density value exists
      if (MaterialValues.find("DENS") != MaterialValues.end() )
	Isotropic.Set( "rho", MaterialValues["DENS"] );
      else
	Isotropic.Set( "rho", 0 );

      // Check if conductivity value exists
      if (MaterialValues.find("KXX") != MaterialValues.end() )
	Isotropic.Set( "K", MaterialValues["KXX"] );
      else
	Isotropic.Set( "K", 0 );

      if (MaterialValues.find("VISC") != MaterialValues.end() )
	Isotropic.Set( "visc", MaterialValues["VISC"] );
      else
	Isotropic.Set( "visc", 0 );
      
      // Add material to PtrVector

      Materials.push_back( &Isotropic );
    }
    
    // Store numbers to assign pointers in elements
    MaterialNr2Ptr[ oldnumber ] = Materials.back();

  }
  
  
  
  //// ---------------------------------------------------------------
  //// Read Real Constants and store em in PropertySets
  //// ---------------------------------------------------------------

  //read once trough rcsets and collect layer data
  //this step has to be performed after the first read in 
  //since it needs informations from that step
  // Go to begin of realconst block in file
  IN.clear();					// Clear all kinda flags of "ifstream"
  pos = Details["RealConstants"].Startloc;	// Go to begin of realconstant block
  IN.seekg(pos, ios::beg);
  
  IN.ignore(300,'\n');
  IN.ignore(300,'\n');
  IN.ignore(300,'\n');

  unsigned layercount = 0;
  unsigned laminatecount = 0;
  
  // Loop over number of RC sets to get the information on the layercount
  for ( int n=0; n< Details["RealConstants"].Count; ++n ){ //loop over realconstantsets
    
    IN >> RCSetNumber;
    IN >> count;  //the number of values within the RCSet
    etype = PropertiesRcNr2Etype[RCSetNumber];
    if ( etype == 91 || etype == 391 ){
      ++laminatecount;
      IN >> value; //read the number of layers this is a double...
      layercount += (unsigned) value;
      for ( unsigned u = 0 ; u < count-1 ; ++u )
        IN >> value;
    }
  }
  Details["Layers"].Count = layercount;
  Details["Laminates"].Count = laminatecount;
  
  //reserve space for these vectors  
  Properties.reserve( Details["RealConstants"].Count );
  Layers.reserve( Details["Layers"].Count );
  Laminates.reserve( Details["Laminates"].Count );

  // Go to begin of realconst block in file
  IN.clear();					// Clear all kinda flags of "ifstream"
  pos = Details["RealConstants"].Startloc;	// Go to begin of realconstant block
  IN.seekg(pos, ios::beg);
  
  IN.ignore(300,'\n');
  IN.ignore(300,'\n');
  IN.ignore(300,'\n');
  
  // Loop over number of RC sets
  for (int n=0; n< Details["RealConstants"].Count; n++){ //loop over realconstantsets
    
    tmpProperty.clear();
    
    
    IN >> RCSetNumber;
    IN >> count;  //the number of values within the RCSet
    
    //all real constant set are read.
    etype = PropertiesRcNr2Etype[RCSetNumber];
    // If it is a RCSet for a layered shell, the real set contains a list of the layers
    // of the laminate.
    if ( etype == 91 || etype == 391 ){
      unsigned i, matnr, nl, flag = 12; //flag skipps the first 12 values of the rcset
 
      //create a laminate
      IN >> value; //read the number of layers this is a double...
      nl = (unsigned) value;
      
      Laminate tmpLaminate;
      tmpLaminate.reserve( nl );
      
      //skip 11 values. most of these values are blank, except he first and fifth
      for ( i = 0 ; i < 11 ; ++i ) IN >> value;
      
      //Loop over number of Layers
      for ( i = 0 ; i < nl  ; ++i ){
        Layer* tmpLayer;
        tmpLayer = new Layer;

        IN >> value; 
        ++flag;
        matnr = (unsigned) value;
        tmpLayer->SetMaterial( MaterialNr2Ptr[matnr] );
        IN >> value;
        ++flag;
        tmpLayer->SetAngle( value );
        IN >> value;
        ++flag;
        tmpLayer->SetThickness( value );
        
        //skip 3 thickness values in the corner nodes
        //if the maximum value of rc numbers is not yet reached
        if ( flag < count ){
          IN >> value; IN >> value; IN >> value;
          flag = flag + 3;
        }
        
        //add the layer to the vector of layers
        Layers.push_back( *tmpLayer );
        
        //add the layer to the laminate
        tmpLaminate.AddTopLayer( & Layers.back() );
      }
      
      //if the laminate is used by any element it is writen to the vectors of laminates
      if ( PropertiesRcNr2Etype.find( RCSetNumber ) != PropertiesRcNr2Etype.end() ){
        Laminates.push_back( tmpLaminate );
        RCSetNr2LaminatePtr[RCSetNumber] = &Laminates.back();
      }
      
    }
    else { //if the RCSet belongs to any element except 91 ans 391
      tmpProperty.SetId( etype );
      // Loop over number of RC values of actual set
      for (unsigned i=0; i< count ; ++i){
	IN >> value;
	// Check if the actual Realconstant i is needed
	if ( ElementNr2RCSetPosition[etype].find( i ) != ElementNr2RCSetPosition[etype].end() )
	  tmpProperty.Set( ElementNr2RCSetPosition[etype][i], value ); 
      }
      // if a RCSet is shorter than the list of Properties for the specific element
      // (etype) the properties asked for by felyx are set equal 0
      if ( count < ElementNr2RCSetPosition[etype].size() ){
	for (unsigned i = count ; i <= ElementNr2RCSetPosition[etype].size() ; i++){
	  if ( ElementNr2RCSetPosition[etype].find( i ) != ElementNr2RCSetPosition[etype].end() )
	    tmpProperty.Set( ElementNr2RCSetPosition[etype][i], 0 );	    
	}
      }
      if ( PropertiesRcNr2Etype.find( RCSetNumber ) != PropertiesRcNr2Etype.end() && etype != 0 ){
	Properties.push_back( tmpProperty  );
	// Store numbers to assign pointers in elements
	PropertyNr2Ptr[RCSetNumber] = & Properties.back();
      }
    }
  }




  //// ---------------------------------------------------------------
  //// Read Boundary Conditions
  //// ---------------------------------------------------------------

  // Reserve storage on STL vector
  TempBoundaryConditions.resize( TempBCNumbers.size() );
  
  // Go to begin of the BC block in file
  IN.clear();					// Clear all kinda flags of "ifstream"
  pos = Details["BoundaryConditions"].Startloc;	// Go to begin of element block
  IN.seekg(pos, ios::beg);
  
  ii=number=strsize=0;
  value=0.0;
  
  // Loop through all lines with BC's
  for (int i=0; i < Details["BoundaryConditions"].Count; i++){
    
    IN.ignore(100, ',');
    IN >> number;
    
    // Find  according number in TempBCNumbers
    vector<int>::iterator it;
    it = std::find(TempBCNumbers.begin(), TempBCNumbers.end(), number);
    FELYX_RUNTIME_ASSERT(it!=TempBCNumbers.end(),"IOAnsysLoadModel: Out of range access for TempBCNumbers" );
    ii = std::distance(TempBCNumbers.begin(), it);
    
    IN.ignore(100,',');
    getline(IN,bclabel,',');
    IN >> value;

    // Remove whitespaces
    strsize = bclabel.find(" ",0);
    bclabel = bclabel.substr(0,strsize);
    
    // Add value at appropriate position in BoundaryConditions
    TempBoundaryConditions[ii].set( BClabel2BCtype[bclabel] ,value);
    
    IN.ignore(200,'\n');	// go to end of line
    
  }


  // Copy TempBoundaryConditions to BoundaryConditions
  // Kick out all duplicates in BoundaryConditions
  typename vector<typename analysis_type::bc_type>::iterator bcit;
  BoundaryConditions = TempBoundaryConditions;
  
  // Loop over all BC sets
  for (unsigned i=0; i < BoundaryConditions.size(); i++){
    // Loop over all BC sets, beginning from i+1
    bcit = BoundaryConditions.begin()+i+1;
    while ( bcit != BoundaryConditions.end() ){
      // If there are identical sets, then...
      if ( *bcit == BoundaryConditions[i] ){
	bcit = BoundaryConditions.erase(bcit);
      }
      else
	bcit++;
    }
  }

  // Set appropriate pointer in "Nodes" :
  // for each node do:
  // - look up index of its boundary conditions set (if any) in TempBCNumbers
  // - tmpset = TempBoundaryCondtions[ index ] --> appropriate bc set
  // - find this set in the BoundaryConditions vector (where duplicates are removed)
  // - set the BoundaryCondition Ptr in Node
  vector<int>::iterator nrit;
  int bc;
  typename analysis_type::bc_type tmpset;
  
  for (unsigned n=0; n < Nodes.size(); n++){

    nrit = find ( TempBCNumbers.begin(), TempBCNumbers.end(), TempNodeNumbers[n] );
    if ( nrit != TempBCNumbers.end() ) {
      bc = distance(TempBCNumbers.begin(), nrit);
      
      tmpset = TempBoundaryConditions[bc];

      bcit = find( BoundaryConditions.begin(), BoundaryConditions.end(), tmpset);
      
      Nodes[n].set ( &(*bcit) );
    }
  }


  //// ---------------------------------------------------------------
  //// Read Coordinate Systems
  //// ---------------------------------------------------------------

  // Reserve storage on STL vector
  //  unsigned csindex = NodeCoordSysList.size();
  //  NodeCoordSysList.resize( csindex + Details["CoordSystems"].Count );
  
  // Go to begin of the Coordinate Systems block in file
  IN.clear();					// Clear all kinda flags of "ifstream"
  pos = Details["CoordSystems"].Startloc;	// Go to begin of element block
  IN.seekg(pos, ios::beg);
  
  number = 0;
  value = 0.0;

  // Loop through lines of a specific coord system
  for (int i=0; i < Details["CoordSystems"].Count; i++){

    //Read the first line of one coordinate system
    //--------------------------------------------
    IN.ignore(100,',');    
    IN.ignore(100,',');    
    IN.ignore(100,',');    // position is now after LOC,

    IN >> number;         //read the coordsys number

    CoordSysNr2Ptr[number] = &NodeCoordSysList[csindex];

    IN.ignore(100,',');
    IN >> number;          //read the coordsystype number

    if ( number == 0 ) NodeCoordSysList[csindex].SetType(CoordSys::cartesian);
    else if ( number == 1 ) NodeCoordSysList[csindex].SetType(CoordSys::cylindrical);
    else {
      cerr << endl << "ERROR in IOAnsysLoadModel.cc - Read Coordinate Systems " 
	   << endl << "The Coordinate Type " << number << " from Ansys is not "
	   << endl << "implemented or supported " << endl << endl;
      exit(1);
    }

    //since the coordinate systems in the ansys db allways use deg....
    NodeCoordSysList[csindex].SetAngleType( CoordSys::deg ); 

    IN.ignore(100, ',');
    IN >> value;          //read the x coordinate of the base point
    NodeCoordSysList[csindex].Set("x0", value);

    IN.ignore(100, ',');
    IN >> value;          //read the x coordinate of the base point
    NodeCoordSysList[csindex].Set("y0", value);

    IN.ignore(100, ',');
    IN >> value;          //read the x coordinate of the base point
    NodeCoordSysList[csindex].Set("z0", value);

    IN.ignore(200,'\n');	// go to end of line

    //Read the second line of one coordinate system 
    //---------------------------------------------
    IN.ignore(100,',');    
    IN.ignore(100,',');    
    IN.ignore(100,',');    
    IN.ignore(100,',');    
    IN.ignore(100,',');        //position in front of first angle

    IN >> value;          //read the first Euler312 angle
    NodeCoordSysList[csindex].Set("Thxy", value);

    IN.ignore(100, ',');
    IN >> value;          //read the second Euler312 angle
    NodeCoordSysList[csindex].Set("Thyz", value);

    IN.ignore(100, ',');
    IN >> value;           //read the third Euler312 angle
    NodeCoordSysList[csindex].Set("Thzx", value);

    IN.ignore(200,'\n');	// go to end of line

    IN.ignore(300,'\n');
    IN.ignore(300,'\n');        // ignore 2 lines

    ++csindex;
  }

  


  //// ---------------------------------------------------------------
  //// Read elements, replace node numbers by pointers and store'em 
  //// ---------------------------------------------------------------
 
  // Syntax of ANSYS EBlOCK command: 
  // -------------------------------
  //	Field 1 - The material number.
  //	Field 2 - The element type number.
  //	Field 3 - The real constant number.
  //	Field 4 - The section ID attribute (beam section) number.
  //		  See elements BEAM188 and BEAM189 for more information.
  //	Field 5 - The element coordinate system number.
  //	Field 6 - The birth/death flag.
  //	Field 7- The solid model reference number.
  //	Field 8 - The element shape flag.
  //	Field 9 - The number of nodes of this element 
  //	Field 10 - The exclude key (p-elements).
  //	Field 11 - The element number.
  //	Field 12-19 - The node numbers.
  //	The next line will have the additional node numbers if there are more than eight. 
  
  // Go to begin of element block in file
  IN.clear();						// Clear all kinda flags of "ifstream"
  pos = Details["Elements"].Startloc;		// Go to begin of element block
  IN.seekg(pos, ios::beg);
  IN.ignore(300,'\n');IN.ignore(300,'\n');		// Ignore the first two lines
  
  // Resize element vector
  Elements.resize( Details["Elements"].Count );

  // Loop over all elements
  typename PtrVector<typename analysis_type::element_type*>::iterator eleit = Elements.begin();

  while ( eleit != Elements.end() ) {
  
    IN >> materialnr;		// field 1
    IN >> typenr;		// field 2
    IN >> realconstantnr;	// field 3
    IN >> zero;			// field 4
    IN >> coordsysnr;	        // field 5
    IN >> zero;			// field 6
    IN >> zero;			// field 7
    IN >> zero;			// field 8
    IN >> nodecount;		// field 9
    				// -> nodes of this element
    				// -> can be different from eleptr->GetNodeCount() (Beam4...)
    IN >> zero;			// field 10
    IN >> zero;			// field 11
    
    // Create appropriate element using the StructElementFactory
    create_element(eleit, ElementNr2Type[ typenr ]);


    // Set element reference number
    (*eleit)->SetRefNumber( typenr );
  
    // Set material pointer or laminate pointer
    if ( dynamic_cast<LayeredShell*>( *eleit ) != NULL )
      (*eleit)->SetLaminatePtr( RCSetNr2LaminatePtr[realconstantnr] );
    else
      (*eleit)->SetMaterialPtr( MaterialNr2Ptr[ materialnr ] );

    // Set real constant pointer only if the element is derived from 
    // BaseBeam or SingleLayerShell
    if ( dynamic_cast<BaseBeam*>( *eleit ) != NULL || 
	 dynamic_cast<SingleLayerShell*>( *eleit ) != NULL || dynamic_cast<LcmElement*>( *eleit ) != NULL  ) 
      (*eleit)->SetPropertiesPtr( PropertyNr2Ptr[ realconstantnr ] );
    
    // Set coordinate system pointer only if the element is derived
    // from BaseShell 
    if ( dynamic_cast<BaseShell*>( *eleit ) != NULL || dynamic_cast<LcmElement*>( *eleit ) != NULL )
      (*eleit)->SetEleCoordSysPtr( CoordSysNr2Ptr[ coordsysnr ] );
    // Reading nodes: field 12 - 12 + number of nodes
    for (unsigned i=0; i < (*eleit)->GetNodeCount()  ; i++){
      // if ((*eleit)->GetId() == 70 && i == (*eleit)->GetNodeCount()-1 ) IN >> nodenr;
      IN >> nodenr;
      nodenumberit = find( TempNodeNumbers.begin(), TempNodeNumbers.end(), nodenr);
      (*eleit)->SetNodeIter(i , Nodes.begin()+(distance(TempNodeNumbers.begin(), nodenumberit )) );
    }
    
    // Go to next line
    IN.seekg(-1,ios::cur);
    IN.ignore(300,'\n');

    // increment eleit
    eleit++;
  }
  //// ---------------------------------------------------------------
  //// Close File and finish Load Function
  //// ---------------------------------------------------------------
  IN.close();
  
  //If there are shells in the model, check for triangular shells
  //and make new felyx elements
  //  if ( ElementNr2Type has 93 or 91 )
  Quad2Tri();
  Brick2Tet();
  
  return path;
  
}


template<class analysis_type>
void IOAnsysFormat<analysis_type>::Quad2Tri(){

  typename vector<typename analysis_type::element_type*>::iterator eleit=Elements.begin();
  typename analysis_type::element_type* myShellPtr = NULL;
  unsigned myId = 0;
  while ( eleit != Elements.end() ){
    myId = (*eleit)->GetId();
    myShellPtr = NULL;
    if ( myId == 93 || myId == 91 ) {
      //if the shell is triangular
      if ( (*eleit)->NodeVec[2] == (*eleit)->NodeVec[3] &&
	   (*eleit)->NodeVec[2] == (*eleit)->NodeVec[6] ) {
	if ( myId == 93 ) {
	  //  myShellPtr = StructElementFactory::Instance().CreateObject( 393 ) ;
	  create_element(myShellPtr,393);	
	  myShellPtr->SetMaterialPtr((*eleit)->GetMaterialPtr());
	  myShellPtr->SetPropertiesPtr((*eleit)->GetPropertiesPtr());
	  myShellPtr->SetRefNumber((*eleit)->GetRefNumber());
	  myShellPtr->SetEleCoordSysPtr((*eleit)->GetEleCoordSysPtr());
	}
	if ( myId == 91 ) {
	  //myShellPtr = StructElementFactory::Instance().CreateObject( 391 ) ;
	  create_element(myShellPtr,391);	
	  myShellPtr->SetLaminatePtr((*eleit)->GetLaminatePtr());
	  myShellPtr->SetRefNumber((*eleit)->GetRefNumber());
	  myShellPtr->SetEleCoordSysPtr((*eleit)->GetEleCoordSysPtr());
	}
	//Copy the nodes to the new shell
	myShellPtr->NodeVec[0] = (*eleit)->NodeVec[0];
	myShellPtr->NodeVec[1] = (*eleit)->NodeVec[1];
	myShellPtr->NodeVec[2] = (*eleit)->NodeVec[2];
	myShellPtr->NodeVec[3] = (*eleit)->NodeVec[4];
	myShellPtr->NodeVec[4] = (*eleit)->NodeVec[5];
	myShellPtr->NodeVec[5] = (*eleit)->NodeVec[7];
	
	(*eleit) = myShellPtr;
      }
/*       else if ( (*eleit)->NodeVec[2] == (*eleit)->NodeVec[3] ) { */
/* 	if ( myId == 55 ) { */
/* 	  create_element(myShellPtr,355); */
/* 	  myShellPtr->SetMaterialPtr((*eleit)->GetMaterialPtr()); */
/* 	  myShellPtr->SetPropertiesPtr((*eleit)->GetPropertiesPtr()); */
/* 	  myShellPtr->SetRefNumber((*eleit)->GetRefNumber()); */
/* 	  myShellPtr->SetEleCoordSysPtr((*eleit)->GetEleCoordSysPtr()); */
/* 	} */
/* 	if ( myId == 57 ) { */
/* 	  create_element(myShellPtr,357); */
/* 	  myShellPtr->SetMaterialPtr((*eleit)->GetMaterialPtr()); */
/* 	  myShellPtr->SetPropertiesPtr((*eleit)->GetPropertiesPtr()); */
/* 	  myShellPtr->SetRefNumber((*eleit)->GetRefNumber()); */
/* 	  myShellPtr->SetEleCoordSysPtr((*eleit)->GetEleCoordSysPtr()); */
/* 	} */
/* 	//Copy the nodes to the new lcm node */
/* 	myShellPtr->NodeVec[0] = (*eleit)->NodeVec[0]; */
/* 	myShellPtr->NodeVec[1] = (*eleit)->NodeVec[1]; */
/* 	myShellPtr->NodeVec[2] = (*eleit)->NodeVec[2]; */
/* 	(*eleit) = myShellPtr; */
/*       } */
    }
    ++eleit;
  } //end of loop over elements
}

template<class analysis_type>
void IOAnsysFormat<analysis_type>::Brick2Tet(){
  
  
}

