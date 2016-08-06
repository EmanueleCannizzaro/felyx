//-----------------------------------------------------------------------------
// IOAnsysFormatSpecifications.h
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
   
#ifndef IOAnsysFormatSpecifications_h
#define IOAnsysFormatSpecifications_h IOAnsysFormatSpecifications_h

#include <string>
#include <map>

#include "StructBoundCon.h" // Needed, cause enum BCtype is used
#include "LcmBoundCon.h"


using namespace std;
using namespace bc;

namespace fe_base{
  
  //! Class to put in all kinda ANSYS format specificiations.
  /*! The IOAnsysFormat class is going to derive from this data container. */
  template<class analysis_type> class IOAnsysFormatSpecifications{
    
  protected:

    ////
    //// Definition of Data structs in this container
    ////
    
    //! Class EntityDetails, storing all kinda Detail information of different entities.
    /*! Storing Id's of ANSYS, location in the file where these ID's occur and counts
        of entities.
    */
    class EntityDetails{
    public:
      // Constructors
      EntityDetails() 
	: Id1("\0"), Id2("\0"), Startloc(0), Count(0) {}
      EntityDetails( string id1_, string id2_ = "\0" ) 
	: Startloc(0), Count(0)
	{ set( id1_, id2_ ); }
      
      //! Set function
      void set( string id1_, string id2_ = "\0" )
	{ Id1 = id1_; Id2 = id2_; }
      
      // Data members
      string Id1;	//! First id for an entity
      string Id2;      	//! Second id for an entity
      int Startloc;    	//! Start location of entity in data file
      int Count;       	//! Count of entities
    };
    
    
    //! Map to store all those entities, key defines the actual entity.
    map<string, EntityDetails> Details;
    
    //! Map that stores information about Real Constant Sets in Ansys.
    /*! Key: Elemnttype of the actual  RealConstant set.
        Inner Map: maps position of item in RealConstant set to name of this item
    */
    map<unsigned, map <unsigned, string > > ElementNr2RCSetPosition;	
    
    //! Map that connects ANSYS BC Labels with enum type "BCtype"
    map<string, typename analysis_type::bc_type::enumType> BClabel2BCtype;		

    //! Define a virtual Destructor
    virtual ~IOAnsysFormatSpecifications() {};
    
    ////
    //// Filling up all the data structs using the default constructor
    ////
    
    // Filling all those maps by default constructor:
    IOAnsysFormatSpecifications() {
      
      // Specify keywords for different entities
      Details["CoordSystems"]       = EntityDetails("LOCAL,R5.0,LOC,","CSYS");
      Details["Nodes"]              = EntityDetails("NBLOCK","N,R5.3");
      Details["ElementTypes"]       = EntityDetails("ET");
      Details["Elements"]           = EntityDetails("EBLOCK","-1");
      Details["Materials"]          = EntityDetails("MPDATA");
      Details["FailureCriteria"]    = EntityDetails("FC");
      Details["RealConstants"]      = EntityDetails("RLBLOCK");
      Details["BoundaryConditions"] = EntityDetails("D,", "F,");
      Details["Parameters"]         = EntityDetails("*SET");
      Details["Layers"]             = EntityDetails("RLBLOCK");
      Details["Laminates"]          = EntityDetails("RLBLOCK");
      
      // Fill ElementNr2RCSetPosition map
      ElementNr2RCSetPosition[91][0] = "LayerCount"; //not shure if this is necessary (nz)

      ElementNr2RCSetPosition[93][0] = "Thickness";
      ElementNr2RCSetPosition[93][4] = "Theta";
      
      ElementNr2RCSetPosition[1][0] = "Area";
      
      ElementNr2RCSetPosition[8][0] = "Area";
  
      ElementNr2RCSetPosition[3][0] = "Area";
      ElementNr2RCSetPosition[3][1] = "Izz";
      ElementNr2RCSetPosition[3][2] = "Height";
      ElementNr2RCSetPosition[3][3] = "ShearZ";
      
      ElementNr2RCSetPosition[4][0] = "Area";
      ElementNr2RCSetPosition[4][1] = "Izz";
      ElementNr2RCSetPosition[4][2] = "Iyy";
      ElementNr2RCSetPosition[4][3] = "ThicknessZ";
      ElementNr2RCSetPosition[4][4] = "ThicknessY";
      ElementNr2RCSetPosition[4][5] = "Theta";
      ElementNr2RCSetPosition[4][7] = "Ip";
      ElementNr2RCSetPosition[4][8] = "ShearZ";
      ElementNr2RCSetPosition[4][9] = "ShearY";

      ElementNr2RCSetPosition[57][0] = "Thickness";
      
      labelmap(BClabel2BCtype);
    }


    void labelmap(map<string, bc::StructBCtype>){

      // Fill BClabel2BCtype map
      BClabel2BCtype["UX"]   = Dx;
      BClabel2BCtype["UY"]   = Dy;
      BClabel2BCtype["UZ"]   = Dz;
      BClabel2BCtype["ROTX"] = Rx;
      BClabel2BCtype["ROTY"] = Ry;
      BClabel2BCtype["ROTZ"] = Rz;
      BClabel2BCtype["FX"]   = Fx;
      BClabel2BCtype["FY"]   = Fy;
      BClabel2BCtype["FZ"]   = Fz;
      BClabel2BCtype["MX"]   = Mx;
      BClabel2BCtype["MY"]   = My;
      BClabel2BCtype["MZ"]   = Mz;

    }


    void labelmap(map<string, bc::LcmBCtype>){

      BClabel2BCtype["TEMP"]   = P;
      BClabel2BCtype["HEAT"]   = V;

    }

  }; // of class
  
} // of namespace
#endif

