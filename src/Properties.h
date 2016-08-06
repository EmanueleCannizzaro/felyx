//-----------------------------------------------------------------------------
// Properties.h
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
   

#ifndef Properties_h
#define Properties_h Properties_h

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "Material.h"

using namespace std;
namespace fe_base{	//begin of namespace fe_base
  

  //! The class PropertySet is a container
  /*! Propertyset contains a pointer to material, an Id and a 
    map<string, double>. There are a GetDouble(string), a GetId() and a GetInt(String)
    implemented. The Get(string) overloading doesnt work since these functions
    vary only in the returntype and not in the argument type.  The Set(...) 
    is overloaded. If ever a new overloaded Set(...) is implemented the function  must 
    check for the Id string. GetSomeType(string) checks,
    if the desired key  exists. If maps to other types than Double  are needed, they have
    to be implemented by the user. A map to int is included as an example but is not 
    used at the moment
    The keys are chosen as strings instead of 
    enums to keep the most flexibility and to avoid the dissemination of
    code through FELyX...
  */
  
  class PropertySet{
  public:
    //Constructors
    //------------
    PropertySet() : MaterialPtr(NULL) { PropertyId = -1; } 
    PropertySet(Material* Material_, int Id_): MaterialPtr(Material_), PropertyId(Id_) {} 
    PropertySet(int Id_): MaterialPtr(NULL), PropertyId(Id_) {}
    ~PropertySet() {}
    //CopyConstructor
    PropertySet( const PropertySet& );

    //Equal operator
    PropertySet operator=( const PropertySet& );
    
    //Set members
    //-----------
    void SetMaterialPtr( Material* Material_ ) { MaterialPtr = Material_; };
    void SetId(int Id_) { PropertyId = Id_; };
    void Set(string, double); //Set DoubleMap
    //void Set(string, int); //Set IntMap
    
    //Get members
    //-----------
    int     GetId() { return PropertyId; };
    double  GetDouble(string); //Get the instances of DoubleMap and the Id
    map<string, double> GetDoubleMap() { return DoubleMap; }; //Get the whole DoubleMap
    //int     GetInt(string); //Get the instances of DoubleMap and the Id

    //Other members
    //-------------  
    bool operator==(const PropertySet&);
    bool operator!=(const PropertySet&);
    void CheckForProperty(map<string, double> Map_, string key_);
    //! Clear erases all contents of property set
    void clear();
    
    friend ostream& operator<<(ostream&, const PropertySet&);

    //Data members
    //------------
    Material* MaterialPtr;

  private:
    int                   PropertyId;
    map<string, double>   DoubleMap;
   //map<string, int>      IntMap;
  };

}	// end of namespace fe_base

#endif
