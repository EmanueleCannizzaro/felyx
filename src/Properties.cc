//-----------------------------------------------------------------------------
// Properties.cc
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
   
#include "Properties.h"

using namespace fe_base;

void PropertySet::CheckForProperty(map<string, double> MyMap, string key_){
  map<string, double>::iterator MapIt;
  MapIt = MyMap.find(key_);
  if (MapIt == MyMap.end()){
    cerr << "ERROR: In PropertySet::CheckForProperty(map, string)" << endl;
    cerr << "ERROR: A PropertySet map has been asked for a non-existing property" << endl;
    cerr << "ERROR: The Property falsely looked for is named: " << key_ << endl << endl;
    exit(1);
  }
}


void PropertySet::Set(string key_, double key_value) {
    DoubleMap[key_] = key_value;
}

double PropertySet::GetDouble(string key_){
  CheckForProperty(DoubleMap, key_);
  return DoubleMap[key_];
}

void PropertySet::clear(){
  DoubleMap.clear();
  PropertyId = -1;
}

//Copyconstructor
PropertySet::PropertySet( const PropertySet& original)
{
  MaterialPtr = original.MaterialPtr;
  PropertyId = original.PropertyId;
  DoubleMap = original.DoubleMap;
}

//Equal operator
PropertySet PropertySet::operator=( const PropertySet& original)
{
  MaterialPtr = original.MaterialPtr;
  PropertyId = original.PropertyId;
  DoubleMap = original.DoubleMap;
  return *this;
}

//These functions are not used yet
//--------------------------------
/*
void PropertySet::Set(string key_, int key_value) {
    IntMap[key_] = key_value;
};

int PropertySet::GetInt(string key_){
  CheckForProperty(DoubleMap, key_);
    return IntMap[key_];
}
*/



//////////////////////////////////////
/// functions not member of any class
//////////////////////////////////////
ostream& fe_base::operator<<( ostream& stream_, const PropertySet& PropertySet_ ) {
  stream_.precision(9);

  string theId = "Id";
  stream_ << "  ";
  stream_.width( 15- theId.size() );
  stream_ << theId << " : " << PropertySet_.PropertyId << endl;

  map<string, double>::const_iterator pos;
  for (pos = PropertySet_.DoubleMap.begin() ; pos != PropertySet_.DoubleMap.end() ; ++pos){
    stream_ << "  ";
    stream_.width( 15 - pos-> first.size() );
    stream_ <<  pos -> first << " : " << pos -> second << endl;
  }
  stream_.unsetf(ios::left);

 return stream_;
}

bool PropertySet::operator==(const PropertySet& a)
{
	return ( PropertyId == a.PropertyId && DoubleMap == a.DoubleMap );
}

bool PropertySet::operator!=(const PropertySet& a)
{
	return !operator==(a);
}
