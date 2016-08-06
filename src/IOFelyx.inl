//-----------------------------------------------------------------------------
// IOFelyx.inl
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
//// Constructors
////
template<class analysis_type>
IOFelyx<analysis_type>::IOFelyx(  vector<typename analysis_type::node_type>& Nodes_,
                                  vector<CoordSys>& NodeCoordSysList_,
                                  PtrVector<typename analysis_type::element_type*>& Elements_,
                                  PtrVector<Material*>& Materials_,
                                  vector<PropertySet>& Properties_,
                                  vector<typename analysis_type::bc_type>& BoundaryConditions_,
                                  vector<Layer>&  Layers_,
                                  vector<Laminate>& Laminates_ )
  : Nodes( Nodes_ ),
    NodeCoordSysList( NodeCoordSysList_ ),
    Elements( Elements_ ),
    Materials( Materials_ ),
    Properties( Properties_ ),
    BoundaryConditions( BoundaryConditions_ ),
    Layers( Layers_ ),
    Laminates( Laminates_)
{}

////
//// Specify some sort of default behavior if virtual functions do not get overloaded
////
template<class analysis_type>
string IOFelyx<analysis_type>::LoadModel (string Fname_, string DataDir_ ) {
  cerr << endl;
  cerr << "\t WARNING in IOFelyx::LoadModel" << endl;
  cerr << "\t The virtual member function LoadModel is not overloaded " << endl;
  cerr << "\t in the calling IO class; no action takes place !!" << endl;
  cerr << endl;
  return "WARNING, NO FILE LOADED!" ;
}
template<class analysis_type>
string IOFelyx<analysis_type>::SaveModel (string Fname_, string DataDir_ ) {
  cerr << endl;
  cerr << "\t WARNING in IOFelyx::SaveModel" << endl;
  cerr << "\t The virtual member function SaveModel is not overloaded " << endl;
  cerr << "\t in the calling IO class; no action takes place !!" << endl;
  cerr << endl;
  return "WARNING, NO FILE CREATED!" ;
}
template<class analysis_type>
string IOFelyx<analysis_type>::SaveResults (string Fname_, string DataDir_ ) {
  cerr << endl;
  cerr << "\t WARNING in IOFelyx::SaveResults" << endl;
  cerr << "\t The virtual member function SaveResults is not overloaded " << endl;
  cerr << "\t in the calling IO class; no action takes place !!" << endl;
  cerr << endl;
  return "WARNING, NO FILE CREATED!" ;
}


////
//// Implement non-virtual member functions of this class
////

//! Member functions to modify filename strings.
/*! If arguments are not filled in, no action is taken ! */
template<class analysis_type>
void IOFelyx<analysis_type>::SetFileNames ( string Fname_,
                             string DataDir_,
                             string LoadSuffix_,
                             string SaveSuffix_,
                             string ResSuffix_ )
{
  // If result suffix is not empty (= default argument), then store it
  if ( !ResSuffix_.empty() )
    ResSuffix = ResSuffix_;

  // If save suffix is not empty (= default argument), then store it
  if ( !SaveSuffix_.empty() )
    SaveSuffix = SaveSuffix_;

  // If load suffix is not empty (= default argument), then store it
  if ( !LoadSuffix_.empty() )
    LoadSuffix = LoadSuffix_;

  // If DataDir_ is not empty (= default argument), then store it
  if ( !DataDir_.empty() )
    DataDir = DataDir_;


  // Fname handling
  if (!Fname_.empty() ){

    // Check if filename includes suffix
    string::size_type idx = Fname_.find_last_of(".");

    if ( idx < Fname_.size() ){
      Fname = Fname_.substr(0,idx);
      LoadSuffix = Fname_.substr( idx );
    }
    else
      Fname = Fname_;


  }

}

template<class analysis_type>
std::string IOFelyx<analysis_type>::ParserFunction( std::string path_ )
{
  std::string info(" not yet implemented ");
  std::cout << info << std::endl;
  return info;

}





