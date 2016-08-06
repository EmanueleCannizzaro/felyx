//-----------------------------------------------------------------------------
// IOFelyx.h
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

#ifndef IOFelyx_h
#define IOFelyx_h IOFelyx_h

#include <string>
#include <fstream>
#include <vector>
#include <map>

#include "PtrVector.h"
#include "ElementHeaders.h"

using namespace std;
namespace fe_base{

  //! Abstract Base Class for storing FE-Data.
  /*! Class defining interfaces to store FE-Data in files.
      Class should be derived to implement Storage functions for different formats
  */
  template<class analysis_type> class IOFelyx{

  public:

    //! Constructor, initializing references to all FE-Data objects
    IOFelyx( vector<typename analysis_type::node_type>&,
              vector<CoordSys>&,
              PtrVector<typename analysis_type::element_type*>&,
              PtrVector<Material*>&,
              vector<PropertySet>&,
              vector<typename analysis_type::bc_type>&,
              vector<Layer>&,
              vector<Laminate>& );

    //! Virtual destructor
    virtual ~IOFelyx() {};

    //! Load FE Model from file
    /*! Return type:        Path of file loaded
        Optional Arguments: Filename, DataDirectory */
    virtual string LoadModel( string Fname_="", string DataDir_="");
    virtual std::string ParserFunction(std::string path_);


    //! Save FE Model to file.
    /*! Return type:        Path of file loaded
        Optional Arguments: Filename, DataDirectory */
    virtual string SaveModel( string Fname_="", string DataDir_="" );

    //! Save FE-Results to file.
    /*! Return type:        Path of file loaded
        Optional Arguments: Filename, DataDirectory */
    virtual string SaveResults( string Fname_="", string DataDir_="" );

    //! Return load path stored
    string GetLoadPath() { return DataDir + "/" + Fname + LoadSuffix; }
    string GetSavePath() { return DataDir + "/" + Fname + SaveSuffix; }
    string GetResPath() { return DataDir + "/" + Fname + ResSuffix; }

  protected:

    //! Member function to assign all kinda filenames, directories, suffixes
    /*! Arguments: Fname, Datadir, LoadSuffix, SaveSuffix, ResSuffix */
    void SetFileNames ( string Fname_ = "",
                        string DataDir_ = "",
                        string LoadSuffix_ = "",
                        string SaveSuffix_ = "",
                        string ResSuffix_  = "" );

    // References to FE-data constructs
    vector<typename analysis_type::node_type>&         Nodes;
    vector<CoordSys>&                                  NodeCoordSysList;
    PtrVector<typename analysis_type::element_type*>&  Elements;
    PtrVector<Material*>&                              Materials;
    vector<PropertySet>&                               Properties;
    vector<typename analysis_type::bc_type>&           BoundaryConditions;
    vector<Layer>&                                     Layers;
    vector<Laminate>&                                  Laminates;

  private:
    // Strings to store filenames / directories / suffixes
    string DataDir;
    string Fname;
    string LoadSuffix;
    string SaveSuffix;
    string ResSuffix;

  };

#include "IOFelyx.inl"

void create_element(PtrVector<StructElement*>::iterator eleit, int id_);
void create_element(PtrVector<LcmElement*>::iterator eleit, int id_);
void create_element(StructElement*& eleit, int id_);
void create_element(LcmElement*& eleit, int id_);

} // of namespace
#endif

