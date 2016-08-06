//-----------------------------------------------------------------------------
// IOAnsysFormat.h
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

#ifndef IOAnsysFormat_h
#define IOAnsysFormat_h IOAnsysFormat_h

#include <string>
#include <iomanip>
#include "IOFelyx.h"
#include "IOAnsysFormatSpecifications.h"

using namespace std;

namespace fe_base{

  //! Class to store FE-Data in ANSYS-Format
  template<class analysis_type> class IOAnsysFormat : public IOFelyx<analysis_type>, public IOAnsysFormatSpecifications<analysis_type>
  {

  public:

    //! Constructor, initializing references to all FE-Data objects.
    /*! Additionally, filename and datadirectory can optionally be specified */
    IOAnsysFormat( vector<typename analysis_type::node_type>& Nodes_,
                   vector<CoordSys>& NodeCoordSysList_,
                   PtrVector<typename analysis_type::element_type*>&   Elements_,
                   PtrVector<Material*>&        _Materials,
                   vector<PropertySet>&         _Properties,
                   vector<typename analysis_type::bc_type>&     _BoundaryConditions,
                   vector<Layer>&               Layers_,
                   vector<Laminate>&            Lamiantes_,
                   map<string, double>& AnsysParameters_,
                   string Fname_   = "",
                   string DataDir_ = "" );

    //! Virtual Destructor
    ~IOAnsysFormat() {};

    //! Load FE Model from file
    /*! Return type:        Path of file loaded
        Optional Arguments: Filename, DataDirectory */
    virtual string LoadModel ( string Fname_ ="", string DataDir_ ="" );

    //! Save FE Model to file.
    /*! Return type:        Path of file loaded
        Optional Arguments: Filename, DataDirectory */
    virtual string SaveModel( string Fname_ ="", string DataDir_ ="" );
    virtual void print_bc(StructBoundCon*, ofstream&, unsigned );
    virtual void print_bc(LcmBoundCon*, ofstream&, unsigned );

    //! Save FE-Results to file.
    /*! Return type:        Path of file loaded
        Optional Arguments: Filename, DataDirectory */
    virtual string SaveResults( string Fname_ ="", string DataDir_ ="" );

    void Quad2Tri();
    void Brick2Tet();

  private:

    //! Additional Ansys Specific Data format to be handled by reference
    map<string, double>& AnsysParameters;

  };


#include "IOAnsysFormat.inl"
#include "IOAnsysLoadModel.inl"
#include "IOAnsysSaveModel.inl"
#include "IOAnsysSaveResults.inl"


} // of namespace
#endif

