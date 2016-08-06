//-----------------------------------------------------------------------------
// IOCAAFormat.h
//
// begin     : Oct 2004
// copyright : (c) 2004 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
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

#ifndef IOCAAFormat_h
#define IOCAAFormat_h IOCAAFormat_h

#include <string>
#include <sstream>
#include "IOFelyx.h"

#include "IOSpiritUtils.h"
#include "boost/spirit/core.hpp"
#include "boost/spirit/utility/loops.hpp"

#define RULE(name, definition) typeof(definition) name = definition
#define bs boost::spirit

namespace fe_base
{

   //! Class to read FE-Data from CAA-Format
   template<class analysis_type> class IOCAAFormat : public IOFelyx<analysis_type>
   {

   public:

      //! Constructor, initializing references to all FE-Data objects.
      /*! Additionally, filename and datadirectory can optionally be specified */
      IOCAAFormat( vector<typename analysis_type::node_type>& Nodes_,
                     vector<CoordSys>& NodeCoordSysList_,
                     PtrVector<typename analysis_type::element_type*>&   Elements_,
                     PtrVector<Material*>&      _Materials,
                     vector<PropertySet>&               _Properties,
                     vector<typename analysis_type::bc_type>&   _BoundaryConditions,
                     vector<Layer>&               Layers_,
                     vector<Laminate>&            Lamiantes_ ,
                     string Fname_   = "",
                     string DataDir_ = "" );

      //! Virtual Destructor
      ~IOCAAFormat() {};

      //! Load FE Model from file
      /*! Return type:      Path of file loaded
         Optional Arguments: Filename, DataDirectory */
      virtual string LoadModel( string Fname_ , string DataDir_ );

      virtual std::string ParserFunction(std::string path_);

      //! Save FE-Results to file.
      /*! Return type:      Path of file loaded
          Optional Arguments: Filename, DataDirectory */
      virtual string SaveResults( string Fname_ ="", string DataDir_ ="" );


   };

#include "IOCAAFormat.inl"
#include "IOCAASpiritFunctors.inl"
#include "IOCAALoadModel.inl"
#include "IOCAASaveResults.inl"

} // of namespace
#endif

