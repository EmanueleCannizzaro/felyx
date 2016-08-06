//-----------------------------------------------------------------------------
// IOLcmBc.h
//
// begin     : Nov 24 2003
// copyright : (c) 2003 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
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
   
#ifndef IOLcmBc_h
#define IOLcmBc_h IOLcmBc_h

#include <string>
#include <iomanip>
#include "IOFelyx.h"
#include "LcmAnalysisType.h"

using namespace std;

namespace fe_base{
  
  //! Class to store BC-Data
  class IOLcmBc : public IOFelyx<LcmAnalysisType>
  {
    
  public:

    //! Constructor, initializing references to all FE-Data objects.
    /*! Additionally, filename and datadirectory can optionally be specified */
    IOLcmBc( vector<LcmNode>&		_Nodes,
		   vector<CoordSys>& 	        NodeCoordSysList_,
		   PtrVector<LcmElement*>&   Elements_,
		   PtrVector<Material*>&	_Materials,
		   vector<PropertySet>&		_Properties,
		   vector<LcmBoundCon>&	_BoundaryConditions,
		   vector<Layer>&               Layers_,
           vector<Laminate>& Laminates_,
		   string 			_Fname   = "",
		   string 			_DataDir = "../fedata" );

    //! Virtual Destructor
    ~IOLcmBc() {};
    
    //! Load temporary Boundary Conditions from file
    virtual string LoadModel ( string, string, vector<TdepBC>& );
   
  };


} // of namespace
#endif

