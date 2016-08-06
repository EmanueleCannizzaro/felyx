//-----------------------------------------------------------------------------
// IOAnsysFormat.cc
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
IOAnsysFormat<analysis_type>::IOAnsysFormat( vector<typename analysis_type::node_type>& Nodes_,
                                             vector<CoordSys>& NodeCoordSysList_,
                                             PtrVector<typename analysis_type::element_type*>&        Elements_,
                                             PtrVector<Material*>&	        Materials_,
                                             vector<PropertySet>&	        Properties_,
                                             vector<typename analysis_type::bc_type>& 		BoundaryConditions_,
                                             vector<Layer>&                    Layers_,
                                             vector<Laminate>& Laminates_,
                                             map<string, double>& AnsysParameters_,
                                             string Fname_, 
                                             string DataDir_ ) 
  : IOFelyx<analysis_type>( Nodes_, NodeCoordSysList_, Elements_, Materials_, Properties_, BoundaryConditions_, Layers_, Laminates_ ),
    AnsysParameters ( AnsysParameters_ )
{
  SetFileNames( Fname_, DataDir_, ".ansys", ".ansysnew", ".felyxres" );
}

////
//// Member functions
////
