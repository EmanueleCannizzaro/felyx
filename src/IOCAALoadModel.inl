/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// IOAnsysLoadModel.cc
//
// begin     : Oct  2004
// copyright : (c) 2004 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
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
//// Implementation of the Member function LoadModel of class IOCAAFormat
////

template<class analysis_type>
std::string IOCAAFormat<analysis_type>::ParserFunction( string path_ ) {

    char * buffer;    // the buffer, where the text file is stored for fast parsing
    unsigned long filelength = fe_base::ReadFileToMemory(path_.c_str(), buffer);

   //// ------------------------------------------------------------------------
   //// Temporary data constructs needed through this procedural member function
   //// ------------------------------------------------------------------------

   bool isLayeredShell = false;

   //objects used and filled by the SetMaterial functor
   std::map<unsigned,fe_base::Material*> MaterialNb2Ptr;
   std::vector<double> materialItems;

   //objects used and filled by the SetLayer functor
   std::map<unsigned,fe_base::Layer*> LayerNb2Ptr;
   unsigned ID =0, MaterialPropertiesID = 0;
   double thickness;

   //objects used and filled by the SetLaminate functor
   std::map<unsigned,fe_base::Laminate*> LaminateNb2Ptr;
   std::vector<unsigned> laminateItems;

   //objects used and filled by the SetNode functor
   std::map<unsigned,typename std::vector<typename analysis_type::node_type>::iterator > NodeNb2Iter;
   std::vector<double> nodeItems;


   //objects used and filled by the SetShellElement functor
   std::vector<double> orientationsVec;
   std::vector<unsigned> elementItems;
   std::string elementtype;
   unsigned NodeCount;

   //objects used and filled by the SetRestraint functor
   std::vector<std::vector<int> > restraintVec;
   std::vector<int> restraintItems;
   std::map<unsigned, CoordSys*> NodeNb2CoordSysPtr;
   std::map<unsigned,typename analysis_type::bc_type*> NodeNb2BCPtr;

   //objects used and filled by the SetShellProperty functor
   unsigned MaterialNb = 0, PropertyNb = 0;
   double ShellThickness = 0.;
   std::map<unsigned, PropertySet*> PropertyNb2Ptr;

   //objects used and filled by the CAALoadLine_p
   std::vector<double> loadItems;
   std::vector<std::vector<double> > loadVec;

   //objects used and filled by the SetCoordSys functor
   std::vector<double> axisItems;
   std::vector<std::vector<double> > axisVec;
   std::map<unsigned, CoordSys*> CoordSysNb2Ptr;



   //===============
   //The LineParsers
   //===============

   materialItems.clear();
   RULE( CAAMaterialLine_p, bs::uint_p[bs::push_back_a(materialItems)]>>
                              *(bs::blank_p >> bs::real_p[bs::push_back_a(materialItems)]) >>
                              bs::eol_p);

   RULE( CAALayerMaterialLine_p, bs::uint_p[bs::assign_a(ID)] >>
                                 bs::blank_p >>
                                 bs::uint_p[bs::assign_a(MaterialPropertiesID)] >>
                                 bs::blank_p >>
                                 bs::real_p[bs::assign_a(thickness)] >>
                                 bs::eol_p );

   laminateItems.clear();
   RULE( CAALaminateLine_p, bs::uint_p[bs::push_back_a(laminateItems)] >>
                              *(bs::blank_p >> bs::uint_p[bs::push_back_a(laminateItems)] ) >>
                              bs::eol_p );

   nodeItems.clear();
   RULE( CAANodeLine_p, bs::uint_p[bs::push_back_a(nodeItems)] >>
                        *(bs::blank_p >> bs::real_p[bs::push_back_a(nodeItems)]) >>
                        bs::eol_p );

   elementItems.clear();
   orientationsVec.clear();
   RULE( CAAShellElementLine_p,( bs::str_p("QD") | bs::str_p("TR") )[bs::assign_a(elementtype)] >>
                                 bs::uint_p[bs::assign_a(NodeCount)] >>
                                 bs::blank_p >>
                                 bs::uint_p[bs::push_back_a(elementItems)] >>
                                 bs::repeat_p(boost::ref(NodeCount))[(bs::blank_p >> bs::uint_p[bs::push_back_a(elementItems)] )] >>
                                 *(bs::blank_p >> bs::real_p[bs::push_back_a(orientationsVec)] ) >>
                                 bs::eol_p );

   restraintItems.clear();
   RULE( CAARestraintLine_p, bs::uint_p[bs::push_back_a(restraintItems)] >>
                              *(bs::blank_p >> bs::int_p[bs::push_back_a(restraintItems)]) >>
                              bs::eol_p);

   RULE( CAAPropertySetLine_p, bs::uint_p[bs::assign_a(PropertyNb)] >>
                                 bs::blank_p >>
                                 bs::uint_p[bs::assign_a(MaterialNb)] >>
                                 bs::blank_p >>
                                 bs::real_p[bs::assign_a(ShellThickness)] >>
                                 bs::eol_p);

  loadItems.clear();
   RULE ( CAALoadLine_p,   bs::uint_p[bs::push_back_a(loadItems)] >> //1,2,3,4 for transl, rot, force, moment
                           bs::blank_p >>
                           bs::uint_p[bs::push_back_a(loadItems)] >>
                           *( bs::blank_p >> bs::real_p[bs::push_back_a(loadItems)] ) >>
                           bs::eol_p );

   axisItems.clear();
   RULE ( CAAAxisLine_p,   bs::uint_p[bs::push_back_a(axisItems)] >>
                           *( bs::blank_p >> bs::real_p[bs::push_back_a(axisItems)] ) >>
                           bs::eol_p );

   RULE( comment_p, bs::ch_p('#') >>
                     *(bs::anychar_p - bs::eol_p) >>
                     bs::eol_p );
   //===============
   //The BlockParsers
   //================
   RULE( CAAIsLayeredShell_p, bs::str_p("ISLAYEREDSHELL") >>
                              bs::eol_p >>
                              bs::uint_p[bs::assign_a(isLayeredShell)] >>
                              bs::eol_p );

   RULE( CAAMaterialsBlock_p, bs::str_p("MATERIALS") >>
                              bs::eol_p >>
                              bs::uint_p[fe_base::reserve_vector_space_a(Materials)] >>
                              bs::eol_p >>
                              * ( CAAMaterialLine_p )[SetMaterial(Materials,materialItems,MaterialNb2Ptr)] );

   RULE( CAALayerMaterialsBlock_p, bs::str_p("LAYERMATERIALS") >>
                                    bs::eol_p >>
                                    bs::uint_p[fe_base::reserve_vector_space_a(Layers)] >>
                                    bs::eol_p >>
                                    *( CAALayerMaterialLine_p[SetLayer(Layers,ID,MaterialPropertiesID,thickness,MaterialNb2Ptr,LayerNb2Ptr)] ) );

   RULE( CAALaminatesBlock_p, bs::str_p("LAMINATES") >>
                              bs::eol_p >>
                              bs::uint_p[fe_base::reserve_vector_space_a(Laminates)] >>
                              bs::eol_p >>
                              *( CAALaminateLine_p[SetLaminate(Laminates,laminateItems,LayerNb2Ptr,LaminateNb2Ptr)] ) );

   RULE( CAANodeBlock_p, bs::str_p("NODES") >>
                        bs::eol_p >>
                        bs::uint_p[fe_base::reserve_vector_space_a(Nodes)] >>
                        bs::eol_p >>
                        *(CAANodeLine_p[fe_base::SetNode<typename analysis_type::node_type, typename analysis_type::bc_type>(Nodes,
                                                                              nodeItems,
                                                                              NodeNb2CoordSysPtr,
                                                                              NodeNb2BCPtr,
                                                                              NodeNb2Iter)]) );

   RULE( CAAElementBlock_p, bs::str_p("ELEMENTS") >>
                              bs::eol_p >>
                              bs::uint_p[fe_base::reserve_vector_space_a(Elements)] >>
                              bs::eol_p >>
                              *( CAAShellElementLine_p[fe_base::SetLayeredShellElement<typename analysis_type::element_type,typename analysis_type::node_type>(Elements, isLayeredShell, elementtype, NodeCount, elementItems, orientationsVec, Properties, MaterialNb2Ptr, NodeNb2Iter, LaminateNb2Ptr)]) );

   RULE( CAARestraintBlock_p, bs::str_p("RESTRAINTS") >>
                              bs::eol_p >>
                              *( CAARestraintLine_p[fe_base::push_back_vector_functor<int>(restraintItems, restraintVec)]) );

   RULE( CAAPropertySetBlock_p, bs::str_p("PROPERTYSET") >>
                                 bs::eol_p >>
                                 bs::uint_p[fe_base::reserve_vector_space_a(Properties)] >>
                                 bs::eol_p >>
                                 *( CAAPropertySetLine_p[fe_base::SetShellProperty(Properties,
                                                                                    PropertyNb,
                                                                                    MaterialNb,
                                                                                    ShellThickness,
                                                                                    MaterialNb2Ptr,
                                                                                    PropertyNb2Ptr)] ) );

   RULE( CAALoadBlock_p, bs::str_p("LOADS") >>
                         bs::eol_p >>
                         *( CAALoadLine_p[fe_base::push_back_vector_functor<double>(loadItems, loadVec)] ) );


   RULE( CAAAxisBlock_p, bs::str_p("COORDINATESYSTEMS") >>
                           bs::eol_p >>
                           bs::uint_p[fe_base::reserve_vector_space_a(NodeCoordSysList)] >>
                           bs::eol_p >>
                           *( CAAAxisLine_p[fe_base::SetCoordSys(NodeCoordSysList,
                                                                 axisItems,
                                                                 CoordSysNb2Ptr)] ) );
   //===================
   //The DocumentParser
   //===================
   RULE( CAAFile_p,
            (    CAAMaterialsBlock_p
               |  CAAIsLayeredShell_p
               |  CAALayerMaterialsBlock_p
               |  CAALaminatesBlock_p
               |  CAANodeBlock_p
               |  CAAElementBlock_p
               |  CAARestraintBlock_p[fe_base::SetAllBoundaryConditions<typename analysis_type::bc_type>(BoundaryConditions,
                                                                                                         restraintVec,
                                                                                                         loadVec,
                                                                                                         CoordSysNb2Ptr,
                                                                                                         NodeNb2CoordSysPtr,
                                                                                                         NodeNb2BCPtr)]
               |  CAAPropertySetBlock_p
               |  CAALoadBlock_p
               |  CAAAxisBlock_p
               |  comment_p
               |  bs::space_p
            )
         );

   //Parse the document
   unsigned long parselength = parse(buffer, *CAAFile_p).length;

//   std::cout << "The parsed length is " << parselength << " and the file length " << filelength << std::endl;

   if (parselength < filelength)
   {
      std::string error("->WARNING: Not the entire input file got parsed.");

      FELYX_RUNTIME_THROW( error.c_str() );
    }

   return path_;

}


template<class analysis_type>
string IOCAAFormat<analysis_type>::LoadModel ( string Fname_, string DataDir_ ) {

   //// -------------
   //// File handling
   //// -------------

   // Modify File name settings, if necessary
   SetFileNames( Fname_, DataDir_ );

   //// Read the file into memory
   ////--------------------------
    std::string path = GetLoadPath();

    return ParserFunction(path);
}

