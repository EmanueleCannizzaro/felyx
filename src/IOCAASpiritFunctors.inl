//-----------------------------------------------------------------------------
// IOCAASpiritFunctors.h
//
// begin     : Nov 2004
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



/*!
   An actor to set the Material
*/
struct SetMaterial
{
   public:
      SetMaterial(fe_base::PtrVector<fe_base::Material*>& materialVec_,
            std::vector<double>& itemsVec_,
            std::map<unsigned, fe_base::Material* >& map_ ):
      materialVec(materialVec_), itemsVec(itemsVec_), map(map_)
      {
      }

   void operator() (const char *first, const char *last) const
   {

   // it is assumed that the following values are written to the itemsVec starting at position itemsVec[1]
   // E1, E2, E3, G12, G31, G23, nu12, nu13, nu23, rho, E_Ratio
   // non-relevant entries are filled with 0
      if ( itemsVec[2] == 0. ) //is isotropic
      {
         fe_base::IsotropicMaterial IsotropicMaterial(itemsVec[1], itemsVec[7], itemsVec[10]);
         materialVec.push_back(&IsotropicMaterial);
      }
      else if ( itemsVec[1] != itemsVec[2] && itemsVec[2] != 0.) //is transversely orthotropic
      {
         fe_base::TransverseIsotropicMaterial23 TransvMaterial(itemsVec[1],itemsVec[2],
                                                               itemsVec[4], itemsVec[6],
                                                               itemsVec[7], itemsVec[10]);
         materialVec.push_back( &TransvMaterial );
      }
      else if ( itemsVec[1] == itemsVec[2] && itemsVec[2] != 0.) //is a weave material
      {
         fe_base::WeaveMaterial WeaveMaterial(itemsVec[1],itemsVec[2],
                                             itemsVec[4], itemsVec[7],
                                             itemsVec[10], itemsVec[11]);
         materialVec.push_back( &WeaveMaterial);
      }
      fe_base::PtrVector<fe_base::Material*>::iterator myIter = materialVec.end();
      --myIter;
      map[(unsigned)itemsVec[0]] = *myIter;

      itemsVec.clear();
      }

   private:
      fe_base::PtrVector< fe_base::Material* >& materialVec;
      std::vector<double>& itemsVec;
      std::map<unsigned, fe_base::Material* >& map;
};


   /*!
   Set the layer
   */
struct SetLayer
{
   public:
      SetLayer(std::vector<fe_base::Layer>& LayerVec_,
               unsigned& ID_,
               unsigned& MaterialPropertiesID_,
               double& thickness_,
               std::map<unsigned, fe_base::Material* >& materialMap_ ,
               std::map<unsigned, fe_base::Layer* >& layerMap_ ):
      LayerVec(LayerVec_), ID(ID_), MaterialPropertiesID(MaterialPropertiesID_),
      thickness(thickness_), materialMap(materialMap_),
      layerMap(layerMap_)
      {
      }

   void operator() (const char *first, const char *last) const
   {
      fe_base::Layer tmpLayer(materialMap[MaterialPropertiesID], 0., thickness);
      LayerVec.push_back(tmpLayer);
      std::vector<fe_base::Layer>::iterator layerIter = LayerVec.end();
      --layerIter;
      layerMap[ID] = &(*layerIter);
   }

   private:
      std::vector<fe_base::Layer>& LayerVec;
      unsigned& ID;
      unsigned& MaterialPropertiesID;
      double& thickness;
      std::map<unsigned,fe_base::Material* >& materialMap;
      std::map<unsigned, fe_base::Layer* >& layerMap;
};

   /*!
   An actor to set the laminates
   */
struct SetLaminate
{
   public:
      SetLaminate(std::vector<fe_base::Laminate>& LaminateVec_,
                  std::vector<unsigned>& laminateItems_,
                  std::map<unsigned, fe_base::Layer* >& layerMap_ ,
                  std::map<unsigned, fe_base::Laminate* >& laminateMap_ ):
      LaminateVec(LaminateVec_), laminateItems(laminateItems_),
      layerMap(layerMap_),laminateMap(laminateMap_)
      {
      }

   void operator() (const char *first, const char *last) const
   {
      fe_base::Laminate tmpLaminate;
      for ( unsigned i = 1 ; i < laminateItems.size() ; ++i )
         tmpLaminate.push_back(layerMap[laminateItems[i]]);

      LaminateVec.push_back(tmpLaminate);

      std::vector<fe_base::Laminate>::iterator laminateIter = LaminateVec.end();
      --laminateIter;
      laminateMap[laminateItems[0]] = &(*laminateIter);

      laminateItems.clear();
   }

   private:
      std::vector<fe_base::Laminate>& LaminateVec;
      std::vector<unsigned>& laminateItems;
      std::map<unsigned, fe_base::Layer* >& layerMap;
      std::map<unsigned, fe_base::Laminate* >& laminateMap;
};

/*!
   An actor to set the coordinates of a node
   the vector contains the nodenumber and the x,y, and z coordinates
   in the order listed
*/
template <class NodeT, class bcT>
struct SetNode
{
   public:
      SetNode(std::vector<NodeT>& nodeVec_,
            std::vector<double>& itemsVec_,
            std::map<unsigned, CoordSys*>& NodeNb2CoordSysPtr_,
            std::map<unsigned, bcT*>& NodeNb2BCPtr_,
            std::map<unsigned,typename std::vector<NodeT>::iterator>& NodeNb2Iter_ ):
      nodeVec(nodeVec_),
      itemsVec(itemsVec_),
      NodeNb2CoordSysPtr(NodeNb2CoordSysPtr_),
      NodeNb2BCPtr(NodeNb2BCPtr_),
      NodeNb2Iter(NodeNb2Iter_)
      {}

   void operator() (const char *first, const char *last) const
   {
      NodeT tmpNode(itemsVec[1],itemsVec[2],itemsVec[3]);

      //look for an axis system attached to that node
      std::map<unsigned, CoordSys*>::iterator csIter;
      csIter = NodeNb2CoordSysPtr.find((unsigned)itemsVec[0]);
      if ( csIter != NodeNb2CoordSysPtr.end() ) // no specified coordsys found
      {
         tmpNode.set(csIter->second);
      }

      //look for a boundary condition and attach to the node if any
      typename std::map<unsigned, bcT*>::iterator bcIter;
      bcIter = NodeNb2BCPtr.find((unsigned)itemsVec[0]);

      if ( bcIter != NodeNb2BCPtr.end() ) // no specified bc found
      {
          tmpNode.set(bcIter->second);
      }

      nodeVec.push_back(tmpNode);
      typename std::vector<NodeT>::iterator myNode = nodeVec.end();
      --myNode;
      NodeNb2Iter[(unsigned)itemsVec[0]] = myNode;

      itemsVec.clear();
   }

   private:
      std::vector<NodeT>& nodeVec;
      std::vector<double>& itemsVec;
      std::map<unsigned, CoordSys*>& NodeNb2CoordSysPtr;
      std::map<unsigned,bcT*>& NodeNb2BCPtr;
      std::map<unsigned,typename std::vector<NodeT>::iterator>& NodeNb2Iter;
};

template <class ElementT, class NodeT>
struct SetLayeredShellElement
{
   public:
      SetLayeredShellElement(fe_base::PtrVector<ElementT*>& elementVec_,
                              bool& isLayeredShell_,
                              std::string& shelltype_,
                              unsigned& NodeCount_,
                              std::vector<unsigned>& itemsVec_,
                              std::vector<double>& orientationsVec_,
                              std::vector<fe_base::PropertySet>& Properties_,
                              std::map<unsigned, fe_base::Material*>& MaterialNb2Ptr_,
                              std::map<unsigned,typename std::vector<NodeT>::iterator>& NodeMap_,
                              std::map<unsigned,fe_base::Laminate*>& LaminateMap_):
      elementVec(elementVec_),
      isLayeredShell(isLayeredShell_),
      shelltype(shelltype_),
      NodeCount(NodeCount_),
      itemsVec(itemsVec_),
      orientationsVec(orientationsVec_),
      Properties(Properties_),
      MaterialNb2Ptr(MaterialNb2Ptr_),
      NodeMap(NodeMap_),LaminateMap(LaminateMap_)
      {}

   void operator() (const char *first, const char *last) const
   {
      std::stringstream eleType;
      eleType << shelltype << NodeCount;

      if ( isLayeredShell )
      {
        unsigned elementID = 0;
        if ( eleType.str() == "TR6" )
          elementID = 391;
        else if ( eleType.str() == "QD8" )
          elementID = 91;

        //std::cout << NodeMap.size() << "\t" << shelltype << std::endl;

        //ElementT* elePtr = StructElementFactory::Instance().CreateObject( elementID ) ;;

        ElementT *elePtr = NULL;
        fe_base::create_element(elePtr, elementID);
        elementVec.push_back(elePtr);
        typename fe_base::PtrVector<ElementT*>::iterator myIter = elementVec.end();
        --myIter;

        //Set the laminate pointer in the element
        (*myIter)->SetLaminatePtr(LaminateMap[itemsVec[0]]);

        //std::cout << " The laminate number is " << itemsVec[0] << std::endl;

        //Set the node pointers in the element
        for (unsigned i=0; i < (*myIter)->GetNodeCount()  ; ++i)
          (*myIter)->SetNodeIter(i , NodeMap[itemsVec[1+i]] );

        //add the orientations to the element
        for ( unsigned mo = 0 ; mo < orientationsVec.size() ; mo+=2 )
        {
          //std::cout << "found some orientations " << orientationsVec[mo] << "\t" << orientationsVec[mo+1] << std::endl;
          MaterialOrientation tmpOrient(orientationsVec[mo], orientationsVec[mo+1]);
          (*myIter)->AppendMaterialOrientation(tmpOrient);
        }
      }//end of is layered shell
      else //it is a nonlayered shell
      {
          //std::cout << "Setting a nonlayered shell " << std::endl;
          unsigned elementID = 0;
          if ( eleType.str() == "TR6" )
            elementID = 393;
          else if ( eleType.str() == "QD8" )
            elementID = 93;

         ElementT *elePtr = NULL;
         fe_base::create_element(elePtr, elementID);
         elementVec.push_back(elePtr);
         typename fe_base::PtrVector<ElementT*>::iterator myIter = elementVec.end();
         --myIter;

          std::vector<fe_base::PropertySet>::iterator propIter = Properties.begin();
          bool found = false;
          while ( propIter != Properties.end() && !found )
          {
              if ( (unsigned)propIter->GetId() == itemsVec[0] )
                found = true;
              ++propIter;
          }
          --propIter;
          if ( propIter != Properties.end() )
          {
            //Set the properties pointer in the element
            (*myIter)->SetPropertiesPtr(&(*propIter));
            (*myIter)->SetMaterialPtr(propIter->MaterialPtr);

            //Set the node pointers in the element
            for (unsigned i=0; i < (*myIter)->GetNodeCount()  ; ++i)
              (*myIter)->SetNodeIter(i , NodeMap[itemsVec[1+i]] );
          }
          else
          std::cout << "ERROR: Some Shell elements have no property assigned " << std::endl;


      }//end of not layeredShell

      shelltype.clear();
      orientationsVec.clear();
      itemsVec.clear();
   }

   private:
      fe_base::PtrVector<ElementT*>& elementVec;
      bool& isLayeredShell;
      std::string& shelltype;
      unsigned& NodeCount;
      std::vector<unsigned>& itemsVec;
      std::vector<double>& orientationsVec;
      std::vector<fe_base::PropertySet>& Properties;
      std::map<unsigned, fe_base::Material*>& MaterialNb2Ptr;
      std::map<unsigned,typename std::vector<NodeT>::iterator>& NodeMap;
      std::map<unsigned,fe_base::Laminate*>& LaminateMap;
};

template <class bcT>
struct SetAllBoundaryConditions
{
   public:
      SetAllBoundaryConditions(std::vector<bcT>& BCVec_ ,
                               std::vector<std::vector<int> >& restraintsVec_,
                               std::vector<std::vector<double> >& loadsVec_,
                               std::map<unsigned, CoordSys*>& CoordSysNb2Ptr_,
                               std::map<unsigned, CoordSys*>& NodeNb2CoordSysPtr_,
                               std::map<unsigned, bcT*>& NodeNb2BCPtr_):
      BCVec(BCVec_),
      restraintsVec(restraintsVec_),
      loadsVec(loadsVec_),
      CoordSysNb2Ptr(CoordSysNb2Ptr_),
      NodeNb2CoordSysPtr(NodeNb2CoordSysPtr_),
      NodeNb2BCPtr(NodeNb2BCPtr_)
      {}

      void operator() (const char *first, const char *last) const
      {

          // looping all loads checking for restraints (and therefore possible coordSys) on the same node
          std::vector<std::vector<double> >::iterator loadsIter = loadsVec.begin();
          std::vector<std::vector<int> >::iterator restraintsIter = restraintsVec.begin();
          std::map< std::vector<std::vector<double> >::iterator, std::vector<std::vector<int> >::iterator > loadIter2restraintIter;
          bool found = false;
          while (loadsIter != loadsVec.end())
          {
            restraintsIter = restraintsVec.begin();
            found = false;
            while ( restraintsIter != restraintsVec.end()  && !found )
            {
              if ( (*restraintsIter)[0] == (*loadsIter)[1] )
              {
                  loadIter2restraintIter[loadsIter] = restraintsIter;
                  found = true; //only one restraint set per load is possible
              }
              ++restraintsIter;
            }
            ++loadsIter;
          }

          BCVec.reserve(restraintsVec.size()+loadsVec.size()- loadIter2restraintIter.size());

          //looping all restraints
          restraintsIter = restraintsVec.begin();
          while (restraintsIter != restraintsVec.end())
          {
            bcT tmpBC;
            //Set the coordinates systemPtr
            if( (*restraintsIter)[1] >= 0 )
                NodeNb2CoordSysPtr[(*restraintsIter)[0]] = CoordSysNb2Ptr[(*restraintsIter)[1]];

            bool hasLoads = false;
            //look for attached loads
            loadsIter = loadsVec.begin();
            while (loadsIter != loadsVec.end())
            {
              if (loadIter2restraintIter[loadsIter] == restraintsIter)
              {
                  if ( (*loadsIter)[0] == 1 ) // a translation is set
                  {
                      for (unsigned tt = 0 ; tt < 3 ; ++tt )
                      {
                          if ( (*restraintsIter)[tt+2] ) //only set, if it is constrained
                          {
                             tmpBC.set((typename bcT::enumType)tt, (*loadsIter)[tt+2]);
                          }
                      }
                      hasLoads = true;
                      (*loadsIter)[0] = -1.; //mark this load as
                  }
                  else if ( (*loadsIter)[0] == 2 ) // a rotation is set
                  {
                    for (unsigned rr = 0 ; rr < 3 ; ++rr )
                    {
                        if ( (*restraintsIter)[rr+5] ) //only set, if it is constrained
                            tmpBC.set((typename bcT::enumType)(rr+3), (*loadsIter)[rr+2]);
                    }
                    hasLoads = true;
                    (*loadsIter)[0] = -2.; //mark this load as
                  }
                  else if ( (*loadsIter)[0] == 3 ) // a force is set
                  {
                    for (unsigned ff = 0 ; ff < 3 ; ++ff )
                          tmpBC.set((typename bcT::enumType)(ff+6), (*loadsIter)[ff+2]);
                    (*loadsIter)[0] = -3.; //mark this load as
                  }
                  else if ( (*loadsIter)[0] == 4 ) // a moment is set
                  {
                    for (unsigned mm = 0 ; mm < 3 ; ++mm )
                          tmpBC.set((typename bcT::enumType)(mm+9), (*loadsIter)[mm+2]);
                    (*loadsIter)[0] = -4.; //mark this load as
                  }
              }
              ++loadsIter;
            }
            if ( !hasLoads ) //no translation or rotation assigned
            {
                for (unsigned kk = 0 ; kk < 6 ; ++kk )
                {
                  if ( (*restraintsIter)[kk+2] ) //only set, if it is constrained
                      tmpBC.set((typename bcT::enumType)kk, 0.);
                }
            }

            typename std::vector<bcT>::iterator findIter;
            // check if the bc is allready present
            // and assign the node to bc ptr
            findIter = find(BCVec.begin(), BCVec.end(), tmpBC);
            if ( findIter == BCVec.end() )
            {
                BCVec.push_back(tmpBC);
                findIter = BCVec.end();
                --findIter;
                NodeNb2BCPtr[(*restraintsIter)[0]] = &(*findIter);
            }
            else
                NodeNb2BCPtr[(*restraintsIter)[0]] = &(*findIter);

            ++restraintsIter;
          }

          //loop over the remaining loads which have no restraint
          //only forces and moments in global coordinates are read here
          loadsIter = loadsVec.begin();
          std::vector<std::vector<double> >::iterator secondloadsIter;
          while ( loadsIter != loadsVec.end() )
          {

            if ( (*loadsIter)[0] > 0 )
            {
              bcT tmpBC;
              if ( (*loadsIter)[0] == 3 )
              {
                for (unsigned ff = 0 ; ff < 3 ; ++ff )
                      tmpBC.set((typename bcT::enumType)(ff+6), (*loadsIter)[ff+2]);
                (*loadsIter)[0] = -3.; //mark this load as used
              }

              secondloadsIter = loadsVec.begin();
              while ( secondloadsIter != loadsVec.end() )
              {
                if ( (*secondloadsIter)[0] == 4 && (*secondloadsIter)[1] == (*loadsIter)[1] )
                {
                    for (unsigned mm = 0 ; mm < 3 ; ++mm )
                          tmpBC.set((typename bcT::enumType)(mm+9), (*secondloadsIter)[mm+2]);
                    (*loadsIter)[0] = -4.; //mark this load as used
                }
                ++secondloadsIter;
              }
              typename std::vector<bcT>::iterator findIter;
              // check if the bc is allready present
              // and assign the node to bc ptr
              findIter = find(BCVec.begin(), BCVec.end(), tmpBC);
              if ( findIter == BCVec.end() )
              {
                  BCVec.push_back(tmpBC);
                  findIter = BCVec.end();
                  --findIter;
                  NodeNb2BCPtr[(unsigned)(*loadsIter)[1]] = &(*findIter);
              }
              else
                  NodeNb2BCPtr[(unsigned)(*loadsIter)[1]] = &(*findIter);
            }
            ++loadsIter;
          }
    }// end of operator()


   private:
      std::vector<bcT>& BCVec;
      std::vector<std::vector<int> >& restraintsVec;
      std::vector<std::vector<double> >& loadsVec;
      std::map<unsigned, CoordSys*>& CoordSysNb2Ptr;
      std::map<unsigned, CoordSys*>& NodeNb2CoordSysPtr;
      std::map<unsigned, bcT*>& NodeNb2BCPtr;
};


struct SetShellProperty
{
   public:

      SetShellProperty( std::vector<PropertySet>& Properties_,
                        unsigned& PropertyNb_,
                        unsigned& MaterialNb_,
                        double& ShellThickness_,
                        std::map<unsigned, fe_base::Material* >& MaterialNb2Ptr_,
                        std::map<unsigned, PropertySet*>& PropertyNb2Ptr_ ):
      Properties(Properties_),
      PropertyNb(PropertyNb_),
      MaterialNb(MaterialNb_),
      ShellThickness(ShellThickness_),
      MaterialNb2Ptr(MaterialNb2Ptr_),
      PropertyNb2Ptr(PropertyNb2Ptr_)
      {}

      void operator() (const char *first, const char *last) const
      {
         PropertySet tmpProperty(MaterialNb2Ptr[MaterialNb], PropertyNb);
         tmpProperty.Set("Thickness", ShellThickness);
         Properties.push_back(tmpProperty);
         std::vector<PropertySet>::iterator iter = Properties.end();
         --iter;
         PropertyNb2Ptr[PropertyNb] = &(*iter);
      }

   private:
      std::vector<PropertySet>& Properties;
      unsigned& PropertyNb;
      unsigned& MaterialNb;
      double& ShellThickness;
      std::map<unsigned, fe_base::Material* >& MaterialNb2Ptr;
      std::map<unsigned, PropertySet*>& PropertyNb2Ptr;
};

struct SetCoordSys
{
   public:
      SetCoordSys( std::vector<CoordSys>& CoordSysVec_,
                     std::vector<double>& axisItems_,
                     std::map<unsigned, CoordSys*>& CoordSysNb2Ptr_):
      CoordSysVec(CoordSysVec_),
      axisItems(axisItems_),
      CoordSysNb2Ptr(CoordSysNb2Ptr_)
      {}

   void operator() (const char *first, const char *last) const
   {
      CoordSys tmpAxis;
      //set the 313 euler angles
      tmpAxis.SetType(CoordSys::cartesian);
      tmpAxis.SetAngleType(CoordSys::rad);
      tmpAxis.Set( CoordSys::Euler313, axisItems[1], axisItems[2], axisItems[3] );

      CoordSysVec.push_back(tmpAxis);
      std::vector<CoordSys>::iterator iter = CoordSysVec.end();
      --iter;
      CoordSysNb2Ptr[(unsigned)axisItems[0]] = &(*iter);
      axisItems.clear();
   }

   private:
      std::vector<CoordSys>& CoordSysVec;
      std::vector<double>& axisItems;
      std::map<unsigned, CoordSys*>& CoordSysNb2Ptr;
};
