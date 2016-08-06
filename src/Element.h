/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// Element.h
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


#ifndef Element_h
#define Element_h Element_h

#include <iostream>
#include <vector>

#include "StructNode.h"
#include "LcmNode.h"
#include "Material.h"
#include "Laminate.h"
#include "Properties.h"
#include "DofSet.h"

//mtl-includes
#include "mtl/mtl_felyx_utils.h"

#include <boost/timer.hpp>


using namespace std;
using namespace mtl;

namespace fe_base{              // Put classes into namespace fe_base

  //! Abstract base class for all elements
  /*!
   * No instances of "Element" can be created
   */
  template<class node_type, class dof_type>
  class Element{
  public:
    // Data members
    // ------------
    //! Pointer to the first node of the elements nodelist
    //node_type** NodePtr;
    std::vector< typename std::vector<node_type>::iterator > NodeVec;
    Dense_Matrix Stresses;

  private:
    //! Flag holding the element type reference number
    unsigned RefNumber;

    //! Flag to check status of element (set to true in constructor)
    bool Active;


  public:
    // Constructors
    // ------------
    Element()                   : NodeVec(0),  RefNumber(1),  Active(1)  {}
    Element(unsigned n)         : RefNumber(1),  Active(1)  { create(n); }
    Element( const Element&);
    virtual ~Element()                  { destroy(); }
    //! Virtual clone function, needed in PtrVector
    virtual Element* Clone() const =0;

    // Create and destroy elements
    // ---------------------------
    void create(unsigned);
    void destroy();

    // Set members
    // -------------

    //! Set Node ptr of ith node, returning true if range check for i succeeded
    bool SetNodeIter( unsigned i, typename std::vector<node_type>::iterator Iter);

    //! Set element to active / not active
    void SetStatus(bool b) { Active = b; }
    //! Set element reference number
    void SetRefNumber( unsigned r ) { RefNumber = r; }


    // Get members
    // -----------
    //!Returns the n-th node.
    typename std::vector<node_type>::iterator    GetNodeIter(int i)    const {return NodeVec[i];}
    //! Check if element is active.
    bool              GetStatus()       const {return Active; }
    //! Check reference number
    unsigned              GetRefNumber()    const {return RefNumber;}

    virtual unsigned      GetNodeCount()  const = 0;
    virtual unsigned      GetId()         const = 0;
    virtual string        GetName()       const = 0;
    virtual dof_type      GetDofSet()     const = 0;
    virtual unsigned      GetEMSize()     const = 0;

    // bool GetNodeCoordSys(Dense_Vector&) const {};
    virtual Material*         GetMaterialPtr() const = 0;
    virtual PropertySet*      GetPropertiesPtr() const = 0;
    virtual CoordSys*         GetEleCoordSysPtr() const = 0;
    virtual Laminate*         GetLaminatePtr() const = 0;


    virtual void SetMaterialPtr( Material* a){};
    virtual void SetLaminatePtr( Laminate* a){};
    virtual void SetPropertiesPtr( PropertySet* a ){};
    virtual void SetEleCoordSysPtr( CoordSys* a){};
    virtual void AppendMaterialOrientation( MaterialOrientation a ) {};


    double GetCx(){
      double Cx=0;
      for (unsigned i=0; i<GetNodeCount(); ++i)
        Cx += (*GetNodeIter(i)).Cx;
      Cx /= GetNodeCount(); return Cx; }

    double GetCy(){
      double Cy=0;
      for (unsigned i=0; i<GetNodeCount(); ++i)
        Cy += (*GetNodeIter(i)).Cy;
      Cy /= GetNodeCount(); return Cy; }

    double GetCz(){
      double Cz=0;
      for (unsigned i=0; i<GetNodeCount(); ++i)
        Cz += (*GetNodeIter(i)).Cz;
      Cz /= GetNodeCount(); return Cz; }

    Dense_Vector GetCOG(){
      Dense_Vector COG(3,0.0);
      unsigned my_node_count = GetNodeCount();

      for(unsigned it = 0; it < my_node_count; ++it)
        {
          mtl::add((*GetNodeIter(it)).Get(), COG, COG);
        }

      mtl::scale(COG, 1.0/double(my_node_count));

      return COG;
    }

    // Some other functions
    //---------------------

    //!Gives back the sum of the first three DOF's
    /*!This gives a feedback on the dimension of the element*/
    //unsigned GetDimension();

 //    // I/O
//     // ---
//     // Provide << operator for all derived classes
//     // --> operator must be public
//     // Make use of the "Virtual Friend Function Idiom"
//     template<class bla>
//     friend ostream& operator<<(ostream&, const bla&);

    // EVAL functions
    // --------------
    //! Eval element matrix (EM)
    virtual Dense_Matrix CalcEM() = 0;

    //! Eval Element Mass Matrix (EMM)
    //virtual Dense_Matrix CalcEmm() = 0;

    //! Eval volume of an element
    //virtual double EvalVolume() const;

    //! Set NodeDofSets based on the attached elements
    //void SetNodeDofSet() const;

    //! Transform element matrix (EM), if any nodes of this element have rotated coord systems
    void TransformEM( Dense_Matrix EM );
    //! Eval influence on the envelope of the GSM that actual element will have, store bandwidth of
    //! involved rows, if their sizes get increased.
    template<class Array>
    void EvalEnvelope( Array& ) const;

    //! Assemble element to Global Stiffness Matrix (GSM)
    template<class mtlMatrix>
    void AssembleElement2GM( mtlMatrix& );

//     //! Assemble element to Global Mass Matrix (GMM)
//     template<class mtlMatrix>
//     void AssembleElement2Gmm( mtlMatrix& );

  };//of class Element


  // Include Element.inl, where template functions are implemented
  #include "Element.inl"

} // of namespace

#endif //Element_h

