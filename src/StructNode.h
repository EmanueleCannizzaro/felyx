/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
 //-----------------------------------------------------------------------------
// StructNode.h
//
// begin     : 12 1 2003
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


#ifndef StructNode_h
#define StructNode_h StructNode_h

#include <iostream>
#include <vector>               // Include Headers from STL <vector>

#include "CoordSys.h"
#include "StructBoundCon.h"
#include "StructDofSet.h"
#include "FelyxException.h"

// Include MTL headers
#include "mtl/mtl_felyx_utils.h"

using namespace std;
using namespace mtl;

namespace fe_base{              // Put classes into namespace fe_base


  //data-struct for a single node
  class StructNode{

  public:
    // Data members
    // --------------------------------------
    //!x-coordinate
    double Cx;
    //!y-coordinate
    double Cy;
    //!z-coordinate
    double Cz;

    //! Pointer to nodal coordinate system
    //! (Pointer is set to NULL, if nodal coord sys = global coord sys)
    CoordSys *NodeCoordSysPtr;
    //! Pointer to the boundary conditions of the node
    //! (Pointer is set to NULL if no BD are present)
    StructBoundCon *BoundConPtr;

    const static int NodeDofVecCount = 6;
    //! Bitset<6> containing the information which nodal degrees of
    //! freedom are activ. (ux, uy, uz, rx, ry, rz)
    StructDofSet  NodeDofSet;
    //! Integer giving the position of the first active DOF in the
    //! Global Stiffness Matrix (GSM)
    unsigned Idx2GM;

  private:
    //! Vector to store nodal deformations
    Dense_Vector Deformations;
    Dense_Vector Stresses;

  public:
    unsigned counter;

    // Constructors
    // -----------------------------------------
    //! Default Constructor
    StructNode();
    //! Constructor using a vector for coordinates
    StructNode( Dense_Vector );
    //! Constructor using doubles for coordinates
    StructNode(double, double, double);
    //! Copy constructor needed to avoid shallow copy of "Deformations"
    StructNode( const StructNode& );

    // Sets & Gets
    // -----------------------------------------
    //! Set coordinates of node using MTL-vector
    void set( Dense_Vector );
    //! Set coordinates of node using doubles
    void set(double, double, double);
    //! Set single coordinates of node using doubles
    void SetX(double value) { Cx = value; };
    void SetY(double value) { Cy = value; };
    void SetZ(double value) { Cz = value; };
    //! Set coordinates and coordsys ptr of node
    void set(double, double, double, CoordSys*);
    //! Set coordsys ptr of node
    void set(CoordSys*);
    //! Set boundary condition ptr
    void set(StructBoundCon*);
    //! Get nodal coordinates as vector
    Dense_Vector Get();
    double GetCx() const {return Cx;}
    double GetCy() const {return Cy;}
    double GetCz() const {return Cz;}

    // Queries
    // ----------------------------------------
    //! Returns true if NodeCoordSysPtr != NULL
    bool ExistNodeCoordSys() const;
    //! Returns true if there are Node-BC
    bool ExistBoundCon() const;
    //! Returns bitset containing active DOF's based on DOF's of node and homogeneous BC's
    StructDofSet GetDofSet() const;
    //! Given a certain DOF (0-5 -> x,y,z,rx,ry,rz) function returns global index in the GSM,
    //! if this DOF is active, else it returns -1
    int GetGMIndex(unsigned) const;

    //! Return coordinates of node in MTL-vector
    Dense_Vector GetCoords() const ;

    //! Return coordinate system of node, if any exists
    bool GetNodeCoordSysEuler312(Dense_Vector& ) const;
    bool GetNodeCoordSysEuler313(Dense_Vector& ) const;

    // Operator overloading
    // --------------------
    //! Overload operator == --> true if coordinates are identical
    bool operator== (const StructNode& ) const;
    //! Overload operator != --> based on coordinate information
    bool operator!= (const StructNode& ) const;

    // Postprocess member functions
    // ----------------------------
    //! Store deformations of this node from the solution vector
    template< class SolVector >
    void SetDeformations( const SolVector& NodalDeformations_ ){

      Deformations.resize(6);
      int GMindex;

      // Copy the appropriate values for this node
      for (unsigned i=0; i < 6; ++i){
        GMindex = GetGMIndex(i);
        if ( GMindex >=0 )
          Deformations[i] = NodalDeformations_[GMindex];
        else
          Deformations[i] = 0;
      }
    }

    template< class SolVector >
    void SetStresses( const SolVector& NodalStresses_ ){

      unsigned vsize = NodalStresses_.size();
      Stresses.resize(vsize+1);

      // Copy the appropriate values for this node
      for (unsigned i=0; i < vsize; ++i){
        Stresses[i] = (Stresses[vsize]*Stresses[i] + NodalStresses_[i]) / (Stresses[vsize]+1);
      }
      Stresses[vsize] = Stresses[vsize]+1;
    }

    float_type GetVMStress() const;

    //! Get nodal deformations of actual node
    /*! Gives back the Deformations stored in the node, even if this vector is empty! */
    Dense_Vector GetDeformations() const { return Deformations; }
    Dense_Vector GetTranslations() const;
    Dense_Vector GetStresses() const { return Stresses; }

  };
  // I/O
  // ---
  ostream& operator<<(ostream&, const StructNode&);

} // of namespace

#endif


