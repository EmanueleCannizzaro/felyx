/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
 //-----------------------------------------------------------------------------
// LcmNode.h
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
   

#ifndef LcmNode_h
#define LcmNode_h LcmNode_h

#include <iostream>
#include <vector>               // Include Headers from STL <vector>

#include "CoordSys.h"		
#include "LcmBoundCon.h"		
#include "LcmDofSet.h"

// Include MTL headers
#include "mtl/mtl_felyx_utils.h"

using namespace std;
using namespace mtl;

namespace fe_base{		// Put classes into namespace fe_base
 

  //data-struct for a single node
  class LcmNode{
    
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
    LcmBoundCon *OrigBoundConPtr;
    LcmBoundCon *BoundConPtr;

    const static int NodeDofVecCount = 6;
    //! Bitset<6> containing the information which nodal degrees of 
    //! freedom are activ. (ux, uy, uz, rx, ry, rz) 
    LcmDofSet  NodeDofSet;          
    //! Integer giving the position of the first active DOF in the
    //! Global Stiffness Matrix (GSM)
    unsigned Idx2GM;
    unsigned number;

  public:
    //! Vector to store nodal deformations 
    double Fill;
    double Pressure;
    double Volume;
    double q;
    double time;
    int entrapment;
    
  public:
    // Constructors
    // -----------------------------------------
    //! Default Constructor
    LcmNode();
    //! Constructor using a vector for coordinates
    LcmNode( Dense_Vector );
    //! Constructor using doubles for coordinates
    LcmNode(double, double, double);	
    //! Copy constructor needed to avoid shallow copy of "Deformations"    
    LcmNode( const LcmNode& );
	
    // Sets & Gets
    // -----------------------------------------
    //! Set coordinates of node using MTL-vector
    void set( Dense_Vector );
    //! Set coordinates of node using doubles
    void set(double, double, double);		
    //! Set coordinates and coordsys ptr of node
    void set(double, double, double, CoordSys*);
    //! Set coordsys ptr of node
    void set(CoordSys*);	
    //! Set boundary condition ptr
    void set(LcmBoundCon*);
    void setTempBCptr(LcmBoundCon*);
    //! Get nodal coordinates as vector
    Dense_Vector Get();
    double GetCx() {return Cx;}
    double GetCy() {return Cy;}		
    double GetCz() {return Cz;}

    void setFill(double);	
    
    // Queries
    // ----------------------------------------
    //! Returns true if NodeCoordSysPtr != NULL
    bool ExistNodeCoordSys() const;		    
    //! Returns true if there are Node-BC
    bool ExistBoundCon() const;	
    bool ExistOrigBoundCon() const;				
    //! Returns bitset containing active DOF's based on DOF's of node and homogeneous BC's 
    LcmDofSet GetDofSet() const;					
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
    bool operator== (const LcmNode& ) const;
    //! Overload operator != --> based on coordinate information
    bool operator!= (const LcmNode& ) const;

    // Postprocess member functions
    // ----------------------------
    //! Store deformations of this node from the solution vector
    template< class SolVector >
    void SetPressures( const SolVector& NodalDeformations_ ){
      
      int GMindex;

      // Copy the appropriate values for this node
	GMindex = GetGMIndex(0);
	if ( GMindex >=0 )
	  Pressure = NodalDeformations_[GMindex];
	else
	  Pressure = 0;

    }
    
    //! Get nodal deformations of actual node
    /*! Gives back the Deformations stored in the node, even if this vector is empty! */
    double GetPressure() const { return Pressure; }
    Dense_Vector GetDeformations() const { Dense_Vector values(6,0.0); return values; }
  };

  // I/O
  // ---
  ostream& operator<<(ostream&, const LcmNode&);   

} // of namespace

#endif


