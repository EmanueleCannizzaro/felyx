/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// FelyxObject.h
//
// begin     : Mai 14 2002
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

#ifndef FelyxObject_h
#define FelyxObject_h FelyxObject_h

#include <sstream>

// Include MTL headers
#include <mtl/matrix.h>
#include <mtl/mtl.h>
#include <mtl/utils.h>
#include <mtl/meta_if.h>
#include "mtl/mtl_felyx_utils.h"

// Include FELyX headers
#include "ElementHeaders.h"
#include "IOAnsysFormat.h"
#include "IOPatranFormat.h"
#include "IOCAAFormat.h"

// ITL headers
#include <itl/interface/mtl.h>
#include "itl/preconditioner/ldl_shift_ict.h"
#include "itl/preconditioner/ldl_robust_ict.h"
#include "itl/preconditioner/ldl_cholesky.h"
#include "itl/preconditioner/ldl_ict.h"
#include "itl/preconditioner/sainv.h"
#include "itl/krylov/cg.h"

// More FELyX headers
#include "Bandwidth.h"
#include "BandwidthWrapper.h"
#include "Fill.h"
#include "SkylineSolver.h"
#include "PreProcessing.h"
#include "ApplyBoundCons.hpp"
#include "PostProcessing.h"
#include "Bandwidth.h"

#include <boost/timer.hpp>

extern "C" int omp_get_max_threads();
/* PARDISO prototype. */
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif
extern "C" int PARDISO
        (void *, int *, int *, int *, int *, int *,
        double *, int *, int *, int *, int *, int *,
        int *, double *, double *, int *);

using namespace std;
using namespace fe_base;
using namespace mtl;

namespace felyx{

  // ! Felyx Object class which summarizes all of a FEM evaluation in FELyX
  template<class AnalysisType>
  class FelyxObject{

  public:
    // Precision of FEM evaluation can be switched using the typedef "float_type"
    // defined in "mtl/mtl_felyx_utils.h"

    //! Typedef to define type of force and DOF solution vector
    typedef mtl::dense1D< float_type > DenseVector;

    //! Typedefs of the different entities
    typedef AnalysisType analysis_type;
    typedef typename AnalysisType::node_type node_type;
    typedef typename AnalysisType::element_type element_type;
    typedef typename AnalysisType::bc_type bc_type;
    typedef typename std::vector<PropertySet>::iterator PropertySetIter;

    // Constructors
    // ------------
    //! Default constructor with optional noise argument
    FelyxObject( unsigned = 0 );

    //! File constructor, reading in FE model
    //! Args: filenam, datadir, noiselevel (optional)
    FelyxObject( string, string, unsigned = 0);
    
    ~FelyxObject();
    
    // IO Functions
    // ------------
    //! Loading ANSYS model; optional args: filename, datadir
    string LoadAnsysModel( string = "", string = "" );
    string LoadPatranModel( string = "", string = "" );
    string LoadCAAModel( string = "", string = "" );

    //! Storing ANSYS model to file,
    string SaveAnsysModel( string = "", string = "" );

    //! Storing nodal deformations to file; optional args: filename, datadir
    string SaveResults( string = "", string = "" );

    // FEM Evaluation functions
    // ------------------------
    //! Nodes reordering function
    /*! Args: BandwidthAlgo: sloan, reversed_sloan, reversed_cuthill_mckee, cuthill_mckee
              W1, W2:             weights for sloan algorithms
    */
    std::vector<unsigned> NodesReordering( string = "sloan", double = 1, double = 2 );

    //! Function which evaluates the FEM model using a direct solver (skyline solver)
    int DirectSolver();

    //! Function which evaluates the FEM model using the pardiso direct sparse solver
    int SparseSolver();

    //! Function which evaluates the FEM model using iterative solvers from ITL
    /*! Arg: tolerance of iterative solver */
    int IterativeSolver( double = 1e-4, string = "none", int = 10, float_type = 1e-6, float_type = 0 );

    // Members to handle information printed out
    // ----------------------------------------------
    unsigned noise;         // specify amount of information printed out
    ostream& OutStream;     // specify a stream object to direct output

    //Member function to print general information on the FEModel
    void info();

    //! Object to evaluate runtimes of different parts
    boost::timer my_timer;

    // Functions to access the data members
    //-------------------------------------------
    std::vector<node_type>&    GetNodes()              { return Nodes; }
    std::vector<CoordSys>&     GetNodeCoordSys()       { return NodeCoordSysList; }
    std::vector<bc_type>&      GetBoundaryConditions() { return BoundaryConditions; }
    PtrVector<element_type*>&  GetElements()           { return Elements; }
    PtrVector<Material*>&      GetMaterials()          { return Materials; }
    std::vector<PropertySet>&  GetProperties()         { return Properties; }
    std::map<std::string,double>&   GetParameters()    { return Parameters; }
    std::vector<Layer>&        GetLayers()             { return Layers; }
    std::vector<Laminate>&     GetLaminates()          { return Laminates; }
    DenseVector&               GetDofSolution()        { return DofSolution; }
    
    unsigned GetDofCount() { return DofCount; }

    //Functions to set the data members
    //---------------------------------
    void SetLayerThickness( unsigned nl_, double value_){
      Layers[nl_].SetThickness(value_);
    }
    

  protected:
    // Data members
    // ------------
    std::vector<node_type> Nodes;
    std::vector<CoordSys> NodeCoordSysList;
    std::vector<bc_type>  BoundaryConditions;
    PtrVector<element_type*> Elements;
    PtrVector<Material*> Materials;
    std::vector<PropertySet> Properties;
    std::vector<Layer>  Layers;
    std::vector<Laminate> Laminates;
    std::map<std::string,double> Parameters;

    //! # of active DOF's
    unsigned DofCount;
    //! Profile of GSM
    unsigned Profile;
    //! MTL vector to store nodal solution
    DenseVector DofSolution;

  protected:
    //! Ansys interface class
    IOFelyx<AnalysisType>* InterfacePtr;

  };

#include "FelyxObject.inl"

}

#endif
