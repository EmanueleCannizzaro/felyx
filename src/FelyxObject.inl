// -*- c++ -*-
//-----------------------------------------------------------------------------
// FelyxObject.inl
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



////
//  Default Constructor
////
template <class AnalysisType>
FelyxObject<AnalysisType>::FelyxObject( unsigned noise_ )
    : noise( noise_ ), OutStream( cout ), my_timer(), DofCount( 0 ), Profile( 0 ), InterfacePtr(NULL) {
  if ( noise > 0 )
    OutStream << "Start time logging" << endl;
}

////
// File Constructor which read ANSYS file
////
template <class AnalysisType>
FelyxObject<AnalysisType>::FelyxObject( string Fname_, string DataDir_, unsigned noise_ )
    : noise( noise_ ), OutStream( cout ), my_timer(), DofCount( 0 ), Profile( 0 ) {
  if ( noise )
    OutStream << "File Constructor of FelyxObject got called; start time logging" << endl;

  // check file-type of input file
  string::size_type s1 = Fname_.find_last_of( "." );
  string::size_type s2 = Fname_.length();
  string ext( Fname_, s1 + 1, s2 );

  if ( ext == "out" ) { // patran
    InterfacePtr = new IOPatranFormat<AnalysisType>( Nodes, NodeCoordSysList, Elements, Materials, Properties, BoundaryConditions, Layers, Laminates, Parameters );;
    LoadPatranModel( Fname_, DataDir_ );
  } else if ( ext == "caa" ) {
    InterfacePtr = new IOCAAFormat<AnalysisType>( Nodes, NodeCoordSysList, Elements, Materials, Properties, BoundaryConditions, Layers, Laminates );
    LoadCAAModel( Fname_, DataDir_ );
  } else { // ansys
    InterfacePtr = new IOAnsysFormat<AnalysisType>( Nodes, NodeCoordSysList, Elements, Materials, Properties, BoundaryConditions, Layers, Laminates, Parameters );
    LoadAnsysModel( Fname_, DataDir_ );
  }

}

template <class AnalysisType>
FelyxObject<AnalysisType>::~FelyxObject(){
  if ( InterfacePtr != NULL )
    delete InterfacePtr;
}

////
//  Loading ANSYS model; optional args: filename, datadir
////
template <class AnalysisType>
string FelyxObject<AnalysisType>::LoadAnsysModel( string Fname_, string DataDir_ ) {
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Loading Ansys Model from : ";

  string path = InterfacePtr->LoadModel( Fname_, DataDir_ );

  if ( noise > 0 )
    OutStream << path << endl;

  return path;
}

////
//  Loading Patran model; optional args: filename, datadir
////
template <class AnalysisType>
string FelyxObject<AnalysisType>::LoadPatranModel( string Fname_, string DataDir_ ) {
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Loading Patran Model from : ";

  string path = InterfacePtr->LoadModel( Fname_, DataDir_ );

  if ( noise > 0 )
    OutStream << path << endl;

  return path;
}

////
//  Loading CAA model; optional args: filename, datadir
////
template <class AnalysisType>
string FelyxObject<AnalysisType>::LoadCAAModel( string Fname_, string DataDir_ ) {
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Loading CAA Model from : ";

  string path = InterfacePtr->LoadModel( Fname_, DataDir_ );

  if ( noise > 0 )
    OutStream << path << endl;

  return path;
}


////
//  Storing ANSYS model; optional args: filename, datadir
////
template <class AnalysisType>
string FelyxObject<AnalysisType>::SaveAnsysModel( string Fname_, string DataDir_ ) {
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Save FE-Model in ANSYS format to : ";

  string path = InterfacePtr->SaveModel( Fname_, DataDir_ );

  if ( noise > 0 )
    OutStream << path << endl;

  return path;
}

////
//  Storing nodal deformations; optional args: filename, datadir
////
template <class AnalysisType>
string FelyxObject<AnalysisType>::SaveResults( string Fname_, string DataDir_ ) {
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Save nodal deformations to : ";

  string path = InterfacePtr->SaveResults( Fname_, DataDir_ );

  if ( noise > 0 )
    OutStream << path << endl;

  return path;
}

////
// Nodes reordering function
////
template <class AnalysisType>
vector<unsigned> FelyxObject<AnalysisType>::NodesReordering( string bwalgo, double W1, double W2 ) {
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Nodes reordering is done, using "
    << bwalgo << " algorithm. " << endl;
  if ( bwalgo == "mmd" ) {
    fill_algos::FillReduction( Nodes, Elements, bwalgo, noise, 0 );
    //conflict signed - unsigned!!
    vector<unsigned> dummy;
    return dummy;
  } else
    return fe_base::BandwidthReduction( Nodes, Elements, bwalgo, noise, W1, W2 );
}


////
// Function which evaluates the FEM model using a direct solver (skyline solver)
////
template <class AnalysisType>
int FelyxObject<AnalysisType>::DirectSolver() {

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Link node DOF's to GSM index" << endl;
  DofCount = LinkNodes2Gsm( Nodes, Elements );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval envelope and profile of GSM : ";

  dense1D<unsigned> Envelope( DofCount );
  Profile = EvalEnvelope( Envelope, Elements );

  if ( noise > 0 )
    OutStream << Profile
    << " -> Memory needs of GSM: "
    << ( int ) ( Profile * sizeof( float_type ) / ( 1024 * 1024 ) ) << " MB " << endl;

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Initialize GSM of size : " << DofCount << endl;
  typedef mtl::matrix< float_type , symmetric<lower>, envelope<>, row_major>::type EnvelopeMatrix;
  EnvelopeMatrix GSM( Envelope, Profile );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Evaluate ESM's and assemble them to GSM" << endl;
  AssembleGM( Elements, GSM );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Resize DofSolution vector" << endl;
  DofSolution.resize( DofCount );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval load vector and apply inhomogeneous BC's" << endl;
  ApplyLoads( Nodes, GSM, DofSolution, Envelope );
  
  //std::cout << "GSM: " << "\n";
  //print_envelope2dense( GSM );
  //std::cout << "Force: " << "\n";
  //print_vector( DofSolution );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Start solving using direct skyline solver ..." << endl;
  int solved = skyline_solve( GSM, DofSolution, Envelope );
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status : " << solved << endl;

  // Store nodal deformations to nodes
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Postprocess: Storing nodal deformations to nodes " << endl;
  StoreNodalDeformations( Nodes, DofSolution );

  return solved;
}
////
// Function which evaluates the FEM model using the pardiso sparse direct solver
////
template <class AnalysisType>
int FelyxObject<AnalysisType>::SparseSolver() {
#ifdef HAVE_PARDISO

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Link node DOF's to GSM index" << endl;
  DofCount = LinkNodes2Gsm( Nodes, Elements );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Initialize GSM of size : " << DofCount << endl;
  // Pardiso needs a matrix where the internal storage indices start with 1
  mtl::matrix< float_type , symmetric<lower>, compressed<int, mtl::internal, mtl::index_from_one>, column_major >::type GM_1( DofCount );

  // Define a force vector
  DenseVector Forces( DofCount, 0 );

  // Insert a new scope, in order to get rid of temporary GM matrix before solving!
  {

    // First a sparse matrix with internal storage indices starting by zero is defined
    // since a matrix with starting index 1 is very slow for assembly (seems to be badly implemented in MTL)
    //mtl::matrix< float_type , symmetric<lower>, compressed<int, mtl::internal, mtl::index_from_one>, column_major >::type GM_0( DofCount );
    mtl::matrix< float_type, symmetric<lower>, array< compressed<> >, column_major>::type GM_0( DofCount );

    if ( noise > 0 )
      OutStream << "--< " << my_timer.elapsed() << " >-- " << "Evaluate ESM's and assemble them to GSM" << endl;
    AssembleGM( Elements, GM_0 );

    if ( noise > 0 )
      OutStream << "--< " << my_timer.elapsed() << " >-- " << "Resize DofSolution vector and set its vals to zero" << endl;
    DofSolution.resize( DofCount );
    mtl::set_value( DofSolution, float_type( 0 ) );

    //std::cout << "GM without BC's: " << "\n";
    //mtl::print_all_matrix( GM_0 );
    //std::cout << "Forces: " << "\n";
    //mtl::print_vector( Forces );

    if ( noise > 0 )
      OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval load vector and apply inhomogeneous BC's" << endl;
    ApplyBoundCons( Nodes, GM_0, Forces );

    //std::cout << "GM: " << "\n";
    //mtl::print_all_matrix( GM_0 );
    //std::cout << "Forces: " << "\n";
    //mtl::print_vector( Forces );

    // Copy 0-indexed matrix to 1-indexed matrix
    GM_1.fast_copy( GM_0 );

  } // GM_0 is not needed anymore



  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Start solving using direct sparse solver ..." << endl;

  //prepare pardiso call parameters
  int pt[ 64 ];  /*Array of 64 ints*/
  int nfact = 1; /* Number of factors kept in memory */
  int mnum = 1; /* Actual factor to compute. Constraint 1<= mnum <= nfact */
  int mtype = 2; /* Matrix type: real, symmtric positive definite matrices */
  int* perm = NULL;  /*Permutaion array of n ints*/
  int phase = 13; /* start phase - end phase */
  int n = DofCount;
  int nrhs = 1; /* number of right hand sides */
  int iparm[ 64 ]; /* Parameter array */
  //initialize these with 0
  for ( int i = 0; i < 64; i++ ) {
    iparm[ i ] = 0;
    pt[ i ] = 0;
  }
  iparm[ 0 ] = 1; /* 0: default values */
  iparm[ 1 ] = 2; /* 0: minimum degree ordering
                               2: nested dissection ordering (recommended)*/
  iparm[ 2 ] = 1; //omp_get_max_threads(); /* number of processors */
  iparm[ 4 ] = 0; /* 0: Apply permutation selected in iparm[1]
                              1: Apply custom permutation given in perm */
  iparm[ 7 ] = 2;
  iparm[ 9 ] = 13; /* Perturb the pivot elements with 1E-13 (Default value)*/
  iparm[ 10 ] = 1; /* Use nonsymmetric permutation and scaling MPS */
  iparm[ 17 ] = -1; /* <0: Report number of nonzeros in Cholesky factor in iparm[17]. Turn off for higher performance */
  iparm[ 18 ] = -1; /* <0: Report MFlops in iparm[18] for Cholesky factorization. Turn off for higher performance */
  iparm[ 19 ] = 0; /* <0: Report Numbers of CG Iterations. Turn off for higher performance */

  int msglvl = ( noise > 1 ) ? 1 : 0; /* message level 1 or 0 */
  int solved = 0; /* error state */

  double ddum; /* Double dummy */


  //The function call, splitted into 3 steps
  phase = 11;  /* first, perform only reordering step ... */
  PARDISO( pt,
           &nfact,
           &mnum,
           &mtype,
           &phase,
           &n,
           ( double* ) GM_1.get_val(),
           GM_1.get_ptr(),
           GM_1.get_ind(),
           perm,
           &nrhs,
           iparm,
           &msglvl,
           &ddum,
           &ddum,
           &solved );

  if ( solved != 0 ) {
    std::ostringstream s;
    s << "PARDISO solver: Reordering failed. Error code: " << solved;
    FELYX_RUNTIME_THROW( s.str().c_str() );
  }

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- "
    << "Reordering completed ... nonzeros in factors = " << iparm[ 17 ] << "\n";


  phase = 22; /* ... second, perform factorization ... */
  PARDISO( pt,
           &nfact,
           &mnum,
           &mtype,
           &phase,
           &n,
           ( double* ) GM_1.get_val(),       //cast float_type explicitly to double
           GM_1.get_ptr(),
           GM_1.get_ind(),
           perm,
           &nrhs,
           iparm,
           &msglvl,
           &ddum,
           &ddum,
           &solved );

  if ( solved != 0 ) {
    std::ostringstream s;
    s << "PARDISO solver: Factorization failed. Error code: " << solved;
    FELYX_RUNTIME_THROW( s.str().c_str() );
  }

  if ( noise > 0 ) {
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Factoriztion completed ... ";
    OutStream << " --> number of factorization MFLOPS = " << iparm[ 18 ];
    OutStream << " --> total peak memory consumption = " << std::max( iparm[ 15 ], iparm[ 16 ] + iparm[ 17 ] ) << std::endl;
  }

  phase = 33; /* ...and third, perform backsubstition. */
  PARDISO( pt,
           &nfact,
           &mnum,
           &mtype,
           &phase,
           &n,
           ( double* ) GM_1.get_val(),        //cast float_type explicitly to double
           GM_1.get_ptr(),
           GM_1.get_ind(),
           perm,
           &nrhs,
           iparm,
           &msglvl,
           ( double* ) Forces.data(),        //cast float_type explicitly to double
           ( double* ) DofSolution.data(),        //cast float_type explicitly to double
           &solved );

  if ( solved != 0 ) {
    std::ostringstream s;
    s << "PARDISO solver: Backsubstitution failed. Error code: " << solved;
    FELYX_RUNTIME_THROW( s.str().c_str() );
  }

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status : " << solved << endl;

  // Store nodal deformations to nodes
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Postprocess: Storing nodal deformations to nodes " << endl;
  StoreNodalDeformations( Nodes, DofSolution );

  phase = -1; /* Last but not least, release internal memory. */
  PARDISO( pt,
           &nfact,
           &mnum,
           &mtype,
           &phase,
           &n,
           &ddum,
           GM_1.get_ptr(),
           GM_1.get_ind(),
           perm,
           &nrhs,
           iparm,
           &msglvl,
           &ddum,
           &ddum,
           &solved );

  return solved;

#else

  FELYX_RUNTIME_THROW( "PARDISO solver not available. Either Intel MKL library or PARDISO libraries must be provided." );
  return -1;

#endif

}

////
// Function which evaluates the FEM model using iterative solvers from ITL
////
template <class AnalysisType>
int FelyxObject<AnalysisType>::IterativeSolver( double tolerance, string precondName, int lfil, float_type droptol, float_type shift ) {

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Link node DOF's to GSM index" << endl;
  DofCount = LinkNodes2Gsm( Nodes, Elements );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Initialize GSM of size : " << DofCount << endl;
  typedef mtl::matrix< float_type, symmetric<lower>, array< compressed<> >, column_major>::type SymSparseMatrix;
  SymSparseMatrix GSM( DofCount );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Evaluate ESM's and assemble them to GSM" << endl;
  AssembleGM( Elements, GSM );

  // Define a force vector
  DenseVector Forces( DofCount, 0 );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Resize DofSolution vector and set its vals to zero" << endl;
  // Resize DofSolution vector and set all vals to 0
  DofSolution.resize( DofCount );
  mtl::set_value( DofSolution, float_type( 0 ) );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval load vector and apply inhomogeneous BC's" << endl;
  ApplyBoundCons( Nodes, GSM, Forces );

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Start solving using iterative ITL solver ..." << endl;

  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Nonzeros: " << GSM.nnz() << endl;

  int max_iter = DofCount;
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Preconditioning ..." << endl;
  int solved = 0;
  if ( precondName == "ic0" ) {
    itl::ldl_cholesky<SymSparseMatrix> precond( GSM );
    if ( noise > 1 ) {
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Starting noisy ccg solver" << endl;
      itl::noisy_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    } else {
      itl::basic_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    }
  } else if ( precondName == "ict" ) {
    itl::ldl_ict<SymSparseMatrix> precond( GSM, lfil, droptol );
    if ( noise > 1 ) {
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Starting noisy ccg solver" << endl;
      itl::noisy_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    } else {
      itl::basic_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    }
  } else if ( precondName == "rict" ) {
    itl::ldl_robust_ict<SymSparseMatrix> precond( GSM, lfil, droptol );
    if ( noise > 1 ) {
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Starting noisy ccg solver" << endl;
      itl::noisy_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    } else {
      itl::basic_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    }
  } else if ( precondName == "sict" ) {
    itl::ldl_shift_ict<SymSparseMatrix> precond( GSM, lfil, droptol, shift );
    if ( noise > 1 ) {
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Starting noisy ccg solver" << endl;
      itl::noisy_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    } else {
      itl::basic_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " shift: " << shift << endl;
    }
  } else if ( precondName == "sainv" ) {
    itl::sainv<SymSparseMatrix> precond( GSM, droptol );

    if ( noise > 1 ) {
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Starting noisy ccg solver" << endl;
      itl::noisy_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    } else {
      itl::basic_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    }
  } else if ( precondName == "none" ) {
    itl::identity_preconditioner precond;
    if ( noise > 1 ) {
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Starting noisy ccg solver" << endl;
      itl::noisy_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    } else {
      itl::basic_iteration<double> iter( Forces, max_iter, tolerance );
      solved = itl::cg( GSM, DofSolution, Forces, precond(), iter );
      if ( noise > 0 )
        OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status: " << solved << " Iterations: " << iter.iterations() << endl;
    }
  } else {
    if ( noise > 0 )
      OutStream << "--< " << my_timer.elapsed() << " >-- " << "No valid preconditioner selected!" << endl;
  }

  // Store nodal deformations to nodes
  if ( noise > 0 )
    OutStream << "--< " << my_timer.elapsed() << " >-- " << "Postprocess: Storing nodal deformations to nodes " << endl;
  StoreNodalDeformations( Nodes, DofSolution );

  return solved;
}

template <class AnalysisType>
void FelyxObject<AnalysisType>::info() {
  std::cout << " The present finite element model contains: " << std::endl;
  std::cout << "\t-> " << Nodes.size() << " nodes " << std::endl;
  std::cout << "\t-> " << NodeCoordSysList.size() << " coordinate systems " << std::endl;
  std::cout << "\t-> " << BoundaryConditions.size() << " boundary conditions " << std::endl;
  std::cout << "\t-> " << Elements.size() << " elements " << std::endl;
  std::cout << "\t-> " << Materials.size() << " materials " << std::endl;
  std::cout << "\t-> " << Properties.size() << " properties " << std::endl;
  std::cout << "\t-> " << Layers.size() << " layers " << std::endl;
  std::cout << "\t-> " << Laminates.size() << " laminates " << std::endl;
  std::cout << "\t-> " << Parameters.size() << " parameters " << std::endl;
}


