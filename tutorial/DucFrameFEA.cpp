//-----------------------------------------------------------------------------
// StructObjectTest.cpp
//
// begin     : Jan 2 2003
// copyright : (c) 2001 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// FELyX Main header
#include "DucFelyx.h"

// Include program options lib from boost
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Include filesystem lib from boost
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
namespace fs = boost::filesystem;



int main( int ac, char* av[] ) {

  try {

    std::cout << "FE analysis of Ducati steel trellis frame using FELyX\n";
   
    // PROGRAM OPTIONS
    // ---------------
    // Declare a group of options that will be allowed only on command line
    po::options_description generic( "Generic options" );
    unsigned noise;
    generic.add_options()
    ( "version,v", "Print version string" )
    ( "help,h", "Produce help message" )
    ( "noise,n", po::value<unsigned>( &noise ) ->default_value( 2 ), "Verbosity level (0-2)" )
    ;

    // Declare groups of options that will be allowed both on command line and in config file
    po::options_description solution( "Solution strategies" );
    std::string bwalgo, solver;
    solution.add_options()
#ifdef HAVE_PARDISO
    ( "bwalgo", po::value<std::string>( &bwalgo ) ->default_value( "mmd" ),
      "Bandwidth algorithm (mmd, sloan, reverse_sloan, cuthill_mckee, reverse_cuthill_mckee)" )
    ( "solver", po::value<std::string>( &solver ) ->default_value( "pardiso" ), "Solver (pardiso, skyline)" )
#else
    ( "bwalgo", po::value<std::string>( &bwalgo ) ->default_value( "sloan" ),
      "Bandwidth algorithm (sloan, reverse_sloan, cuthill_mckee, reverse_cuthill_mckee)" )
    ( "solver", po::value<std::string>( &solver ) ->default_value( "skyline" ), "Solver (pardiso, skyline)" )
#endif
    ;

    // Files to manage
    po::options_description files( "Files" );
    files.add_options()
    ( "model-path", po::value<std::string>() ->default_value( "DucatiFrameFEModel.ansys" ),
      "Finite Element model to be used" )
    ( "res-path", po::value<std::string>() ->default_value( "DucatiFrame.objectives" ),
      "Path to result file where objective values are stored" )
    ( "param-format", po::value<std::string>() ->default_value( "DoubleList" ),
      "Format of param file to read" )
    ( "param-path", po::value<std::string>() ->default_value( "DoubleList.params" ),
      "Path to parameter file")
    ;

    // Put the different option containers together
    po::options_description options( "\nFELyX usage: program [options] model_path \nAllowed options" );
    options.add( generic ).add( solution ).add( files );

    // Store options in variables_map
    po::variables_map vm;
    store( po::parse_command_line( ac, av, options ), vm );
    notify( vm );

    // Define reactions for the different options
    if ( vm.count( "help" ) ) {
      std::cout << options << "\n";
      return 0;
    }

    if ( vm.count( "version" ) ) {
      std::cout << "\nFELyX " << VERSION << " - The Finite Element LibrarY eXperiment \n\n";
      return 0;
    }

    fs::path ModelPath( vm[ "model-path" ].as<std::string>(), fs::native );
    fs::path ResPath  ( vm[ "res-path"   ].as<std::string>(), fs::native );
    std::string ParamFormat = vm[ "param-format"   ].as<std::string>();
    fs::path ParamPath( vm[ "param-path"   ].as<std::string>(), fs::native );

    
    // RUN FEA ANALYSIS
    //StructObject FEM(ModelPath.leaf(), fs::complete(ModelPath).branch_path().native_file_string(),noise );

    DucFelyxObject FEM( noise, 15, ModelPath );
   
    FEM.NodesReordering( bwalgo );

    FEM.updateTubeProperties(ParamPath,ParamFormat);
    
    // Save new model, if needed
    // FEM.SaveAnsysModel();
    
    FEM.ApplyTorsionLoadcase();

    if ( solver == "skyline" )
      FEM.DirectSolver();
    else if ( solver == "pardiso" )
      FEM.SparseSolver();
    else
      std::cerr << "WARNING: No valid solver type specified, evaluating nothing!\n ";

    double Stiffness= FEM.EvalTorsionStiffness();

    FEM.ApplyBrakingLoadcase();

    if ( solver == "skyline" )
      FEM.DirectSolver();
    else if ( solver == "pardiso" )
      FEM.SparseSolver();
    else
      std::cerr << "WARNING: No valid solver type specified, evaluating nothing!\n ";

    double MaxStress= FEM.EvalMaxStress( 1, 2 ) ;

    double Mass= FEM.EvalMass(1,2);
    
    exportObjectives( ResPath, Mass, Stiffness, MaxStress);
     
  }

  // Catch exceptions
  catch ( std::exception & e ) {
    cerr << e.what() << endl;
    return 1;
  }

  std::cout << "FELyX analysis done\n";

  return 0;
}
