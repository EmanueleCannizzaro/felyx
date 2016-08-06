/***************************************************************************
*   Copyright (C) 2005 by Oliver Koenig                                   *
*   Email: <koenig@even-ag.ch>                                            *
*                                                                         *
*   This file is part of FELyX (Finite Element Library eXperiment).       *
*   See library home page at http://felyx.sourceforge.net/                *
*                                                                         *
*   FELyX is free software; you can redistribute it and/or modify         *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   FELyX is distributed in the hope that it will be useful,              *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along FELyX; if not, write to the                                     *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/

// At the moment the following function to parse ansys results as well as the whole test is implemented for a single FE problem, very brute force...
// Should be generalized to check all kinda ANSYS FEM models!!

// For accuracy reasons, ANSYS result should be printed using something like:
// -->/FORMAT,12,E,25,15,1000,1000
// -->PRNSOL,U

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//#define BOOST_SPIRIT_DEBUG
// BOOST headers
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite_ex.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/bind.hpp>
#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>
using namespace boost::unit_test;

// Spirit headers
#include <boost/spirit.hpp>
#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/clear_actor.hpp>

#include "StructObject.h"

typedef std::vector<double> node_deformationT;
typedef std::vector< node_deformationT > node_deformationsT;

// Define an actor to store values in a container
template < class valueT, class containerT>
class push_back_value {

  public:
    push_back_value( valueT& value_, containerT& container_ )
        : value( value_ ), container( container_ ) {}

    template <class IteratorT>
    void operator() ( const IteratorT first, const IteratorT last ) const {
      container.push_back( value );
    }

  private:
    valueT& value;
    containerT& container;
};

void parse_ansys_nodal_solution( const std::string& path, node_deformationsT& node_deformations ) {
  using namespace boost::spirit;

  //! Define types of SPIRIT constructs
  typedef char * iterator_t;
  typedef scanner<iterator_t> scanner_t;
  typedef rule<scanner_t> rule_t;

  std::ifstream infile( path.c_str() );
  BOOST_CHECK( infile.good() );

  // get size of file
  infile.seekg( 0, std::ifstream::end );
  unsigned long filesize = infile.tellg();
  infile.seekg( 0 );

  // allocate memory for file content
  char * buffer;    // the buffer, where the text file is stored for fast parsing
  buffer = new char [ filesize ];

  // Read content of infile to buffer
  infile.read( buffer, filesize );

  // close the infilestream, it is no longer used
  infile.close();

  // Store the deformations in a simple vector
  node_deformationT node_deformation;

  // Rule to read nodal displacements
  rule_t nodal_displacement_r = *blank_p >> uint_p
                                >> repeat_p( 3, 6 ) [ *blank_p >> real_p[ push_back_a( node_deformation ) ] ]
                                >> *( anychar_p - eol_p )
                                >> eol_p
                                [ push_back_value<node_deformationT, node_deformationsT>( node_deformation, node_deformations ) ]
                                [ clear_a( node_deformation ) ];
  BOOST_SPIRIT_DEBUG_RULE( nodal_displacement_r );

  iterator_t first( buffer );
  iterator_t last = first + filesize;
  parse_info<iterator_t> info = parse( first, last, *( nodal_displacement_r | ( *( anychar_p - eol_p ) >> eol_p ) ) );
  BOOST_CHECK( info.full );
  delete[] buffer;

}

void linear_struct_analysis_test( const std::string& dir, const std::string& fname, const double& percent_tolerance ) {

  BOOST_MESSAGE( "START analysis" );

  felyx::StructObject FEM( fname + ".ansys", dir, 2 );
  FEM.SparseSolver();
  //FEM.DirectSolver();
  FEM.PrintGlobalStatus();
  FEM.SaveResults( fname + ".felyxres" );

  // Read ansys nodal solution for this model
  BOOST_MESSAGE( "Compare evaluated deformations with ANSYS results" );
  node_deformationsT ansys_deformations;
  parse_ansys_nodal_solution( dir + "/" + fname + ".ansysres", ansys_deformations );

  // Check nodal deformations
  std::vector<felyx::StructObject::node_type>& nodes = FEM.GetNodes();
  BOOST_CHECK_EQUAL( nodes.size(), ansys_deformations.size() );

  for ( unsigned i = 0; i < ansys_deformations.size(); ++i ) {
    Dense_Vector felyx_deformation = nodes[ i ].GetDeformations();
    for ( unsigned j = 0; j < 3; ++j ) {
      // Use check_is_close instead of BOOST_CHECK_CLOSE to customize error message printed
      if ( !boost::test_tools::check_is_close( ansys_deformations[ i ][ j ], felyx_deformation[ j ], percent_tolerance ) ) {
        BOOST_ERROR( "Comparison of nodal deformation of DOF "
                    + boost::lexical_cast<std::string>( j ) + " for node "
                    + boost::lexical_cast<std::string>( i + 1 ) + " failed: "
                    + boost::lexical_cast<std::string>( ansys_deformations[ i ][ j ] ) + " != "
                    + boost::lexical_cast<std::string>( felyx_deformation[ j ]) + "\n" );
      }
    }

  }
  BOOST_MESSAGE( "Comparison done." );
}

test_suite*
init_unit_test_suite( int argc, char* argv[] ) {

  test_suite * test = BOOST_TEST_SUITE( "Test FELyX linear structural analysis" );

  std::string dir = "../../fedata";
  std::vector<std::string> fnames;
  std::vector<double> percent_tolerances;

  fnames.push_back( std::string( "StructModel-PrescribedDisplacements" ) );
  percent_tolerances.push_back( 1e-5 );
  fnames.push_back( std::string( "StructModel-Winschhebel" ) );
  percent_tolerances.push_back( 1e-3 );

  // The 3Dmix model seems to reveil some errors, results in FELYX and ANSYS are different.
  // --> Should be investigated!!
  //fnames.push_back( std::string( "StructModel-3Dmix-1131e" ) );
  //percent_tolerances.push_back( 1 );

  for ( unsigned i = 0; i < fnames.size(); ++i ) {
    test->add
    ( BOOST_TEST_CASE( boost::bind( &linear_struct_analysis_test, dir, fnames[ i ], percent_tolerances[ i ] ) ) );
  }

  return test;

}
