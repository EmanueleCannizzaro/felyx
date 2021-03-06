//
// C++ Implementation: spirit_exercise
//
// Description:
//
//
// Author: B.Meier, U.Mennel, O.Koenig, R.Roos, M.Wintermantel, N.Zehnder <boris.meier@geberit.com, umennel@student.ethz.ch, koenig@even-ag.ch, roos@imes.mavt.ethz.ch, wintermantel@even-ag.ch, zehnder@even-ag.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#define BOOST_SPIRIT_DEBUG

#include <iostream>
#include <fstream>
#include <string>

// Spirit headers
#include <boost/spirit.hpp>
#include <boost/spirit/core.hpp>

using namespace boost::spirit;

int main()
{

  //! Define types of SPIRIT constructs
  typedef char * iterator_t;
  typedef scanner<iterator_t> scanner_t;
  typedef rule<scanner_t> rule_t;


  std::ifstream infile( "../../tutorial/myFirstNodeList.txt" );
  if ( !infile.good() )
    std::cout << "ERROR: In ReadFileToMemory: Could not open the file " << std::endl;


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

  std::vector<double> x_vec( 0 ), y_vec( 0 ), z_vec( 0 );
  std::vector<unsigned> nodeIndex, elementIndex, firstNode, secondNode, thirdNode, fourthNode;


  // Here goes the rule to read a single line consisting of 3 comma delimited real numbers
  rule_t read_coordinates_line
  = real_p[ push_back_a( x_vec ) ] >> ch_p(',') >> *blank_p >>
    real_p[ push_back_a( y_vec ) ] >> ch_p(',') >> *blank_p >>
    real_p[ push_back_a( z_vec ) ] >> eol_p;

  BOOST_SPIRIT_DEBUG_RULE( read_coordinates_line );

  iterator_t first( buffer );
  iterator_t last = first + filesize;
  parse_info<iterator_t> info = parse( first,last, *( read_coordinates_line ) );


  //Check the input
  if ( x_vec.size() == y_vec.size() && x_vec.size() == z_vec.size() )
  {
    for (unsigned i = 0 ; i < x_vec.size() ; ++i)
      std::cout << "[ " << x_vec[i] << ", " << y_vec[i] << ", " << z_vec[i] << " ]" << std::endl;
    std::cout << std::endl;
  }


  delete[] buffer;
  return 0;
}
