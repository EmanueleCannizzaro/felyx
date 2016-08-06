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

//
// Template for future test programs in FELyX
//
// Using on BOOST unit test framework
// Documentation of BOOST test library:
// Main page:               http://www.boost.org/libs/test/doc/index.html
// Quick start docu:        http://www.boost.org/libs/test/doc/components/utf/getting_started/index.html
// Reference of test tools: http://www.boost.org/libs/test/doc/components/test_tools/reference/index.html
 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
using namespace boost::unit_test;

void TestTemplate() {
  
  BOOST_CHECK( "bli"=="bla");
  BOOST_CHECK_EQUAL( 1.234, 1.2345 );
  BOOST_CHECK_CLOSE( 1.234, 1.2345, 0.001 ); 
  BOOST_CHECK_CLOSE( 1.234, 1.2345, 0.1 ); 
  //...
}

test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    test_suite* test= BOOST_TEST_SUITE( "Template for FELyX tests using BOOST unit test framework" );

    test->add( BOOST_TEST_CASE( &TestTemplate ) );

    return test;
}
