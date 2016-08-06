//
// Header file for MyStructObject
//

#ifndef MyStructObject_h
#define MyStructObject_h MyStructObject_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//#define BOOST_SPIRIT_DEBUG

#include <fstream>

#include "boost/spirit.hpp"
#include <boost/spirit/core.hpp>
#define RULE(name, definition) typeof(definition) name = definition

// include FELyX object
#include "StructObject.h"

//namespace bs = boost::spirit;
using namespace std;
using namespace felyx;

//! Definition of class MyStructObject
/*! Class is derived from StructObject */
class MyStructObject : public StructObject {

  public:

    //! Empty constructor, initialize StructObject with noise parameter
    MyStructObject() : StructObject(2) {}

    //! Create FE-model from scratch
    void create_model();
    
    //! Read nodes coordinates and elements from a file
    void create_model_from_file(std::string , std::string);

    //! Print model data
    void print_model();

    //! Write the max_deformation to a file
    void write_max_deformation(std::string , std::string);
     
    //! Write the nodes coordinates and the coincidence table to a file
    void write_nodes_and_elements(std::string , std::string);



};



#endif
