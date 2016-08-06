//
// Implementation file for MyStructObject
//

#include "MyStructObject.h"

// Implementation of function to create FE model
void MyStructObject::create_model()
{

  std::cout << "Create FE model\n";

  // Data of rectangle to create
  // ---------------------------
  double length = 100.0;
  unsigned n_ele_length = 5;

  double width = 30.0;
  unsigned n_ele_width = 3;

  unsigned n_ele = n_ele_length * n_ele_width;
  unsigned n_nodes = ( n_ele_length + 1 ) * ( n_ele_width + 1 );

  // Create nodes of model
  // ----------------------

  // Resize node container
  Nodes.resize( n_nodes );

  // We need a node iterator
  std::vector<node_type>::iterator nit = Nodes.begin();

  // Fill in coordinates of nodes
  for ( unsigned j = 0; j < n_ele_width + 1; ++j )
  {
    for ( unsigned i = 0; i < n_ele_length + 1; ++i )
    {

      // Set coords of node
      nit->set
      ( i * length / n_ele_length, j * width / n_ele_width, 0.0 );

      // Go to next node
      ++nit;
    }
  }

  // Create a material
  // -----------------
  Materials.push_back( new fe_base::IsotropicMaterial() );
  Materials[ 0 ] ->Set( "E", 70000 );
  Materials[ 0 ] ->Set( "nu", 0.3 );

  // Create elements of model
  // ------------------------
  //Resize element ptrvector
  Elements.resize( n_ele );

  // We need a element container iterator
  PtrVector<element_type*>::iterator ele_it = Elements.begin();
  unsigned n0, n1, n2, n3;

  // Fill in the nodes of the elements
  for ( unsigned j = 0; j < n_ele_width; ++j )
  {

    n0 = j * ( n_ele_length + 1 );

    for ( unsigned i = 0; i < n_ele_length; ++i )
    {

      // Create new plane element
      ( *ele_it ) = new fe_base::Plane182;

      // Compute node indices
      n1 = n0 + 1;
      n2 = n0 + n_ele_length + 2;
      n3 = n0 + n_ele_length + 1;

      // Add node iterator to plane elements
      ( *ele_it ) ->SetNodeIter( 0, Nodes.begin() + n0 );
      ( *ele_it ) ->SetNodeIter( 1, Nodes.begin() + n1 );
      ( *ele_it ) ->SetNodeIter( 2, Nodes.begin() + n2 );
      ( *ele_it ) ->SetNodeIter( 3, Nodes.begin() + n3 );

      // Set material
      ( *ele_it ) ->SetMaterialPtr( Materials[ 0 ] );

      // Increment node and element indices/iterators
      ++ele_it;
      ++n0;

    }
  }

  // Create and assign boundary conditions
  // ---------------------------
  BoundaryConditions.resize( 3 );

  BoundaryConditions[ 0 ].set( Dx, 0.0 );
  BoundaryConditions[ 0 ].set( Dy, 0.0 );

  Nodes[ 0 ].set( &BoundaryConditions[ 0 ] );
  Nodes[ n_ele_width * ( n_ele_length + 1 ) ].set( &BoundaryConditions[ 0 ] );
  BoundaryConditions[ 1 ].set( Fx, 100.0 );
  BoundaryConditions[ 2 ].set( Fy, 100.0 );
  Nodes[ n_ele_length ].set( &BoundaryConditions[ 1 ] );
  Nodes[ n_nodes - 1 ].set( &BoundaryConditions[ 2 ] );

}

void MyStructObject::create_model_from_file(std::string fname_, std::string fpath_)
{

  using namespace boost::spirit;

  //! Define types of SPIRIT constructs
  typedef char * iterator_t;
  typedef scanner<iterator_t> scanner_t;
  typedef rule<scanner_t> rule_t;


  // The file path
  std::string filepath((fpath_ + "/" + fname_).c_str());

  std::ifstream infile(filepath.c_str());

  if (!infile.good())
  {
    std::cout << "ERROR: In ReadFileToMemory: Could not open the file " << std::endl;
      }

  // get size of file
  infile.seekg(0,std::ifstream::end);
  unsigned long filesize=infile.tellg();
  infile.seekg(0);

  // allocate memory for file content
  char * buffer;    // the buffer, where the text file is stored for fast parsing
  buffer = new char [filesize];

  // Read content of infile to buffer
  infile.read(buffer,filesize);

  // close the infilestream, it is no longer used
  infile.close();

  std::vector<double> x_vec(0), y_vec(0), z_vec(0);
  std::vector<unsigned> nodeIndex, elementIndex, firstNode, secondNode, thirdNode, fourthNode;


  // Here goes the rule to read a single line consisting of 3 comma delimited real numbers
  rule_t read_coordinates_line
  = uint_p[ push_back_a( nodeIndex ) ] >> blank_p >>
    real_p[ push_back_a( x_vec ) ] >> blank_p >>
    real_p[ push_back_a( y_vec ) ] >> blank_p >>
    real_p[ push_back_a( z_vec ) ] >> eol_p;

  BOOST_SPIRIT_DEBUG_RULE( read_coordinates_line );

  rule_t read_elements_line
  = uint_p[ push_back_a( elementIndex ) ] >> blank_p >>
    uint_p[ push_back_a( firstNode ) ] >> blank_p >>
    uint_p[ push_back_a( secondNode ) ] >> blank_p >>
    uint_p[ push_back_a( thirdNode ) ] >> blank_p >>
    uint_p[ push_back_a( fourthNode ) ] >> eol_p;
  BOOST_SPIRIT_DEBUG_NODE( read_elements_line );

  rule_t comment_parser = comment_p('#'); //ch_p( '#' ) >> *( anychar_p - eol_p ) >> eol_p ;
  BOOST_SPIRIT_DEBUG_NODE( comment_parser );

  rule_t read_any_line = read_coordinates_line | read_elements_line | comment_parser ;
  BOOST_SPIRIT_DEBUG_NODE( read_any_line );

  iterator_t first( buffer );
  iterator_t last = first + filesize;
  unsigned long psize = parse( first,last, *(read_coordinates_line | read_elements_line | comment_parser |eol_p ) ).length;



  if ( psize < filesize ) {
    std::cout <<  "WARNING: Not the entire input file " << filepath << " has been parsed." << std::endl;
    std::cout <<  "Only " << psize << " characters of total " << filesize << " are read" << std::endl;
  }
  

  // Create a material
  // -----------------
  Materials.push_back( new fe_base::IsotropicMaterial() );
  Materials[ 0 ] ->Set( "E", 70000 );
  Materials[ 0 ] ->Set( "nu", 0.3 );

  // Fill the nodes of the model
  // ----------------------

  // Resize node container
  Nodes.resize( 24 );

  // Fill in coordinates of nodes
    for ( unsigned i = 0; i < 24 ; ++i )
    {
      // Set coords of node
      Nodes[nodeIndex[i]].set( x_vec[i], y_vec[i], z_vec[i] );
    }


  // Create elements of model
  // ------------------------
  //Resize element ptrvector
  Elements.resize( firstNode.size() );

  // We need a element container iterator
  PtrVector<element_type*>::iterator ele_it = Elements.begin();

  // Fill in the nodes of the elements
    for ( unsigned i = 0; i < Elements.size() ; ++i )
    {
      // Create new plane element
      ( *ele_it ) = new fe_base::Plane182;

      // Add node iterator to plane elements
      ( *ele_it ) ->SetNodeIter( 0, Nodes.begin() + firstNode[i] );
      ( *ele_it ) ->SetNodeIter( 1, Nodes.begin() + secondNode[i] );
      ( *ele_it ) ->SetNodeIter( 2, Nodes.begin() + thirdNode[i] );
      ( *ele_it ) ->SetNodeIter( 3, Nodes.begin() + fourthNode[i] );

      // Set material
      ( *ele_it ) ->SetMaterialPtr( Materials[ 0 ] );

      ++ele_it;

    }

  // Create and assign boundary conditions
  // ---------------------------
  BoundaryConditions.resize( 3 );

  BoundaryConditions[ 0 ].set( Dx, 0.0 );
  BoundaryConditions[ 0 ].set( Dy, 0.0 );

  Nodes[ 0 ].set( &BoundaryConditions[ 0 ] );
  Nodes[ 18 ].set( &BoundaryConditions[ 0 ] );
  BoundaryConditions[ 1 ].set( Fx, 100.0 );
  BoundaryConditions[ 2 ].set( Fy, 100.0 );
  Nodes[ 5 ].set( &BoundaryConditions[ 1 ] );
  Nodes[ 23 ].set( &BoundaryConditions[ 2 ] );

  //free the memory allocated to buffer
  delete[] buffer;
}


// Implementation of function to print out model data
void MyStructObject::print_model()
{

  std::cout << "Node list:\n";
  unsigned index = 0;
  for ( std::vector<node_type>::const_iterator nit = Nodes.begin(); nit != Nodes.end(); ++nit, ++index )
  {
    std::cout << index << " " << *nit << "\n";
  }

  std::cout << "Element list:\n";
  index = 0;
  for ( PtrVector<element_type*>::const_iterator eit = Elements.begin(); eit != Elements.end(); ++eit, ++index )
  {
    std::cout << index << " " << **eit << "\n";
  }

  std::cout << "Boundary condition list:\n";
  index=0;
  for ( std::vector<bc_type>::const_iterator bit = BoundaryConditions.begin(); bit!=BoundaryConditions.end(); ++bit, ++index)
  {
    std::cout << index << " " << *bit << "\n";
  }

}

void MyStructObject::write_nodes_and_elements(std::string fname_, std::string fpath_)
{

  ofstream outfile;
  outfile.open((fpath_+ "/" + fname_).c_str());
  if (outfile.is_open())
  {
    for (unsigned i = 0 ; i < Nodes.size() ; ++i )
    {
      outfile.setf(std::ios::showpoint);
      outfile << i << " " << Nodes[i].GetCx() << " "
      << Nodes[i].GetCy() << " "
      << Nodes[i].GetCz() << std::endl;
    }
    outfile << std::endl;

    for (unsigned k = 0 ; k < Elements.size() ; ++k )
    {
      outfile << k ;
      for (unsigned p = 0 ; p < Elements[k]->NodeVec.size() ; ++p )
      {
        std::vector<node_type>::iterator n_begin = Nodes.begin();
        outfile << " " << distance(n_begin, Elements[k]->GetNodeIter(p));
      }
      outfile << std::endl;
    }

    outfile.close();
  }
  else
  {
    cout << "Error opening file";
  }

}

void MyStructObject::write_max_deformation(std::string fname_, std::string fpath_)
{

  ofstream outfile;
  outfile.open((fpath_+ "/" + fname_).c_str());
  if (outfile.is_open())
  {
    outfile << "The Max Deformation is : " << GetMaximumDeformation(Nodes) << std::endl;
    outfile.close();
  }
  else
  {
    cout << "Error opening file";
  }

}



