//
// Standard Template Library (STL) exercise
//

// Everything needed can be found in the STL documentation on www.sgi.com/tech/stl !!

#include <iostream>


// We need some container classes
#include <vector>
#include <map>

// We need algorithms
#include <algorithm>
#include <numeric>

// New namespace to define some simplified FEM classes / structs
namespace mini_fem {

  //! Mini node data struct
  /*! A struct is a simple class containing only public data member */
  struct node {

  public:

    //! Empty constructor
    node()
        : x( 0.0 ), y( 0.0 ), z( 0.0 ) {}

    //! Constructor which initializes coordinates
    node( const double& x_, const double& y_, const double& z_ )
        : x( x_ ), y( y_ ), z( z_ ) {}

    //! Coordinates
    double x, y, z;

  };

  //! Print mini node
  std::ostream& operator<<( std::ostream& os, const node& n ) {
    os << "[" << n.x << "," << n.y << "," << n.z << "]";
    return os;
  }

  //! Equality operators
  bool operator==( const node& n1, const node& n2 ) {
    return n1.x == n2.x && n1.y == n2.y && n1.z == n2.z;
  }
  bool operator!=( const node& n1, const node& n2 ) {
    return !( n1 == n2 );
  }


}

// Main function
int main() {
  std::cout << "START stl_exercise\n\n";

  std::cout << "Dealing with node vectors\n";
  std::cout << "=========================\n";

  // Define a node vector type
  typedef std::vector< mini_fem::node > node_vector_type;

  // Create a vector with 10 nodes
  node_vector_type nodes( 5 );


  // Change values in vector using operator[]
  for ( unsigned i = 0; i < nodes.size(); ++i ) {
    nodes[ i ].x = i * 2.0;
    nodes[ i ].y = i * 4.0;
    nodes[ i ].z = i * 6.0;
  }

  // Print vector using iterators
  std::cout << "nodes=( ";
  for ( node_vector_type::const_iterator it = nodes.begin(); it != nodes.end(); ++ it ) {
    std::cout << *it << ", ";
  }
  std::cout << ")\n";

  // Insert two additional nodesat the end of the vector
  nodes.push_back( mini_fem::node( 12.0, 14.0, 16.0 ) );
  nodes.push_back( mini_fem::node( 15.0, 16.0, 18.0 ) );

  // Resize vector to 20 ints, new values shall be initialized with [100,100,100]
  nodes.resize( 10, mini_fem::node( 1.0, 1.0, 1.0 ) );


  // Find specific node in vector
  std::cout << "Find node [8,16,24] in vector: ";
  node_vector_type::iterator my_it = std::find( nodes.begin(), nodes.end(), mini_fem::node( 8.0, 16.0, 24.0 ) );
  std::cout << " index=" << std::distance( nodes.begin(), my_it );
  std::cout << " node=" << *my_it << "\n";

  std::cout << "Find (non-existing) node in vector: ";
  my_it = std::find( nodes.begin(), nodes.end(), mini_fem::node( 12.0, 1.0, 3.0 ) );
  std::cout << " index=" << std::distance( nodes.begin(), my_it );
  std::cout << " val=" << *my_it << "\n";

  // Print final node vector using copy algorithm and ostream_iterator
  std::cout << "nodes_mod=( ";
  std::copy( nodes.begin(), nodes.end(), std::ostream_iterator<mini_fem::node>( std::cout, ", " ) );
  std::cout << ")\n";


  std::cout << "\nEND stl_exercise\n";


  return 0;
}
