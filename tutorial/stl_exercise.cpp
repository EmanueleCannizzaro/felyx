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

int main() {
  std::cout << "START stl_exercise\n\n";

  std::cout << "Handling vectors\n";
  std::cout << "================\n";

  // Create a vector with 10 int values, initialized with 0
  std::vector<int> v( 10, 0 );

  // Change values in vector using operator[]
  for ( unsigned i = 0; i < v.size(); ++i ) {
    v[ i ] = i * 10;
  }

  // Insert two additional values at the end of the vector, and another one at the begin
  v.push_back( 110 );
  v.push_back( 120 );


  // Resize vector to 20 ints, new values shall be initialized with 200
  v.resize( 20, 200 );

  // Print vector using iterators
  std::cout << "v=( ";
  for ( std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++ it ) {
    std::cout << *it << ", ";
  }
  std::cout << ")\n";

  // Find elements in vector
  std::cout << "Find element in vector with value = 80: ";
  std::vector<int>::iterator my_it = std::find( v.begin(), v.end(), 80 );
  std::cout << " index=" << std::distance( v.begin(), my_it );
  std::cout << " val=" << *my_it << "\n";

  std::cout << "Find (non-existing) element in vector with value = 79: ";
  my_it = std::find( v.begin(), v.end(), 79 );
  std::cout << " index=" << std::distance( v.begin(), my_it );
  std::cout << " val=" << *my_it << "\n";


  // Count the number of elements in vector with value =200
  std::cout << "# elements in vector with value = 200: " << std::count( v.begin(), v.end(), 200 ) << "\n";

  // Find max/min values in vector
  my_it = std::max_element( v.begin(), v.end() );
  std::cout << "Max element in vector: index=" << std::distance( v.begin(), my_it ) << " val=" << *my_it << "\n";

  // Add 5 to all values < 200
  for ( my_it = v.begin(); my_it != v.end(); ++my_it ) {
    if ( *my_it < 200 ) {
      *my_it += 5;
    }
  }

  // Print modified vector using copy algorithm and ostream_iterator
  std::cout << "v_mod=( ";
  std::copy( v.begin(), v.end(), std::ostream_iterator<int>( std::cout, ", " ) );
  std::cout << ")\n";
  // Sum up vector values
  int sum = std::accumulate( v.begin(), v.end(), 0 );
  std::cout << "Sum(v_mod) = " << sum << "\n";



  std::cout << "\nEND stl_exercise\n";


  return 0;
}
