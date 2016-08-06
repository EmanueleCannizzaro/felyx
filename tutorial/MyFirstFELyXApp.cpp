//
// MyFELyXApp.cpp
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

// FELyX Main header
#include "StructObject.h"
#include "export_tecplot.hpp"

using namespace std;
using namespace felyx;

int main() {

  try {

    std::cout << "START MyFELyXApp\n";

    StructObject FEM( "StructModel-Winschhebel.ansys", "../../fedata", 1 );

    FEM.SparseSolver();

    FEM.EvalStresses();

    //export_nodes("/home/nino/projects/felyx/tutorial/felyx-winschhebel.dat", FEM.GetNodes() );

    export_tetrahedron("../../tutorial/felyx-winschhebel.dat", FEM.GetNodes(), FEM.GetElements() );

    FEM.PrintGlobalStatus();

    FEM.SaveResults();

    std::cout << "END MyFELyXApp\n";

    // Catch exceptions
  } catch ( std::exception & e ) {
    cerr << e.what() << endl;
  }

  return 0;
}

