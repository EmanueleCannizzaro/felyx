//
// Application to instanciate and call MyStructObject
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// FELyX Main header
#include "MyStructObject.h"

int main( int ac, char* av[] ) {

  try {

    std::cout << "FE analysis with MyStructObject\n";

    // Instanciate MyStructObject
    MyStructObject FEM;

    // Create model
    //FEM.create_model();
    FEM.create_model_from_file("SimpleFEmodel.txt","../../tutorial");

    //FEM.print_model();

    FEM.SaveAnsysModel("RectangleFEModel.ansys", "../../tutorial");

    // Solve FE problem
    FEM.SparseSolver();
  
    FEM.PrintGlobalStatus();
  
    FEM.write_max_deformation("MaxDeformation.txt" , "../../tutorial");
  }

  // Catch exceptions
  catch ( std::exception & e ) {
    cerr << e.what() << endl;
    return 1;
  }

  std::cout << "FELyX analysis done\n";

  return 0;
}

