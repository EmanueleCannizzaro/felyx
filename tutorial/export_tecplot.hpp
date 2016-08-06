 ///
/// Export meshes and results to TECPLOT files
///
#ifndef FELYX_EXPORT_TECPLOT
#define FELYX_EXPORT_TECPLOT

#include <string>
#include <fstream>
#include <sstream>
#include <map>

#include "FelyxException.h"
#include "mtl/mtl_felyx_utils.h"

namespace felyx {

  //! Export nodal results
  template <class NodeContainerT>
  void export_nodes( const std::string& fpath, const NodeContainerT& nodes ) {

    std::cout << "Export TecPlot file " + fpath + "\n";

    std::ostringstream os;
    os << "title = \"FELyX analysis results\"\n";
    os << "variables = \"x\", \"y\", \"z\",\"ux\",\"uy\",\"uz\",\"usum\",\"vonMises\"\n";
    os << "ZONE T = \"NODAL SOLUTION\"\n";
    for ( typename NodeContainerT::const_iterator it = nodes.begin(); it != nodes.end(); ++it ) {
      os << it->GetCx() << " " << it->GetCy() << " " << it->GetCz() << " ";
      mtl::Dense_Vector def = it->GetTranslations();
      os << def[ 0 ] << " " << def[ 1 ] << " " << def[ 2 ] << " " << mtl::two_norm( def ) << " ";
      os << it->GetVMStress() << "\n";
    }

    // Print everything to file
    std::ofstream outfile( fpath.c_str() );
    if ( !outfile.is_open() ) {
      std::string error = "File " + fpath + " can not be opened for writing.";
      FELYX_RUNTIME_THROW( error.c_str() );
    }
    outfile << os.str();
    outfile.close();

    std::cout << "finished\n";

  }

//! Export element results
template <class NodeContainerT, class ElementContainerT>
void export_tetrahedron(const std::string& fpath,
                        NodeContainerT& nodes,
                        const ElementContainerT& elements ) {

    std::cout << "Export TecPlot file " + fpath + "\n";

    std::ostringstream os;
    os << "title = \"FELyX analysis results\"\n";
    os << "variables = \"x\", \"y\", \"z\",\"ux\",\"uy\",\"uz\",\"usum\",\"vonMises\"\n";
    os << "ZONE T = \"NODAL SOLUTION\", N=" << nodes.size() << ", E=" << elements.size() << ", ET=TETRAHEDRON, F=FEPOINT" << std::endl;
    for ( typename NodeContainerT::const_iterator it = nodes.begin(); it != nodes.end(); ++it ) {
      os << it->GetCx() << " " << it->GetCy() << " " << it->GetCz() << " ";
      mtl::Dense_Vector def = it->GetTranslations();
      os << def[ 0 ] << " " << def[ 1 ] << " " << def[ 2 ] << " " << mtl::two_norm( def ) << " ";
      os << it->GetVMStress() << "\n";
    }



    typename NodeContainerT::iterator nodes_begin = nodes.begin();

    for ( typename ElementContainerT::const_iterator eleiter = elements.begin(); eleiter != elements.end(); ++eleiter){
        for (unsigned i = 0 ; i < 4 ; ++i)
           os <<  ( distance(nodes_begin,(*eleiter)->GetNodeIter(i)) + 1 ) << " ";
         os << "\n";
    }


    // Print everything to file
    std::ofstream outfile( fpath.c_str() );
    if ( !outfile.is_open() ) {
      std::string error = "File " + fpath + " can not be opened for writing.";
      FELYX_RUNTIME_THROW( error.c_str() );
    }
    outfile << os.str();
    outfile.close();

    std::cout << "finished\n";


}


}
#endif
