//-----------------------------------------------------------------------------
// Bandwidth.h
//
// begin     : Dec 6 2001
// copyright : (c) 2001 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {okoenig, wintermantel, nzehnder}@imes.mavt.ethz.ch
// www       : www.structures.ethz.ch
/* 
   This file is part of FELyX (Finite Element Library eXperiment).
   
   FELyX is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   FELyX is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with FELyX; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
//-----------------------------------------------------------------------------
   
#ifndef Bandwidth_h
#define Bandwidth_h Bandwidth_h

#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <numeric>

#include "ElementHeaders.h"
#include "PtrVector.h"

#include "boost/config.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/sloan_ordering.hpp"
#include "boost/graph/cuthill_mckee_ordering.hpp"
#include "boost/graph/properties.hpp"
#include "boost/graph/bandwidth.hpp"
#include "boost/graph/profile.hpp"

using namespace std;
using namespace fe_base;
using namespace boost;

namespace fe_base{


//////////////////////////////////////////////
//// Bandwith Reduction - Reordering nodes 
//////////////////////////////////////////////
  template<class element_type, class node_type>
    vector<unsigned> BandwidthReduction(vector<node_type> &Nodes,
			    PtrVector<element_type*> &Elements,
			    string algorithm, 
			    int noise = 1,
			    double sloan_weight1 = 1.0 ,
			    double sloan_weight2 = 2.0)
    {

      typedef adjacency_list<
	vecS, 
	vecS, 
	undirectedS, 
	property<
	vertex_color_t, 
	default_color_type,
	property<
	vertex_degree_t,
	int,
	property<
	vertex_priority_t,
	double > > > > Graph;
      
  typedef graph_traits<Graph>::vertex_descriptor Vertex;
  typedef graph_traits<Graph>::vertices_size_type size_type;
  typedef std::pair<std::size_t, std::size_t> Pair;

  Graph G( Nodes.size() );
  
  int orgBandwidth, newBandwidth, orgProfile, newProfile = 0;
  bool improved = false;
  vector<node_type> dummyNodes;
  
  // Filling up a set of pairs storing all edges 
  // of the model --> automatically avoid duplicate edges
  // Finally copy edges to Graph
  // ---------------------------------------------------------------------
  
  Pair Edge;
  std::set< Pair > Edges;
  std::set< Pair >::iterator itEdges;
  
  //vector<Node>::iterator itNode;
  typename vector<node_type>::iterator FirstNode, SecondNode;
  unsigned i, j;
  typename PtrVector<element_type*>::iterator itElement;
  
  // loop through elements and the nodes of each element to find all edges
  for(itElement = Elements.begin(); itElement!=Elements.end(); itElement++) {   
    for(i = 0; i < (*itElement)->GetNodeCount(); ++i) {
      for(j = i + 1; j < (*itElement)->GetNodeCount(); ++j ) {
	
	// Find indices of nodes of this certain edge
	FirstNode     = (*itElement)->GetNodeIter(i);
	Edge.first    = distance( Nodes.begin(), FirstNode  );
	SecondNode    = (*itElement)->GetNodeIter(j);
	Edge.second   = distance( Nodes.begin(), SecondNode );

	// Ensure that smaller index is stored in PosFirstNode
	if ( Edge.first > Edge.second )
	  std::swap ( Edge.first, Edge.second );
	
	// Insert pair into set
	Edges.insert( Edge );
      }
    }
  }

  // Copy set of Edges to Graph
  for ( itEdges = Edges.begin(); itEdges != Edges.end() ; ++itEdges ){
    //cout << itEdges->first << " \t " << itEdges->second << endl;
    add_edge( itEdges->first, itEdges->second  ,G );
  }

  // -------------------------------------------------------------------------

  //Creating two iterators over the vertices
  graph_traits<Graph>::vertex_iterator ui, ui_end;

  //Creating a property_map with the degrees of the degrees of each vertex
  property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
  for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
    deg[*ui] = degree(*ui, G);
  
  //Creating a property_map for the indices of a vertex
  property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, G);
  
  //Creating a property_map for the the current_degree of a vertex
  //property_map<Graph, vertex_current_degree_t>::type current_degree_map = get(vertex_current_degree, G);
  
  //Creating a property_map for the the priority of a vertex
  // property_map<Graph, vertex_priority_t>::type priority_map = get(vertex_priority, G);
  
  //Calculating and printing the original bandwidth
  orgBandwidth = bandwidth(G);
  orgProfile = profile(G);

  if (noise >0){
  std::cout << "\t \t - Original bandwidth  : " << bandwidth(G) << std::endl;
  std::cout << "\t \t - Original profile    : " << profile(G) << std::endl;  
  }

  std::vector<Vertex> inv_perm(num_vertices(G));
  std::vector<unsigned> perm(num_vertices(G));
  
    
  // Switch between different algorithms
  if(algorithm == "cuthill_mckee")
    {
      //reverse cuthill_mckee_ordering
      cuthill_mckee_ordering(G, inv_perm.begin(), 
			     get(vertex_color, G), 
			     make_degree_map(G)
			     );
      
      //calculating the permutation vector
      for (size_type c = 0; c != inv_perm.size(); ++c)
	perm[index_map[inv_perm[c]]] = c;
      
      //calculating bandwidth and profile
      newBandwidth = bandwidth(G, make_iterator_property_map(&perm[0], index_map));
      newProfile = profile(G, make_iterator_property_map(&perm[0], index_map));

      if(newBandwidth < orgBandwidth) improved = true;
    }
  
 
  
  else if( algorithm == "sloan" )
    {
      //sloan ordering 
      sloan_ordering(G, inv_perm.begin(), 
	 	     get(vertex_color, G), 
		     make_degree_map(G), 
		     get(vertex_priority, G), 
		     sloan_weight1, 
		     sloan_weight2
		     );
      
      //calculating the permutation vector
      for (size_type c = 0; c != inv_perm.size(); ++c)
	perm[index_map[inv_perm[c]]] = c;
      
      //calculating bandwidth and profile
      newBandwidth = bandwidth(G, make_iterator_property_map(&perm[0], index_map));    
      newProfile = profile(G, make_iterator_property_map(&perm[0], index_map));

      if(newProfile < orgProfile) improved = true;
    }
  

  else if(algorithm == "reversed_cuthill_mckee")
    {
      //reverse cuthill_mckee_ordering
      cuthill_mckee_ordering(G, inv_perm.rbegin(), 
			     get(vertex_color, G), 
			     make_degree_map(G)
			     );
      
      //calculating the permutation vector
      for (size_type c = 0; c != inv_perm.size(); ++c)
	perm[index_map[inv_perm[c]]] = c;
      
      //calculating bandwidth and profile
      newBandwidth = bandwidth(G, make_iterator_property_map(&perm[0], index_map));
      newProfile = profile(G, make_iterator_property_map(&perm[0], index_map));

      if(newBandwidth < orgBandwidth) improved = true;
    }
    
  else if(algorithm == "reversed_sloan" )
    {
      //sloan ordering 
      sloan_ordering(G, inv_perm.rbegin(), 
		     get(vertex_color, G), 
		     make_degree_map(G), 
		     get(vertex_priority, G), 
		     sloan_weight1, 
		     sloan_weight2 
		     );
      
      //calculating the permutation vector
      for (size_type c = 0; c != inv_perm.size(); ++c)
	perm[index_map[inv_perm[c]]] = c;
      
      //calculating the bandwidth
      newBandwidth = bandwidth(G, make_iterator_property_map(&perm[0], index_map));
      newProfile = profile(G, make_iterator_property_map(&perm[0], index_map));

      if(newProfile < orgProfile) improved = true;
    }
  
    else {
      cout << endl << "\t WARNING in Bandwidth.cc : " << endl;
      cout << "\t - No valid bandwidth algorithm specified - no optimization is done !" << endl << endl;
      newBandwidth = orgBandwidth;
    }

    if(improved)
      {
	if (noise > 0 ){
	cout << "\t\t - Optimized bandwidth : " << newBandwidth << endl;
	cout << "\t\t - Optimized profile   : " << newProfile << endl;
	}
	
	// Permutation of the Node-List
	dummyNodes=Nodes;
	for(size_t i = 0; i < Nodes.size(); ++i)
	  {
	    Nodes[i] = dummyNodes[index_map[inv_perm[i]]];
	  }
	
	// Permutation of the Nodes in the Elements
	typename vector<node_type>::iterator NodeIter;
	size_t pos;
	
	for(itElement = Elements.begin(); itElement!=Elements.end(); itElement++)
	  {
	    for(unsigned i = 0; i < (*itElement)->GetNodeCount(); ++i)
	      {
		NodeIter = (*itElement)->GetNodeIter(i);
		pos = distance(Nodes.begin(), NodeIter );
		(*itElement)->SetNodeIter(i, Nodes.begin()+(perm[pos]));
	      }	
	  }
	
      }
    else
      cout << "\t - Original ordering better than new one!" << endl;
    return perm;
}
  
} // of namespace

#endif

