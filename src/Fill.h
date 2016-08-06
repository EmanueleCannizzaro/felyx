//
// C++ Interface: Fill
//
// Description: 
//
//
// Author: Oliver Koenig, Boris Meier, Marc Wintermantel, Nino Zehnder <{okoenig,meierbo,manitou,nzehnder}@users.sourceforge.net>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//

/*
  This file is to demo how to use minimum_degree_ordering algorithm.
  
  Important Note: This implementation requires the BGL graph to be
  directed.  Therefore, nonzero entry (i, j) in a symmetrical matrix
  A coresponds to two directed edges (i->j and j->i).

  The bcsstk01.rsa is an example graph in Harwell-Boeing format,
  and bcsstk01 is the ordering produced by Liu's MMD implementation.
  Link this file with iohb.c to get the harwell-boeing I/O functions.
  To run this example, type:

  ./minimum_degree_ordering bcsstk01.rsa bcsstk01

*/
#ifndef Fill_h
#define Fill_h Fill_h

#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <numeric>

#include "ElementHeaders.h"
#include "PtrVector.h"
#include <boost/config.hpp>
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_utility.hpp"
#include "boost/graph/minimum_degree_ordering.hpp"

namespace fill_algos{

using namespace std;
using namespace fe_base;
using namespace boost;


//////////////////////////////////////////////
//// Fill Reduction - Reordering nodes 
//////////////////////////////////////////////
  template<class element_type, class node_type>
    vector<int> FillReduction(vector<node_type> &Nodes,
			    PtrVector<element_type*> &Elements,
			    string algorithm, 
			    int noise,
                            int delta)
    {

      
    //must be BGL directed graph now
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
        
    //typedef adjacency_list<vecS, vecS, directedS>  Graph;
    
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
 
    typedef std::vector<int> Vector;
    /*    
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::vertices_size_type size_type;
    */
    typedef std::pair<std::size_t, std::size_t> Pair;
    unsigned int n = Nodes.size();
    Graph G( n );
  
    std::vector<node_type> dummyNodes;
  
    // Filling up a set of pairs storing all edges 
    // of the model --> automatically avoid duplicate edges
    // Finally copy edges to Graph
    // ---------------------------------------------------------------------
  
    Pair Edge;
    std::set< Pair > Edges;
    std::set< Pair >::iterator itEdges;
  
    //vector<Node>::iterator itNode;
    typename std::vector<node_type>::iterator FirstNode, SecondNode;
    unsigned i, j, elemCount;
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
    elemCount = 0;
    // Copy set of Edges to Graph
    for ( itEdges = Edges.begin(); itEdges != Edges.end() ; ++itEdges ){
      elemCount++;
      add_edge( itEdges->first, itEdges->second, G );
      //add_edge( itEdges->second, itEdges->first, G );
    }
  
  Vector inverse_perm(n, 0);
  Vector perm(n, 0);
  
  Vector supernode_sizes(n, 1); // init has to be 1

  boost::property_map<Graph, vertex_index_t>::type 
    id = get(vertex_index, G);

  Vector degree(n, 0);
  // Switch between different algorithms
  if(algorithm == "mmd") {
    minimum_degree_ordering
    (G,
     make_iterator_property_map(&degree[0], id, degree[0]),
     &inverse_perm[0],
     &perm[0],
     make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]), 
     delta, id);
     
   } 
     
    //permutation of the nodes 
    dummyNodes = Nodes; 
    for(size_t i = 0; i < n; ++i)
      {
	Nodes[i] = dummyNodes[id[perm[i]]];
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
	    (*itElement)->SetNodeIter(i, Nodes.begin()+(id[inverse_perm[pos]]) );
	  }	
    }
    
    return perm;
  }
}

#endif

