/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// PreProcessing.h
//
// begin     : Novembre 8 2002
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

#ifndef SelectionSet_h
#define SelectionSet_h SelectionSet_h

//   Example
//   -------

//   cout << "sizevec " << Nodes.size() << endl;
//   SelectionSet<LcmNode> SelectedNodes(Nodes);
//   cout << "bla " << SelectedNodes[1]->Cx << endl;
//   SelectedNodes.ReselectByLoc(x,-0.02,0);
//   SelectedNodes.ReselectByLoc(y,0,0.02);
//   cout << "size " << SelectedNodes.size() << endl;
//   SelectedNodes.print();

//   SelectionSet<LcmNode> SelectedNodes2(Nodes);
//   SelectedNodes2.ReselectByLoc(x,-0.01,0.05);
//   SelectedNodes2.ReselectByLoc(y,0.01,0.05);
//   cout << "size " << SelectedNodes2.size() << endl;
//   SelectedNodes2.print();

//   SelectedNodes.subtract( SelectedNodes2 );
//   cout << "size " << SelectedNodes.size() << endl;
//   SelectedNodes.print();

//   SelectionSet<LcmElement> SelectedElements(Elements);
//   SelectionSet<LcmNode> SelectedNodes4(Nodes);
//   SelectedElements.ReselectByLoc(x,-0.04,0);
//   SelectedElements.ReselectByLoc(y,0,0.04);
//   cout << "size " << SelectedElements.size() << endl;
//   SelectedElements.ReselectAttachedTo(SelectedNodes2);
//   cout << "size " << SelectedElements.size() << endl;

#include <iostream>
#include <cstdlib>
#include <vector>
#include <iterator>

using namespace std;

namespace fe_base{		// Put classes into namespace fe_base

template<class Entity>
class SelectionSet : public std::vector<Entity*>{

 public:
 
  enum Coordinate {x, y, z};

  SelectionSet( ) : std::vector<Entity*>( ){}

  SelectionSet( std::vector<Entity>& container_ ) {
    typename std::vector<Entity>::iterator it = container_.begin();
    while ( it != container_.end() ) {
      push_back( &(*it) ); ++it; }
    sort(this->begin(), this->end());
  }

  SelectionSet( std::vector<Entity*>& container_ ) {
    typename std::vector<Entity*>::iterator it = container_.begin();
    while ( it != container_.end() ) {
      push_back( *it ); ++it; }
    sort(this->begin(), this->end());
  }

  SelectionSet( SelectionSet<Entity>& other_ ){
    *this = other_;
  }

  void ReselectByLoc( Coordinate c, double low, double high ) {
    typename std::vector<Entity*>::iterator it = this->begin();
    
    switch (c)
      {
      case x : {
	while ( it != this->end() ) {
	  if ( low>(*it)->GetCx() || high<(*it)->GetCx() ) it=erase(it);
	  else ++it;
	}
      }
      case y : {
	while ( it != this->end() ) {
	  if ( low>(*it)->GetCy() || high<(*it)->GetCy() ) it=erase(it);
	  else ++it;
	}
      } 
      case z : {
	while ( it != this->end() ) {
	  if ( low>(*it)->GetCz() || high<(*it)->GetCz() ) it=erase(it);
	  else ++it;
	}
      } 
      }
  }

  void ReselectSphere( double x, double y, double z, double radius ) {
    typename std::vector<Entity*>::iterator it = this->begin();
    while ( it != this->end() ) {
      if ( ((*it)->Cx-x)*((*it)->Cx-x) + ((*it)->Cy-y)*((*it)->Cy-y) + ((*it)->Cz-z)*((*it)->Cz-z)
	   > radius*radius ) it=erase(it);
      else ++it;
    }
  }

  void add(SelectionSet& x) {
    SelectionSet<Entity> temp;
    insert_iterator<SelectionSet<Entity> > temp_ins(temp, temp.begin());
    set_union(this->begin(), this->end(), x.begin(), x.end(), temp_ins );
    *this = temp;
  }

  void subtract(SelectionSet& x) {
    SelectionSet<Entity> temp;
    insert_iterator<SelectionSet<Entity> > temp_ins(temp, temp.begin());
    set_difference(this->begin(), this->end(), x.begin(), x.end(), temp_ins );
    *this = temp;
  }

  void intersect(SelectionSet& x) {
    SelectionSet<Entity> odd;
    insert_iterator<SelectionSet<Entity> > odd_ins(odd, odd.begin());
    set_intersection(this->begin(), this->end(), x.begin(), x.end(), odd_ins);
    *this = odd;
  }

  // calling entity must be element
  void ReselectAttachedTo( SelectionSet<LcmNode>& x ) {
    typename std::vector<LcmNode*>::iterator it2;
    typename std::vector<Entity*>::iterator it = this->begin();
    while ( it != this->end() ) {
      bool del_it = false;
      for (unsigned i=0; i<(*it)->GetNodeCount(); ++i) {
	it2=find(x.begin(), x.end(), (*it)->(*GetNodeIter(i)));
	if (it2 != x.end()) del_it = true;
      }
      if ( del_it ) it=erase(it);
      else ++it;
    }
  }
  
  // calling entity must be node
  void print() {
    typename std::vector<Entity*>::iterator it = this->begin();
    while ( it != this->end() ) {
      cout << "Node " << **it << endl;
      ++it;
    }
  }

bool compare(Entity* a, Entity* b) {
  if (a->Cz == b->Cz) {
    if (a->Cy == b->Cy) {
      return (a->Cx < b->Cx);
    }
    else
      return (a->Cy < b->Cy);
  }
  else
    return (a->Cz < b->Cz);
}

}; // end of class

} // end of namespace fe_base

#endif
