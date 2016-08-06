//-----------------------------------------------------------------------------
// IOSpiritUtils.h
//
// begin     : Nov 2004
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
#ifndef IOSpiritUtils_h
#define IOSpiritUtils_h IOSpiritUtils_h

#include <vector>
#include <iostream>
#include <fstream>

#include "FelyxException.h"

#include "boost/spirit/core.hpp"
#include <boost/spirit/actor/ref_value_actor.hpp>

namespace fe_base{

  /*!
      Function to read the file content to the memory
      The first argument is the filepath and the second argument
      the char array where to store the stream.
  */
  unsigned long ReadFileToMemory(const char* filepath, char* &buffer);

  /*!
      An actor to reserve enough memory for a vector
  */
  struct reserve_vector_space
  {
    template <class vecT>
    void act( vecT& vec_, unsigned const& value_) const
    {
        vec_.reserve(value_);
    }
  };
  template<class vecT>
  inline boost::spirit::ref_value_actor<vecT,reserve_vector_space>
  reserve_vector_space_a(vecT& ref_)
  {
     return boost::spirit::ref_value_actor<vecT,reserve_vector_space>(ref_);
  }


  /*!
     An actor to push back an arbitrary Object into the appropriate std::vector
     The actor takes the element to push back and the vector as arguments
  */
  template<typename T>
  struct push_back_functor
  {
  public:
    push_back_functor(T& so_, std::vector<T>& soVec_) :
        so(so_), soVec(soVec_)
    {
    }

    void operator() (const char *first, const char *last) const
    {
        soVec.push_back(so);

    }

  private:
    T& so;
    std::vector<T>& soVec;
  };

   /*!
     An actor to push back a vector into the appropriate std::vector<std::vector<Type> >
     The actor takes the vector to push back and the vector as arguments
     The vector gets cleared after push_back
  */
  template<typename T>
  struct push_back_vector_functor
  {
  public:
    push_back_vector_functor(std::vector<T>& so_, std::vector<std::vector<T> >& soVec_) :
        so(so_), soVec(soVec_)
    {
    }

    void operator() (const char *first, const char *last) const
    {
        soVec.push_back(so);
        so.clear();
    }

  private:
    std::vector<T>& so;
    std::vector<std::vector<T> >& soVec;
  };


  /*!
      Actors setting fe_base::Node members
  */
  struct SetNodeX
  {
    template <class NodeT>
    void act( NodeT& node_, double const& value_) const
    {
        node_.SetX( value_ );
    }
  };
  template<class NodeT>
  inline boost::spirit::ref_value_actor<NodeT,SetNodeX>
  SetNodeX_a(NodeT& ref_)
  {
     return boost::spirit::ref_value_actor<NodeT,SetNodeX>(ref_);
  }

} // of namespace

#endif
