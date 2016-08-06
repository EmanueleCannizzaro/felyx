////////////////////////////////////////////////////////////////////////////////
// The Loki Library
// Copyright (c) 2001 by Andrei Alexandrescu
// This code accompanies the book:
// Alexandrescu, Andrei. "Modern C++ Design: Generic Programming and Design 
//     Patterns Applied". Copyright (c) 2001. Addison-Wesley.
// Permission to use, copy, modify, distribute and sell this software for any 
//     purpose is hereby granted without fee, provided that the above copyright 
//     notice appear in all copies and that both that copyright notice and this 
//     permission notice appear in supporting documentation.
// The author or Addison-Welsey Longman make no representations about the 
//     suitability of this software for any purpose. It is provided "as is" 
//     without express or implied warranty.
////////////////////////////////////////////////////////////////////////////////

// Last update: June 20, 2001

// Modified by okoenig in order to make it run on g++
// - replaced AssocVector by map
// - simplified Exception handling cause of compilation errors on g++ 3.01.
// - deleted CloneFactory cause it's not needed right now!

#ifndef FACTORY_INC_
#define FACTORY_INC_

#include <map>
#include <iostream>

namespace Loki
{

////////////////////////////////////////////////////////////////////////////////
// class template DefaultFactoryError
// Manages the "Unknown Type" error in an object factory
////////////////////////////////////////////////////////////////////////////////

// The following function has been modified by okoenig, cause of compilation
// problems on gcc 3.0.1
  
    template <typename IdentifierType, class AbstractProduct>
    struct DefaultFactoryError
    {
        static AbstractProduct* OnUnknownType(IdentifierType)
        {
//            throw Exception();
	  cerr << "ERROR in DefaultFactoryError: Unknown Type" << endl;
	  exit(1);
        }
    };

////////////////////////////////////////////////////////////////////////////////
// class template Factory
// Implements a generic object factory
////////////////////////////////////////////////////////////////////////////////

    template
    <
        class AbstractProduct, 
        typename IdentifierType,
        typename ProductCreator = AbstractProduct* (*)(),
        template<typename, class>
            class FactoryErrorPolicy = DefaultFactoryError
    >
    class Factory 
        : public FactoryErrorPolicy<IdentifierType, AbstractProduct>
    {
    public:
        bool Register(const IdentifierType& id, ProductCreator creator)
        {
            return associations_.insert(
                IdToProductMap::value_type(id, creator)).second;
        }
        
        bool Unregister(const IdentifierType& id)
        {
            return associations_.erase(id) == 1;
        }
        
        AbstractProduct* CreateObject(const IdentifierType& id)
        {
            typename IdToProductMap::const_iterator i = associations_.find(id);
            if (i != associations_.end())
            {
                return (i->second)();
            }
            return OnUnknownType(id);
        }
        
    private:
      typedef map<IdentifierType, ProductCreator> IdToProductMap;
      
      //typedef AssocVector<IdentifierType, ProductCreator> IdToProductMap;
        IdToProductMap associations_;
    };
}

  
////////////////////////////////////////////////////////////////////////////////
// Change log:
// June 20, 2001: ported by Nick Thurn to gcc 2.95.3. Kudos, Nick!!!
////////////////////////////////////////////////////////////////////////////////

#endif // FACTORY_INC_
