//-----------------------------------------------------------------------------
// FELyX_Exception.h
//
// begin     : March 2004
// copyright : (c) 2004 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
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
#ifndef FELyX_EXCEPTION_H
#define FELyX_EXCEPTION_H

#include <exception>
#include <string>

namespace fe_base {

/**
  The assertion class defines the type of exception that FELyX functions
  will throw.
  Basically std::logic_error and std::runtime_error are used, 
  where std::logic_error can be disabled, if the macro FELYX\_EXCEPTIONS is not defined.
  
  Note that you need not put a {\tt try} and {\tt catch} around
  each FELyX function call. A better style is to make your whole
  function (perhaps with many calls to FElLyX) a {\tt try} clause.
  This way the exception handling code does not interfere with
  the readability of your algorithm. In addition, this will
  result in your having to write fewer {\tt try} and {\tt catch}
  clauses.
*/
class logic_error : public std::exception {
public:
  logic_error(const char* error){
    desc_ = "FELyX logic error: ";
    desc_ += error;
  }
  
  logic_error(const char* assertion, const char* function) {
    desc_ = "FELyX logic assertion: ";
    desc_ += assertion;
    desc_ += " failed in ";
    desc_ += function;
  }

  virtual ~logic_error() throw() { }

  virtual const char* what() const throw() { return desc_.c_str(); }

protected:
  std::string desc_;

};

class runtime_error : public std::exception {
public:
  runtime_error(const char* error){
    desc_ = "FELyX runtime error: ";
    desc_ += error;
  }
  
  runtime_error(const char* assertion, const char* function) {
    desc_ = "FELyX runtime assertion: ";
    desc_ += assertion;
    desc_ += " failed in ";
    desc_ += function;
  }

  virtual ~runtime_error() throw() { }

  virtual const char* what() const throw() { return desc_.c_str(); }

protected:
  std::string desc_;

};

} /* namespace fe_base */

// Runtime exceptions are always thrown!
#define FELYX_RUNTIME_THROW(X)           throw( fe_base::runtime_error(X) )
#define FELYX_RUNTIME_ASSERT(X,Y)        if (!(X)) { throw fe_base::runtime_error(#X,Y); }


// Logic exceptions can be disabled with macro FELYX/_EXCEPTIONS
#if FELYX_EXCEPTIONS
#define FELYX_LOGIC_THROW(X)           throw( fe_base::logic_error(X) )
#define FELYX_LOGIC_ASSERT(X,Y)        if (!(X)) { throw fe_base::logic_error(#X,Y); }
#else
#define FELYX_LOGIC_THROW(X)           /* nothing */
#define FELYX_LOGIC_ASSERT(X,Y)        /* nothing */
#endif


#endif
