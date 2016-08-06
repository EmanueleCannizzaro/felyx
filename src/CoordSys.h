//-----------------------------------------------------------------------------
// CoordSys.h
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


#ifndef CoordSys_h
#define CoordSys_h CoordSys_h

#include <cmath>
#include <iostream>
#include <map>

#include "mtl/mtl_felyx_utils.h"

using namespace std;
using namespace mtl;



namespace fe_base{      // introduces the namespace fe_base

  /*! Class for all coordinate systems used in FELyX
    Implemented are only cartesian and cylindrical coordinate systems.
    The location of the coordinate system is not always stored (e.g
    if a coordinate system defines a rotation). The coordinate system
    or rotation is stored in two sets of three doubles. One set for
    the base point and one set for the rotation. There are functions
    to transform a rotational vector and euler angles to rotational
    matrices.
    Euler312 represents the ansys rotation sequence about local axis 3-1-2 (about Thxy, Thyz, and Thzx)
    Euler313 represents the rotation sequence about local 3-1-3 axis (about Phi, Theta, and Psi)
  */

  class CoordSys{

  public:
    enum AngleType { not_set, deg, rad };
    enum CoordSysType { not_determined, cartesian, cylindrical };
    enum DataType { Euler312, Rotvec, Base, Euler313 };

    // Constructors
    // ------------
    CoordSys();                         //default constructor
    CoordSys( CoordSysType, AngleType );                 //sets the CoordSys type
    CoordSys( CoordSysType, AngleType, DataType, double, double, double );
    CoordSys( CoordSysType, AngleType, DataType, double, double, double,
              double, double, double );
    CoordSys( CoordSysType, AngleType, DataType, Dense_Vector );

    // Sets and gets
    // -------------
    Dense_Vector Get( DataType );
    double Get( string );
    string GetType();
    string GetAT();
    AngleType GetAngleType() { return AType; }
    CoordSysType GetCSType() { return CSType; }

    //!checks for the coordinates stored and returns the rotational matrix
    Dense_Matrix GetRotMat();
    /*!Gives back the local coordsys for a certain point. the point must
      be given in global coordinates
    */
    Dense_Matrix GetLocalCS( Dense_Vector );


    void Set( string, double ); //
    void Set( DataType, double, double, double ); //for rotations only or base only
    void Set( DataType, double, double, double,
              double, double, double ); //for base and rotation
    void SetType( CoordSysType );
    void SetAngleType( AngleType );

    // Operator overloading
    // --------------------
    bool operator== (const CoordSys&) const;    // overload "equal" and "not equal" operator
    bool operator!= (const CoordSys&) const;    // to apply find functions of STL-vector


    // Other functions
    // ---------------
    bool CheckForCoord( string );
    bool BaseKeys( string );
    bool Euler312Keys( string );
    bool Euler313Keys( string );
    bool RotvecKeys( string );
    bool IsKeySet( string );

    // overload <<-operator
    friend ostream& operator<<(ostream&, CoordSys&);

    // Transformations
    // -----
    Dense_Matrix Euler2Mat();           //transforms Euler-angels to a rotational matrix
    Dense_Matrix Rotvec2Mat();          //transforms the Rotvec to  a rotational matrix

  private:
    // Data members
    // ------------

    CoordSysType CSType; //cartesian, cylindrical
    AngleType AType; // deg, rad
    map<string, double> DoubleMap;

  public:
  static Dense_Matrix Vec2Mat(double, double, double);

  //produces the rotational vector from the matrix of direction cosines
  static Dense_Vector Mat2Rotvec( Dense_Matrix );

  //!Converts degrees to radiants
  static double Deg2Rad(double);

  static double Rad2Deg(double);

  //a function to solve singular simultaneaous equations
  //check Numerical Recipes for more information
  static void svdcmp(Dense_Matrix, Dense_Vector, Dense_Matrix);

  //Some functions used in svdcmp please check Numerical Recipes for more
  //information
  static double nr_pythag(const double, const double );
 };
  template<class T>
    T nr_sign(const T& a, const T& b) {
    return b >= 0.0 ? (a >= 0.0 ? a : -a) : (a >= 0.0 ? -a : a);
  }
  template<class T>
  T nr_max(const T& a, const T& b){
    return b > a ? (b) : (a);
  }
  template<class T>
  T nr_min(const T& a, const T& b){
    return b < a ? (b) : (a);
  }

  template<class T>
    T nr_sqr(const T a){
    return a*a;
  }

  // overload <<-operator
  ostream& operator<<(ostream&, const CoordSys&);

}

#endif


