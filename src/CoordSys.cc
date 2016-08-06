//-----------------------------------------------------------------------------
// CoordSys.cc
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


#include "CoordSys.h"

extern const double Pi = 3.14159265359;

using namespace fe_base;

//---------------------------------------------
// functions member of the class 'CoordSys'
//---------------------------------------------

// Constructors
//! default constructor
CoordSys::CoordSys():CSType(not_determined), AType(not_set) {};


//! Specialized constructor
CoordSys::CoordSys( CoordSysType cstype_, AngleType atype_ )
  :CSType(cstype_), AType(atype_){}


//! Specialized constructor.
CoordSys::CoordSys( CoordSysType sysType_, AngleType atype_, DataType dataType_,
                    double x_, double y_, double z_ )
  :CSType(sysType_), AType(atype_){

  Set(dataType_, x_, y_, z_);
}

//!Specialized constructor
CoordSys::CoordSys( CoordSysType sysType_, AngleType atype_, DataType dataType_, Dense_Vector vec_)
  :CSType(sysType_), AType(atype_) {
  double x = vec_[0], y = vec_[1], z = vec_[2];
  Set(dataType_, x, y, z);
}

//! Specialized constructor.
CoordSys::CoordSys( CoordSysType sysType_, AngleType atype_, DataType dataType_,
                    double x0_, double y0_, double z0_,
                    double x_, double y_, double z_ )
  :CSType(sysType_), AType(atype_){
  Set(dataType_,x0_, y0_, z0_, x_, y_, z_);
}


// ////
// //// Sets and gets
// ////

//Get a single value from the map
double CoordSys::Get( string key ) {
  if (CheckForCoord( key))
    // if ( true )
    return DoubleMap[key];
  else return 0.0;
}

Dense_Vector CoordSys::Get( DataType dataType_ ){

  Dense_Vector res(3,0.0);

  switch (dataType_)
    {
    case Euler312 :
      {
        res[0] =  Get("Thxy");
        res[1] =  Get("Thyz");
        res[2] =  Get("Thzx");
        break;
      }
    case Euler313 :
      {
        res[0] =  Get("Phi");
        res[1] =  Get("Theta");
        res[2] =  Get("Psi");
        break;
      }
    case Rotvec :
      {
        res[0] =  Get("x_Rv");
        res[1] =  Get("y_Rv");
        res[2] =  Get("z_Rv");
        break;
      }
    case Base :
      {
        res[0] =  Get("x0");
        res[1] =  Get("y0");
        res[2] =  Get("z0");
        break;
      }
    default :
      {
        cerr << endl << "ERROR in CoordSys::Get( string ) "
             << endl << "The dataType " << dataType_ << " is not implemented or supported" << endl;
        exit(1);
      }
    }// end of switch
  return res;

}

string CoordSys::GetType(){

  string bla;

  if ( CSType == 0 ) bla = "not_determined";
  else if ( CSType == 1 ) bla = "cartesian";
  else if ( CSType == 2 ) bla = "cylindrical";

  return bla;

}

string CoordSys::GetAT(){

  string bla;

  if ( AType == 0 ) bla = "not_set";
  else if ( AType == 1 ) bla = "deg";
  else if ( AType == 2 ) bla = "rad";

  return bla;
}


Dense_Matrix CoordSys::GetRotMat(){

  Dense_Matrix rotmat(3,3);
  mtl::set_value(rotmat,0.0);

  if ( (IsKeySet("Thxy") && IsKeySet("Thyz") && IsKeySet("Thzx") ) ||
       (IsKeySet("Phi") && IsKeySet("Theta") && IsKeySet("Psi") ))
  {
    mtl::copy( Euler2Mat(), rotmat );
  }
  else if (RotvecKeys("x_Rv") && RotvecKeys("y_Rv") && RotvecKeys("z_Rv")){
    mtl::copy( Rotvec2Mat(), rotmat );
  }
  else {
    cerr << endl << "ERROR: In CoordSys::GetRotMat()"
         << endl << "ERROR: one or more parameters are missing"
         << endl << "ERROR: to calculate the rotational matrix" << endl;
    exit(1);
  }

  return rotmat;
}



Dense_Matrix CoordSys::GetLocalCS( Dense_Vector localBase_ ){

  if (CSType == cartesian){
    return GetRotMat();
  }
  else if ( CSType == cylindrical ){

    double norm2 = 0.0;

    Dense_Vector localX(3,0.0), localY(3,0.0), localZ(3,0.0);

    Dense_Matrix RotMat(3,3);
    mtl::set_value(RotMat,0.0);

    mtl::copy(GetRotMat(),RotMat);


//     cout << endl << "The point where the local coordsys is built" << endl;
//     print_vector(localBase_);
//     cout << endl << "The base point of the local coordsys" << endl;
//     print_vector(Get(Base));

    //get the localY as crossprod
    //--------------------------
    // calculate the connection globalBase, localBase
    mtl::copy(RotMat[2],localZ);
    mtl::add(localBase_,scaled(Get( Base ),-1.0),localX);

//     cout << endl << "connection from global to local base" << endl;
//     print_vector(localX);

    cross_prod_3d(localZ,localX,localY);
    norm2 = two_norm(localY);
    mtl::scale(localY,(1.0/norm2));

//     cout << endl << "localY" << endl;
//     print_vector(localY);


    //Check if the local base lies near the local coordsys
    //z-axis. if true give back the local coordsys
    //not sure if ansys does this as well...
    if ( norm2 < 0.1 )
      return RotMat;
    //else return the coordsys
    else {
      mtl::set_value(localX,0.0);
      cross_prod_3d(localY,localZ,localX);

      //write the coordsys
      mtl::copy(localX,RotMat[0]);
      mtl::copy(localY,RotMat[1]);


      return RotMat;
    }
  }

  else {
    cerr << endl << "ERROR: The fuction CoordSys::GetLocalCylCS(...) is not supported"
         << endl << "ERROR: or not yet implemented for the CoordSysTyp " << GetType() << endl;
    exit(1);
  }

}





void CoordSys::Set( string key_, double value_ ){
  if (CSType == not_determined){
    cerr << endl << "ERROR in CoordSys::Set( string, double ) "
         << endl << "The coordinate system type is not yet set"
         << endl << "Please set  a coordsystype " << endl;
    exit(1);
  }
  if (AType == not_set){
    cerr << endl << "ERROR in CoordSys::Set( string, double ) "
         << endl << "The unit for the angles is not yet set"
         << endl << "Please set these units" << endl;
    exit(1);
  }
  else if ( BaseKeys(key_) || Euler312Keys(key_) || Euler313Keys(key_) || RotvecKeys(key_) )
    DoubleMap[key_] = value_;
  else {
    cerr << endl << "ERROR in CoordSys::Set( string, double ) "
         << endl << "The value type " << key_ << " is not implemented or supported" << endl;
    exit(1);
  }
}


void CoordSys::Set( DataType dataType_,  double x_, double y_, double z_ ){
  if (CSType == not_determined){
    cerr << endl << "ERROR in CoordSys::Set( DataType, double, double, double ) "
         << endl << "The coordinate system type is not yet set"
         << endl << "Please set always a coordsystype " << endl;
    exit(1);
  }
  if (AType == not_set){
    cerr << endl << "ERROR in CoordSys::Set( string, double^3 ) "
         << endl << "The unit for the angles is not yet set"
         << endl << "Please set these units" << endl;
    exit(1);
  }
  switch (dataType_)
    {
    case Euler312 :
      {
        DoubleMap["Thxy"] = x_;
        DoubleMap["Thyz"] = y_;
        DoubleMap["Thzx"] = z_;
        break;
      }
    case Euler313 :
      {
        DoubleMap["Phi"] = x_;
        DoubleMap["Theta"] = y_;
        DoubleMap["Psi"] = z_;
        break;
      }
    case Rotvec :
      {
        DoubleMap["x_Rv"] = x_;
        DoubleMap["y_Rv"] = y_;
        DoubleMap["z_Rv"] = z_;
        break;
      }
    case Base :
      {
        DoubleMap["x0"] = x_;
        DoubleMap["y0"] = y_;
        DoubleMap["z0"] = z_;
        break;
      }
    default :
      {
        cerr << endl << "ERROR in CoordSys::Set( DataType, double, double, double) "
             << endl << "The DataType " << dataType_ << " is not implemented or supported" << endl;
        exit(1);
      }
    }//end of switch

}


void CoordSys::Set( DataType dataType_, double x0_, double y0_, double z0_,
                    double x_, double y_, double z_ ){
  if (CSType == not_determined){
    cerr << endl << "ERROR in CoordSys::Set( DataType, double, double, double ) "
         << endl << "The coordinate system type is not yet set"
         << endl << "Please set always a coordsystype " << endl;
    exit(1);
  }
  if (AType == not_set){
    cerr << endl << "ERROR in CoordSys::Set( string, double^6) "
         << endl << "The unit for the angles is not yet set"
         << endl << "Please set these units" << endl;
    exit(1);
  }

  switch (dataType_)
    {
    case Euler312 :
      {
        DoubleMap["x0"] = x0_;
        DoubleMap["y0"] = y0_;
        DoubleMap["z0"] = z0_;
        DoubleMap["Thxy"] = x_;
        DoubleMap["Thyz"] = y_;
        DoubleMap["Thzx"] = z_;
        break;
      }
    case Euler313 :
      {
        DoubleMap["x0"] = x0_;
        DoubleMap["y0"] = y0_;
        DoubleMap["z0"] = z0_;
        DoubleMap["Phi"] = x_;
        DoubleMap["Theta"] = y_;
        DoubleMap["Psi"] = z_;
        break;
      }
    case Rotvec :
      {
        DoubleMap["x0"] = x0_;
        DoubleMap["y0"] = y0_;
        DoubleMap["z0"] = z0_;
        DoubleMap["x_Rv"] = x_;
        DoubleMap["y_Rv"] = y_;
        DoubleMap["z_Rv"] = z_;
        break;
      }
    case Base :
      {
        cerr << endl << "ERROR in CoordSys::Set( DataType, double, double, double...) "
             << endl << "The DataType " << dataType_ << " is set immediate with this "
             << endl << "function.... use the 3 double version to set the Base independently " << endl;
        exit(1);
      }
    default :
      {
        cerr << endl << "ERROR in CoordSys::Set( DataType, double, double, double) "
             << endl << "The DataType " << dataType_ << " is not implemented or supported" << endl;
        exit(1);
      }
    }//end of switch

}



void CoordSys::SetType( CoordSysType  type_ ){
    CSType = type_;
}

void CoordSys::SetAngleType( AngleType type_ ){
  AType = type_;
}

// ////
// //// Operator overloading
// ////

// /*! Overload "equal" and "not equal" operator
//   in order to apply find functions of STL-vector . */
bool CoordSys::operator==(const CoordSys& a_) const{
  return (a_.CSType == CSType && a_.DoubleMap == DoubleMap);
}

bool CoordSys::operator!=(const CoordSys& a_) const{
    return ( a_.CSType != CSType || a_.DoubleMap != DoubleMap );
}


// ////
// //// Other
// ////

bool CoordSys::CheckForCoord(string key){
  map<string, double>::iterator MapIt;
  MapIt = DoubleMap.begin();
  MapIt = DoubleMap.find(key);
  if (MapIt == DoubleMap.end()){
    cerr << "ERROR: In CoordSys::CheckForCoord(map, string)" << endl;
    cerr << "ERROR: A CoordSys map has been asked for a non-existing entry" << endl;
    cerr << "ERROR: The entry falsely looked for is named: " << key << endl << endl;
    return 0;
  }
  else return 1;
}

bool CoordSys::BaseKeys(string key){
  return ( key == "x0" || key == "y0" || key == "z0");
}

bool CoordSys::Euler312Keys(string key){
  return ( key == "Thxy" || key == "Thyz" || key == "Thzx");
}

bool CoordSys::Euler313Keys(string key){
  return ( key == "Phi" || key == "Theta" || key == "Psi" );
}

bool CoordSys::RotvecKeys(string key){
  return ( key == "x_Rv" || key == "y_Rv" || key == "z_Rv");
}

bool CoordSys::IsKeySet(std::string key)
{
  std::map<std::string, double>::iterator findIter;
  findIter = DoubleMap.find(key);
  bool found = false;
  if ( findIter != DoubleMap.end() )
    found = true;

  return found;
}


/*! Function that gives back 3x3 rotational matrix from the three Euler-Coordinate angles.
The angles are given as defined in ansys... The rotations take place in the order given.
The angles are given in degrees or radiants. The stored angular values will not be changed.
 */
Dense_Matrix CoordSys::Euler2Mat(){

  if ( IsKeySet("Thxy") &&  IsKeySet("Thyz") &&  IsKeySet("Thzx") )
  {
    Dense_Matrix rot1(3, 3), rot2(3,3), rot3(3,3);
    mtl::set_value(rot1,0.0);
    mtl::set_value(rot2,0.0);
    mtl::set_value(rot3,0.0);

    if ( AType == deg )
    {
      double tmpX, tmpY, tmpZ;

      tmpX = Deg2Rad(DoubleMap["Thxy"]);
      tmpY = Deg2Rad(DoubleMap["Thyz"]);
      tmpZ = Deg2Rad(DoubleMap["Thzx"]);

      copy(Vec2Mat(0., 0., tmpX),rot1);
      // rot2 = Vec2Mat(Thyz, 0, 0) * rot1;
      mult(Vec2Mat(tmpY,0.,0.),rot1,rot2);
      // rot3 = Vec2Mat(0, Thzx, 0) * rot2;
      mult(Vec2Mat(0.,tmpZ,0.),rot2,rot3);
    }
    else if ( AType == rad )
    {

      copy(Vec2Mat(0., 0., DoubleMap["Thxy"]),rot1);
      // rot2 = Vec2Mat(Thyz, 0, 0) * rot1;
      mult(Vec2Mat(DoubleMap["Thyz"],0.,0.),rot1,rot2);
      // rot3 = Vec2Mat(0, Thzx, 0) * rot2;
      mult(Vec2Mat(0.,DoubleMap["Thzx"],0.),rot2,rot3);
    }
    else
    {
      cerr << endl << "ERROR in CoordSys::Euler2Mat() "
          << endl << "The unit for the angles is not yet set"
          << endl << "Please set these units" << endl;
    }
    return rot3;
  } // end of the  calculation of euler312 angles
  else if (  IsKeySet("Phi") &&  IsKeySet("Theta") &&  IsKeySet("Psi") )
  {
    Dense_Matrix rot1(3, 3), rot2(3,3), rot3(3,3);
    mtl::set_value(rot1,0.0);
    mtl::set_value(rot2,0.0);
    mtl::set_value(rot3,0.0);

    if ( AType == deg )
    {
      double tmpPhi, tmpTheta, tmpPsi;

      tmpPhi = Deg2Rad(DoubleMap["Phi"]);
      tmpTheta = Deg2Rad(DoubleMap["Theta"]);
      tmpPsi = Deg2Rad(DoubleMap["Psi"]);

      //rot1 = Vec2Mat(0,0,Phi)
      copy(Vec2Mat(0., 0., tmpPhi),rot1);
      // rot2 = Vec2Mat(Theta, 0, 0) * rot1;
      mult(Vec2Mat(tmpTheta,0.,0.),rot1,rot2);
      // rot3 = Vec2Mat(0, 0, Psi) * rot2;
      mult(Vec2Mat(0.,0.,tmpPsi),rot2,rot3);
      //rot3 = Psi*Theta*Phi

    }
    else if ( AType == rad )
    {

      copy(Vec2Mat(0., 0., DoubleMap["Phi"]),rot1);
      // rot2 = Vec2Mat(Thyz, 0, 0) * rot1;
      mult(Vec2Mat(DoubleMap["Theta"],0.,0.),rot1,rot2);
      // rot3 = Vec2Mat(0, Thzx, 0) * rot2;
      mult(Vec2Mat(0.,0.,DoubleMap["Psi"]),rot2,rot3);
      //rot3 = Psi*Theta*Phi
    }
    else
    {
      cerr << endl << "ERROR in CoordSys::Euler2Mat() "
          << endl << "The unit for the angles is not yet set"
          << endl << "Please set these units" << endl;
    }
    return rot3;

  }
  else
  {
    cerr << endl << "ERROR : In CoordSys::Euler2Mat "
         << endl << "The calling CoordSys has one or more of the Euler-angles not set" << endl;
    exit(1);
  }
}

/*! Function that gives back 3x3 rotational matrix. The three arguments provided to the
  function building a vector of whom the direction gives the rotational axis and the
  length provides the rotation angle in radiants. */
Dense_Matrix CoordSys::Rotvec2Mat(){

  if (RotvecKeys("x_Rv") && RotvecKeys("y_Rv") && RotvecKeys("z_Rv") ){
    return Vec2Mat(DoubleMap["x_Rv"], DoubleMap["y_Rv"], DoubleMap["z_Rv"]);
  }
  else {
    cerr << endl << "ERROR : In CoordSys::Rotvec2Mat "
         << endl << "The calling CoordSys has one or more of the Rotvec coordinates  not set" << endl;
    exit(1);
  }
}

//This function transformes a vector given by 3 doubles
//to a rotational matrix. The vector direction gives the
//axis of the rotation and its length the angle of rotation
//in radiants
Dense_Matrix CoordSys::Vec2Mat(double x_, double y_, double z_){

    Dense_Matrix rot(3, 3);
    mtl::set_value(rot,0.0);
    if ( x_ == 0.0 && y_ == 0.0 && z_ == 0.0 ){
      mtl::set_diagonal(rot,1.0);
    }
    else {
      double dummy1 = sqrt(x_*x_ + y_*y_ + z_*z_);
      double dummy2 = sin(dummy1)/dummy1;
      double dummy3 = (1-cos(dummy1))/dummy1/dummy1;

      rot(0, 0) = x_*x_*dummy3 + cos(dummy1);
      rot(1, 1) = y_*y_*dummy3 + cos(dummy1);
      rot(2, 2) = z_*z_*dummy3 + cos(dummy1);

      rot(0, 1) = z_*dummy2 + x_*y_*dummy3;
      rot(1, 0) = -z_*dummy2 + x_*y_*dummy3;

      rot(0, 2) = -y_*dummy2 + x_*z_*dummy3;
      rot(2, 0) = y_*dummy2 + x_*z_*dummy3;

      rot(1, 2) = x_*dummy2 + y_*z_*dummy3;
      rot(2, 1) = -x_*dummy2 + y_*z_*dummy3;
    }
    return rot;
}


double CoordSys::Deg2Rad( double a_){
  double PI = acos(-1.0);
  return (a_*PI)/180;
}

double CoordSys::Rad2Deg( double a_){
  double dummy = 1.0/(Deg2Rad(a_));
  return dummy;
}


//Gives back a vector parallel to the rotation axis with a length
//equal to the angle of rotation in radiants. The parameter D_ is
//changed by this function....
Dense_Vector CoordSys::Mat2Rotvec(Dense_Matrix D_){
  //3D check
  int rows = D_.nrows(), cols = D_.ncols();

  if ( rows != 3 || cols != 3 ){
    cerr << endl << "******************************************"
         << endl << "*****  ERROR in Mat2Rotvec(matrix) *******"
         << endl << "****The matrix dimension is not (3,3)*****"
         << endl << "******************************************" << endl;
  }

  //orthogonality-check for the matrix D_
  orthogonalityCheck(D_, 0);


  Dense_Vector globalX(3,0.0), localX(3,0.0), eigenvector(3,0.0), result(3,0.0);

  //The axis of rotation is given as the only invariant vector which is an eigenvector
  //since the matrix is not symmetric it is not so easy to compute it.
  //That is why we implemented a very special form of eigenvector evaluation specially
  //tailored to this problem. For a rotation matrix, the corresponding diskriminant D  is
  // D > 0. Means, the matrix has two conjugte complex and one real eigenvalue
  //The real eigenvalue is the one we have to find the eigenvector for!

  //The determinante of D_ is built as
  //det(D_) = x^3+ax^2+bx+c = 0
  //maple gave the following polynomial
  /*
  -x^3
  + (t11 + t33 + t22) x^2
  + (t31 t13 - t11 t33 - t22 t33 + t21 t12 - t11 t22 + t23 t32) x
  + t31 t12 t23 - t11 t23 t32 - t31 t13 t22 - t21 t12 t33 + t11 t22 t33 + t21 t13 t32 = 0
  */
  double a = (-1.0)*(D_(0,0) + D_(1,1) + D_(2,2)),
    b = (-1.0)*(D_(2,0)*D_(0,2)-D_(0,0)*D_(2,2)-D_(1,1)*D_(2,2)+D_(1,0)*D_(0,1)-D_(0,0)*D_(1,1)+D_(1,2)*D_(2,1)),
    c = (-1.0)*(D_(2,0)*D_(0,1)*D_(1,2)-D_(0,0)*D_(1,2)*D_(2,1)-D_(2,0)*D_(0,2)*D_(1,1)-D_(1,0)*D_(0,1)*D_(2,2)+D_(0,0)*D_(1,1)*D_(2,2)+D_(1,0)*D_(0,2)*D_(2,1));

//   cout << endl << "The incredible coefficients" << endl
//        << a << endl << b << endl << c << endl ;

  //The diskriminante is computed according to Papula (Mathematische Formelsammlung)
  double p = (3.0*b - a*a)/3.0,
    q = (2.0*a*a*a/27.0)-(a*b/3.0)+c;

  double diskr = (pow((p/3.0),3)) + (pow((q/2.0),2));

  //cout << endl << " The coefficients are " << endl << p << endl << q << endl;
  // cout << endl << " The diskriminante is : " << diskr << endl;
  double realEigenvalue = 0.0;

  if ( diskr > 0 ){ // exactly one real solution for the polynomial
    double onethird = 1.0/3.0, sqrtdiskr = sqrt(diskr);

    //array of the two variables u, v
    double u[2];

    u[0] =  -q/2.0+sqrtdiskr;
    u[1] = -q/2.0-sqrtdiskr;

    for (int i = 0 ; i < 2 ; ++i ){
      if ( u[i] <= 0 ){
        u[i] = pow(abs(u[i]),onethird);
        u[i] = -u[i];
      }
      else {
        u[i] = pow(u[i],onethird);
      }
    }
    realEigenvalue = u[0]+u[1]-a/3;
  }
  else {
    cerr << endl << " Error in Mat2Rotvec: The diskriminant of the matrix is not > 0 " << endl;
  }

  //cout << endl << " The real eigen value is : " << realEigenvalue << endl;

  //Find the eigenvector to the real eigenvalue...
  //the matrix D_-realEigenvalue*Identity is singular
  //We use svdecomposition to solve the singular system of equations in its nullspace

  Dense_Matrix Identity(3,3);
  mtl::set_value(Identity, 0.0); mtl::set_diagonal(Identity,-realEigenvalue);

  Dense_Matrix U(3,3);
  mtl::set_value(U,0.0);
  copy(D_,U);
  add(Identity,U); //build the singular matrix

  //Set some Matrices and vectors for svdecomposition
  Dense_Matrix V(3,3);
  mtl::set_value(V,0.0);
  Dense_Vector w(3,0.0);


  svdcmp(U,w,V);

//  cout << endl << "U " << endl;
//  print_all_matrix(U);
//  cout << endl << " w " << endl;
//  print_vector(w);
//  cout << endl << " V " << endl;
//  print_all_matrix(V);


  //  look for the smallest value of w and check, if it is zero
  double thesmallestvalue = nr_min(nr_min(abs(w[0]),abs(w[1])),abs(w[2]));
  int i = 0;
  while ( abs(w[i]) != thesmallestvalue ) ++i;

  if (abs(thesmallestvalue) > 1e-6){
    cout << endl << "ERROR in Rotvec2Mat: The smallest value of the V vector is > 1e-6 (it should be zero: check the singularity of the Matrix)"
         << endl << "The error only occured for the coordsys (1,0,1). For (1,0,1.0000000000001) it disapears...There may be other CoordSys where"
         << endl << "it appears. The problem could be within svdcmp()...... Check Numerical Recipes for more information on the function " << endl ;
    exit(0);
 }

  for (int k = 0 ; k < 3 ; ++k) eigenvector[k] = V(k,i);
  //normalizing the eigenvector and multypling it with the sign of its z coordinate
  //This mirrors the vector to get only vectors in the first four octants
  double norm2 = 1.0 / two_norm(eigenvector);
  scale(eigenvector,norm2);

  //The angle of rotation is given in the normal plane to the eigenvector
  globalX[0] = 1.0;
  //We use the scalar product to determine the footpoint
  double scalarprod = globalX[0]*eigenvector[0];

  //check if the axis is near to the global x-axis
  //if it is we use y-axis for the rotation
  if (1.0-scalarprod < 0.01){
    globalX[0] = 0.0;
    globalX[1] = 1.0;
    scalarprod = globalX[1]*eigenvector[1];
    copy(D_[1],localX);
  }
  else{
    copy(D_[0],localX);
  }

  //  cout << "scalr " << scalarprod << endl;

  Dense_Vector footpoint(3,0.0);
  copy(scaled(eigenvector,-1.0*scalarprod),footpoint);

  Dense_Vector dummyVec1(3,0.0), dummyVec2(3,0.0);


  add(footpoint,globalX,dummyVec1);
  add(footpoint,localX,dummyVec2);

  scale(dummyVec1,two_norm(dummyVec1));
  scale(dummyVec2,two_norm(dummyVec2));


  //since both vectors are normalized this gives the sin of their angle
  cross_prod_3d(dummyVec1, dummyVec2, result);

  // The angle has to be unique, so we need the atan2(x,y)
  double cos_result, sin_result;
  sin_result = two_norm(result);
  cos_result = scalar_prod_3d(dummyVec1, dummyVec2);

  norm2 = atan2(sin_result,cos_result);


  // cout << endl << "The vector length is : " << norm2 << endl;

  //in order to get the right sense of rotation, one has to reach
  //dummyVec2 from dummyVec1  turning right....
  if ( scalar_prod_3d(result,eigenvector) < 0.0 ){
     norm2 *= -1.0;
   }
  scale(eigenvector,norm2);

  return eigenvector;
}// of fe_base::Mat2Rotvec



double CoordSys::nr_pythag(const double a, const double b){
  double absa, absb;

  absa = fabs(a);
  absb = fabs(b);

  if ( absa > absb ) return absa*sqrt(1.0+nr_sqr(absb/absa));
  else return ( absb == 0.0 ? 0.0 : absb*sqrt(1.0+nr_sqr(absa/absb)));
}




void CoordSys::svdcmp(Dense_Matrix a, Dense_Vector w, Dense_Matrix v)
{
  bool flag;
  int i=0, its=0, j=0, jj=0, k=0, l=0, nm=0;
  float_type anorm,c,f,g,h,s,scale,x,y,z;
  int m=a.nrows();
  int n=a.ncols();
  Dense_Vector rvl(n,0.0); //**
  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rvl[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a(k,i));
      if (scale != 0.0) {
        for (k=i;k<m;k++) {
          a(k,i) /= scale;
          s += a(k,i)*a(k,i);
        }
        f=a(i,i);
        g = -nr_sign(sqrt(s),f);
        h=f*g-s;
        a(i,i)=f-g;
        for (j=l-1;j<n;j++) {
          for (s=0.0,k=i;k<m;k++) s += a(k,i)*a(k,j);
          f=s/h;
          for (k=i;k<m;k++) a(k,j) += f*a(k,i);
        }
        for (k=i;k<m;k++) a(k,i) *= scale;
      }
    }
    w[i]=scale*g;

    g=s=scale=0.0;
    if (i+1 <= m && i != n) {
      for (k=l-1;k<n;k++) scale += fabs(a(i,k));
      if (scale != 0.0) {
        for (k=l-1;k<n;k++) {
          a(i,k) /= scale;
          s += a(i,k)*a(i,k);
        }
        f=a(i,l-1);
        g = -nr_sign(sqrt(s),f );
        h=f*g-s;
        a(i,l-1)=f-g;
        for (k=l-1 ;k<n;k++) rvl[k]=a(i,k)/h;
        for (j=l-1;j<m;j++) {
          for (s=0.0,k=l-1;k<n;k++) s += a(j,k)*a(i,k);
          for (k=l-1;k<n;k++) a(j,k) += s*rvl[k];
        }
        for (k=l-1;k<n;k++) a(i,k) *= scale;
      }
    }
    anorm=nr_max(anorm,(fabs(w[i])+fabs(rvl[i])));
  }
  for (i=n-1;i>=0;i--) {        //Accumulation of right-hand transformations.
    if (i < n-1) {
      if (g != 0.0) {
        for (j=l;j<n;j++)       //Double division to avoid possible underflow.
          v(j,i)=(a(i,j)/a(i,l))/g;
        for (j=l;j<n;j++) {
          for (s=0.0,k=l;k<n;k++) s += a(i,k)*v(k,j);
          for (k=l;k<n;k++) v(k,j) += s*v(k,i);
        }
      }
      for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
    }
    v(i,i) =1.0;
    g=rvl[i];
    l=i;
  }
  for (i=nr_min(m,n)-1;i>=0;i--) {      //Accumulation of left-hand transformations.
    l=i+1;
    g=w[i];
    for (j=l ; j<n ; j++) a(i,j) =0.0 ;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
        for (s=0.0,k=l;k<m;k++) s += a(k,i)*a(k,j);
        f=(s/a(i,i))*g;
        for (k=i;k<m;k++) a(k,j) += f*a(k,i);
      }
      for (j=i;j<m;j++) a(j,i) *= g;
    } else for (j=i;j<m;j++) a(j,i)=0.0;
    ++a(i,i);
  }
  for (k=n-1;k>=0;k--) {        //Diagonalizationofthebidiagonalform: Loopover
    for (its=0; its<30; its++) {        //singular values, and over allowed iterations.
      flag=true;
      for (l=k;l>=0;l--) {      //Test for splitting.
        nm=l-1; //Note that rvl [0] is always zero.
        if (fabs(rvl[l])+anorm == anorm) {
          flag=false;
          break;
        }
        if (fabs(w[nm])+anorm == anorm) break;
      }

      if (flag) {
        c=0.0; //Cancellation of rvl[l], if 1 > 0.
        s=1.0;
        for (i=l-1;i<k+1;i++) {
          f=s*rvl[i];
          rvl[i]=c*rvl[i];
          if (fabs(f)+anorm == anorm) break;
          g=w[i];
          h=nr_pythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s = -f*h;
          for (j=0;j<m;j++) {
            y=a(j,nm);
            z=a(j,i);
            a(j,nm)=y*c+z*s;
            a(j,i)=z*c-y*s;
          }
        }
      }
      z=w[k];
      if (l == k) {     //Convergence.
        if (z < 0.0) {  //Singular value is made nonnegative.
          w[k] = -z;
          for (j=0;j<n;j++) v(j,k) = -v(j,k);
        }
        break;
      }
      //if (its == 29) cerr << endl << "no convergence in 30 svdcmp iterations" << endl;
      x=w[l];   //Shift from bottom 2-by-2 minor.
      nm=k-1;
      y=w[nm];
      g=rvl[nm];
      h=rvl[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=nr_pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+nr_sign(g,f)))-h))/x;
      c=s=1.0;//        Next QR transformation:
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rvl[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=nr_pythag(f,h);
        rvl[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y *= c;
        for (jj=0;jj<n;jj++) {
          x=v(jj,j);
          z=v(jj,i);
          v(jj,j)=x*c+z*s;
          v(jj,i)=z*c-x*s;
        }
        z=nr_pythag(f,h);
        w[j]=z; //Rotation can be arbitrary if z = 0.
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=0;jj<m;jj++) {
          y=a(jj,j);
          z=a(jj,i);
          a(jj,j)=y*c+z*s;
          a(jj,i)=z*c-y*s;

        }
      }
      rvl[l]=0.0;
      rvl[k]=f;
      w[k]=x;
    }
  }
}






//overloaded '<<'-operator to easily print single node-coordinate systems
ostream& fe_base::operator<< ( ostream& stream_, CoordSys& a_ )
{
  if ( &a_ != NULL ){
    stream_.precision(8);

    stream_ << "  ";
    stream_.width( 15 );
    stream_ << a_.GetType() << " CoordSys in " <<  a_.GetAT() << endl;

    map<string, double>::const_iterator pos;
    for (pos = a_.DoubleMap.begin() ; pos != a_.DoubleMap.end() ; ++pos){
      stream_ << "  ";
      stream_.width( 15 - pos -> first.size() );
      stream_ <<  pos -> first << " = " << pos -> second ;
    }
    stream_.unsetf(ios::left);
  }
  else{
    stream_ << endl << " ERROR: In operator<< (ostream, CoordSys& )"
            << endl << " ERROR: The calling pointer to CoordSys is set to NULL " << endl;
  }
  return stream_;
}



