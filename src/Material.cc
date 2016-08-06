
//-----------------------------------------------------------------------------
// Material.cc
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

#include "Material.h"

using namespace fe_base;

/////////////////////////////////////////
/// functions member of class 'Material'
/////////////////////////////////////////
void Material::Set( string mystring_, double value_ ){
  cerr << endl << "ERROR: In Material::Set(string, double) the material calling"
       << endl << "ERROR: has no Set function defined...." << endl;
  exit(1);
}

double Material::Get( string mystring_ ) const {
  cerr << endl << "ERROR: In Material::Get(string) the material calling"
       << endl << "ERROR: has no Get function defined...." << endl;
  exit(1);
  return 0.0;
}


// void Material::LayeredMatError() const  {
//   cerr << endl << "ERROR: In Material the material calling this function"
//        << endl << "ERROR: has no layers defined. It is a " << ClassName() << endl;
//   exit(1);
// }



///////////////////////////////////////////////////
/// functions member of class 'OrthotropicMaterial'
///////////////////////////////////////////////////

OrthotropicMaterial::OrthotropicMaterial(){

  Set(0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0,
        0.0, 0.0, 0.0,
        0.0);

  PlaneStress12_flag = not_inited;
  PlaneStrain12_flag = not_inited;
  RotationSymmetric3_flag = not_inited;
  ThreeDSolid_flag = not_inited;
  ThreeDShell_flag = not_inited;
}

OrthotropicMaterial::OrthotropicMaterial(double E1, double E2, double E3,
                                        double G12, double G31, double G23,
                                        double nu12, double nu13, double nu23,
                                        double rho,
                                        double K1, double K2, double K3,
                                        double visc)
{

  Set(E1, E2, E3,
        G12, G31, G23,
        nu12, nu13, nu23,
        rho,
        K1, K2, K3,
        visc);

  PlaneStress12_flag = not_inited;
  PlaneStrain12_flag = not_inited;
  RotationSymmetric3_flag = not_inited;
  ThreeDSolid_flag = not_inited;
  ThreeDShell_flag = not_inited;
}

void OrthotropicMaterial::Set(double E1, double E2, double E3,
                                double G12, double G31, double G23,
                                double nu12, double nu13, double nu23,
                                double rho,
                                double K1, double K2, double K3,
                                double visc)
{

  EngConstants["E1"] = E1;
  EngConstants["E2"] = E2;
  EngConstants["E3"] = E3;

  EngConstants["G12"] = G12;
  EngConstants["G31"] = G31;
  EngConstants["G23"] = G23;

  EngConstants["nu12"] = nu12;
  EngConstants["nu13"] = nu13;
  EngConstants["nu23"] = nu23;

  EngConstants["rho"] = rho;

  EngConstants["K1"] = K1;
  EngConstants["K2"] = K2;
  EngConstants["K3"] = K3;

  EngConstants["visc"] = visc;
  
  

  if(PlaneStress12_flag == up_to_date) PlaneStress12_flag = inited;
  if(PlaneStrain12_flag == up_to_date) PlaneStrain12_flag = inited;
  if(RotationSymmetric3_flag == up_to_date) RotationSymmetric3_flag = inited;
  if(ThreeDSolid_flag == up_to_date) ThreeDSolid_flag = inited;
  if(ThreeDShell_flag == up_to_date) ThreeDShell_flag = inited;
}

void OrthotropicMaterial::Set( string EngConst_, double value_){

  if ( EngConst_ == "E1" || EngConst_ == "E2" || EngConst_ == "E3" ||
       EngConst_ == "G12" || EngConst_ == "G31" || EngConst_ == "G23" ||
       EngConst_ == "nu12" || EngConst_ == "nu13" || EngConst_ == "nu23" ||
       EngConst_ == "rho" ||
       EngConst_ == "K1" || EngConst_ == "K2" || EngConst_ == "K3" ||
       EngConst_ == "visc" ||
       EngConst_ == "Xt" || EngConst_ == "Xc" || EngConst_ == "Yt" || EngConst_ == "Yc" ||
       EngConst_ == "S" )
    EngConstants[EngConst_] = value_;
  else if ( EngConst_ == "nu21" || EngConst_ == "nu31" || EngConst_ == "nu32" ||
        EngConst_ == "G21" || EngConst_ == "G13" || EngConst_ == "G32")
  {
    cerr << "ERROR: In OrthotropicMaterial::Set(string, double)" << endl
     << "ERROR: OrthotropicMaterial has the following data" << endl
     << "ERROR: members: E1, E2, E3, G12, G31, G23, nu12, nu13, nu23, " << endl
     << "ERROR: rho, K1, K2, K3, visc, Xt, Xc, Yt, Yc, S. " << endl
     << "ERROR: Set has been called with the following string: " << EngConst_ << endl
     << "ERROR: Although the inserted value might define the material propperly " << endl
     << "ERROR: the conversion is not yet implemented... feel free to do so !!!" << endl << endl;
    exit(1);
  }
  else {
    cerr << "ERROR: In OrthotropicMaterial::Set(string, double)" << endl
     << "ERROR: OrthotropicMaterial has the following data" << endl
     << "ERROR: members: E1, E2, E3, G12, G31, G23, nu12, nu13, nu23, " << endl
     << "ERROR: rho, K1, K2, K3, visc, Xt, Xc, Yt, Yc, S. " << endl
     << "ERROR: Set has been called with the following string: " << EngConst_ << endl << endl;
    exit(1);
  }
  if(PlaneStress12_flag == up_to_date) PlaneStress12_flag = inited;
  if(PlaneStrain12_flag == up_to_date) PlaneStrain12_flag = inited;
  if(RotationSymmetric3_flag == up_to_date) RotationSymmetric3_flag = inited;
  if(ThreeDSolid_flag == up_to_date) ThreeDSolid_flag = inited;
  if(ThreeDShell_flag == up_to_date) ThreeDShell_flag = inited;

}

void OrthotropicMaterial::SetUltimate( double Xt, double Xc, double Yt, double Yc, double S){

  if (Xt < 0. || Xc < 0. || Yt < 0. || Yc < 0. || S < 0.) {
    cerr << "ERROR: In OrthotropicMaterial::Set(string, double)" << endl
     << "ERROR: At least one ultimate load is negative! All ultimate loads" << endl
     << "ERROR: for tension and compression must be positive" << endl << endl;
    exit(1);
  }
  else {
  EngConstants["Xt"] = Xt;
  EngConstants["Xc"] = Xc;
  EngConstants["Yt"] = Yt;
  EngConstants["Yc"] = Yc;
  EngConstants["S"] = S;
  }
}

double OrthotropicMaterial::Get(string EngConst_) const{
  map<string, double>::const_iterator mapIter;
  mapIter = EngConstants.find(EngConst_);
  if ( mapIter == EngConstants.end() ) {
    cerr << "ERROR: In OrthotropicMaterial::Set(string, double)" << endl
     << "ERROR: OrthotropicMaterial has the following data" << endl
     << "ERROR: members: E1, E2, E3, G12, G31, G23, nu12, nu13, nu23, " << endl
     << "ERROR: rho, K1, K2, K3, visc, Xt, Xc, Yt, Yc, S. " << endl
     << "ERROR: Set has been called with the following string: " << EngConst_ << endl << endl;
    exit(1);
  }
  else
    return mapIter -> second;
}


//Get material properties functions
mtl::Dense_Matrix  OrthotropicMaterial::Get( StressStrainRelation type )
{
  switch (type)
    {
    case ThreeDSolid :
      {
        switch ( ThreeDSolid_flag )
          {
          case up_to_date :
            {
              return ThreeDSolid_mat;
              break;
            }
          case inited :
            {
              if( OrthotropicMaterial::Get("E1") == 0. || OrthotropicMaterial::Get("E2") == 0. || OrthotropicMaterial::Get("E3") == 0. ||
                  OrthotropicMaterial::Get("G12") == 0. || OrthotropicMaterial::Get("G31") == 0. || OrthotropicMaterial::Get("G23") == 0. ||
                  OrthotropicMaterial::Get("nu12") == 0. || OrthotropicMaterial::Get("nu13") == 0. || OrthotropicMaterial::Get("nu23") == 0.)
                  {
                        cerr << "ERROR: In OrthotropicMaterial::Get(ThreeDSolid)" << endl
                                << "ERROR: OrthotropicMaterial needs the following data" << endl
                                << "ERROR: members: E1, E2, E3, G12, G31, G23, nu12, nu13, nu23, " << endl
                                << "ERROR: Get has been called with one of the above values equal to zero" << endl << endl;
                        exit(1);
                  }

              ThreeDSolid_flag  = up_to_date;
              mtl::set_value(ThreeDSolid_mat,0.0);

              double E1 = OrthotropicMaterial::Get("E1"), E2 = OrthotropicMaterial::Get("E2"), E3 = OrthotropicMaterial::Get("E3"),
                        G12 = OrthotropicMaterial::Get("G12"), G31 = OrthotropicMaterial::Get("G31"), G23 = OrthotropicMaterial::Get("G23"),
                        nu12 = OrthotropicMaterial::Get("nu12"), nu13 = OrthotropicMaterial::Get("nu13"), nu23 = OrthotropicMaterial::Get("nu23");

              double divisor = 2*nu13*nu12*nu23*E2*E3+nu13*nu13*E2*E3+nu12*nu12*E2*E2-E1*E2+E1*nu23*nu23*E3;

              ThreeDSolid_mat(0,0) = E1*E1*(-E2+nu23*nu23*E3)/divisor;
              ThreeDSolid_mat(0,1) = -(nu12*E2+nu13*nu23*E3)*E1*E2/divisor;
              ThreeDSolid_mat(0,2) = -(nu12*nu23+nu13)*E1*E2*E3/divisor;

              ThreeDSolid_mat(1,0) = ThreeDSolid_mat(0,1);
              ThreeDSolid_mat(1,1) = (-E1+nu13*nu13*E3)*E2*E2/divisor;
              ThreeDSolid_mat(1,2) = -(nu23*E1+nu13*nu12*E2)*E2*E3/divisor;

              ThreeDSolid_mat(2,0) = ThreeDSolid_mat(0,2);
              ThreeDSolid_mat(2,1) = ThreeDSolid_mat(1,2);
              ThreeDSolid_mat(2,2) = (-E1+nu12*nu12*E2)*E2*E3/divisor;

              ThreeDSolid_mat(3,3) = G23;
              ThreeDSolid_mat(4,4) = G31;
              ThreeDSolid_mat(5,5) = G12;

              return ThreeDSolid_mat;
              break;
            }
          case not_inited :
            {
              ThreeDSolid_mat = mtl::Dense_Matrix(6, 6);
              ThreeDSolid_flag  = inited;
              return Get( ThreeDSolid );
              break;
            }
          }//end of inner switch

      }//end of case ThreeDSolid

     case PlaneStrain12 :
      {

        switch (PlaneStrain12_flag)
          {
          case up_to_date :
            {
              return PlaneStrain12_mat;
              break;
            }
          case inited :
            {

                if( OrthotropicMaterial::Get("E1") == 0. || OrthotropicMaterial::Get("E2") == 0. || OrthotropicMaterial::Get("E3") == 0. ||
                  OrthotropicMaterial::Get("G12") == 0. ||
                  OrthotropicMaterial::Get("nu12") == 0. || OrthotropicMaterial::Get("nu13") == 0. || OrthotropicMaterial::Get("nu23") == 0.)
                  {
                        cerr << "ERROR: In OrthotropicMaterial::Get(PlaneStrain12)" << endl
                                << "ERROR: OrthotropicMaterial needs the following data" << endl
                                << "ERROR: members: E1, E2, E3, G12, nu12, nu13, nu23, " << endl
                                << "ERROR: Get has been called with one of the above values equal to zero" << endl << endl;
                        exit(1);
                  }
              PlaneStrain12_flag = up_to_date;
              mtl::set_value(PlaneStrain12_mat,0.0);

              double E1 = OrthotropicMaterial::Get("E1"), E2 = OrthotropicMaterial::Get("E2"), E3 = OrthotropicMaterial::Get("E3"),
                        G12 = OrthotropicMaterial::Get("G12"),
                        nu12 = OrthotropicMaterial::Get("nu12"), nu13 = OrthotropicMaterial::Get("nu13"), nu23 = OrthotropicMaterial::Get("nu23");

              double divisor = 2*nu13*nu12*nu23*E2*E3+nu13*nu13*E2*E3+nu12*nu12*E2*E2-E1*E2+E1*nu23*nu23*E3;

              PlaneStrain12_mat(0,0) = E1*E1*(-E2+nu23*nu23*E3)/divisor;
              PlaneStrain12_mat(0,1) = -(nu12*E2+nu13*nu23*E3)*E1*E2/divisor;

              PlaneStrain12_mat(1,0) = PlaneStrain12_mat(0,1);
              PlaneStrain12_mat(1,1) = (-E1+nu13*nu13*E3)*E2*E2/divisor;

              PlaneStrain12_mat(2,2) = G12;

              return PlaneStrain12_mat;
              break;
            }
          case not_inited :
            {
              PlaneStrain12_mat = mtl::Dense_Matrix(3,3);
              PlaneStrain12_flag = inited;
              return Get(PlaneStrain12);
              break;
            }
          }//end of inner switch
      } //end of case planestress12

     case PlaneStress12 :
      {

        switch (PlaneStress12_flag)
          {
          case up_to_date :
            {
              return PlaneStress12_mat;
              break;
            }
          case inited :
            {
                if( OrthotropicMaterial::Get("E1") == 0. || OrthotropicMaterial::Get("E2") == 0. ||
                OrthotropicMaterial::Get("G12") == 0. || OrthotropicMaterial::Get("nu12") == 0. )
                {
                cerr << "ERROR: In OrthotropicMaterial::Get(PlaneStress12)" << endl
                        << "ERROR: OrthotropicMaterial needs the following data" << endl
                        << "ERROR: members: E1, E2, G12, nu12 " << endl
                        << "ERROR: Get has been called with one of the above values equal to zero" << 
                        OrthotropicMaterial::Get("E1") << " " << OrthotropicMaterial::Get("E2") << " " <<
                        OrthotropicMaterial::Get("G12") << " " <<
                        OrthotropicMaterial::Get("nu12") << endl << endl;
                exit(1);
                }
              PlaneStress12_flag = up_to_date;
              mtl::set_value(PlaneStress12_mat,0.0);

              double E1 = OrthotropicMaterial::Get("E1"), E2 = OrthotropicMaterial::Get("E2"),
                        G12 = OrthotropicMaterial::Get("G12"), nu12 = OrthotropicMaterial::Get("nu12");

              double divisor = E1-nu12*nu12*E2;

              PlaneStress12_mat(0,0) = E1*E1/divisor;
              PlaneStress12_mat(0,1) = nu12*E1*E2/divisor;

              PlaneStress12_mat(1,0) = PlaneStress12_mat(0,1);
              PlaneStress12_mat(1,1) = E2*E1/divisor;

              PlaneStress12_mat(2,2) = G12;

              return PlaneStress12_mat;
              break;
            }
          case not_inited :
            {
              PlaneStress12_mat = mtl::Dense_Matrix(3,3);
              PlaneStress12_flag = inited;
              return Get(PlaneStress12);
              break;
            }
          }//end of inner switch
      } //end of case planestress12

    case ThreeDShell :
      {
        switch ( ThreeDShell_flag )
          {
          case up_to_date :
            {
              return ThreeDShell_mat;
              break;
            }
          case inited :
            {
              //This stressstrainrelation is not to be used as computed here
              //The diagonal values 4 and 5 are to be scaled depending on the
              //area and the thickness of the element. These parameters can not
              //be accessed in Material but only from the element..

              ThreeDShell_flag  = up_to_date;

              mtl::set_value(ThreeDShell_mat,0.0);

              mtl::Dense_Matrix PlaneStress12Matrix(3,3);
              PlaneStress12Matrix = Get(PlaneStress12);

              ThreeDShell_mat(0,0) = PlaneStress12Matrix(0,0);
              ThreeDShell_mat(1,1) = PlaneStress12Matrix(1,1);

              ThreeDShell_mat(0,1) = PlaneStress12Matrix(0,1);
              ThreeDShell_mat(1,0) = PlaneStress12Matrix(1,0);

              ThreeDShell_mat(3,3) = PlaneStress12Matrix(2,2);
              ThreeDShell_mat(4,4) = PlaneStress12Matrix(2,2);
              ThreeDShell_mat(5,5) = PlaneStress12Matrix(2,2);


              return ThreeDShell_mat;
              break;
            }
          case not_inited :
            {
              ThreeDShell_mat = mtl::Dense_Matrix(6, 6);
              ThreeDShell_flag  = inited;
              return Get( ThreeDShell );
              break;
            }
          }//end of inner switch

      }//end of case ThreeDShell

      
      default :
      {
          cerr << "ERROR: In OrthotropicMaterial::Get(StressStrainRelation)" << endl;
          cerr << "ERROR: The requested stress-strain-relation is not available" << endl;
          cerr << "ERROR: It is either not implemented or makes no sense to use...." << endl;
          cerr << "ERROR: A 3x3 dummy matrix is returned" << endl;

          mtl::Dense_Matrix dummy;
          mtl::set_value(dummy,0.0);
          return dummy;
      }


    }//end of outer switch
}


void OrthotropicMaterial::print(ostream& stream) const
{
  stream << " Orthotropic: ";
  stream.setf(ios::right, ios::adjustfield);
  stream.precision(8);

  stream.width(15); stream << Get("E1");
  stream.width(15); stream << Get("E2");
  stream.width(15); stream << Get("E3");

  stream.width(15); stream << Get("G12");
  stream.width(15); stream << Get("G31");
  stream.width(15); stream << Get("G23");

  stream.width(15); stream << Get("nu12");
  stream.width(15); stream << Get("nu13");
  stream.width(15); stream << Get("nu23");

  stream.width(15); stream << Get("rho");

  stream.width(15); stream << Get("K1");
  stream.width(15); stream << Get("K2");
  stream.width(15); stream << Get("K3");

  stream.width(15); stream << Get("visc");
  
  stream.width(15); stream << Get("Xt");
  stream.width(15); stream << Get("Xc");
  stream.width(15); stream << Get("Yt");
  stream.width(15); stream << Get("Yc");
  stream.width(15); stream << Get("S");

}



/////////////////////////////////////////////////
/// constructors of class 'IsotropicMaterial'
/////////////////////////////////////////////////

IsotropicMaterial::IsotropicMaterial(){
  Set(0., 0., 0., 0., 0.);

  PlaneStress12_flag = not_inited;
  PlaneStrain12_flag = not_inited;
  RotationSymmetric3_flag = not_inited;
  ThreeDSolid_flag = not_inited;
  ThreeDShell_flag = not_inited;
}

IsotropicMaterial::IsotropicMaterial(double E, double nu, double rho, double K, double visc){
  Set( E, nu, rho, K, visc );

  PlaneStress12_flag = not_inited;
  PlaneStrain12_flag = not_inited;
  RotationSymmetric3_flag = not_inited;
  ThreeDSolid_flag = not_inited;
  ThreeDShell_flag = not_inited;
}

/////////////////////////////////////////////////
/// functions member of class 'IsotropicMaterial'
/////////////////////////////////////////////////

void IsotropicMaterial::Set(double E, double nu, double rho, double K, double visc){

  Set("nu", nu);
  Set("E", E);
  Set("rho", rho);
  Set("K", K);
  Set("visc", visc);

  if(PlaneStress12_flag == up_to_date) PlaneStress12_flag = inited;
  if(PlaneStrain12_flag == up_to_date) PlaneStrain12_flag = inited;
  if(RotationSymmetric3_flag == up_to_date) RotationSymmetric3_flag = inited;
  if(ThreeDSolid_flag == up_to_date) ThreeDSolid_flag = inited;
  if(ThreeDShell_flag == up_to_date) ThreeDShell_flag = inited;
}

void IsotropicMaterial::Set(string EngConst, double value ){

  if ( EngConst == "E" )
        EngConstants["E1"] = EngConstants["E2"] = EngConstants["E3"] = value;
  else if ( EngConst == "nu" )
        EngConstants["nu12"] = EngConstants["nu13"] = EngConstants["nu23"] = value;
  else if ( EngConst == "K" )
        EngConstants["K1"] = EngConstants["K2"] = EngConstants["K3"] = value;
  else if ( EngConst == "rho" || EngConst == "visc" )
        EngConstants[EngConst] = value;
  else if ( EngConst == "G" )
  {
    cerr << "ERROR: In IsotropicMaterial::Set(string, double)" << endl
     << "ERROR: IsotropicMaterial has the following data" << endl
     << "ERROR: members: E, nu, rho, K, visc." << endl
     << "ERROR: Set has been called with the following string: " << EngConst << endl
     << "ERROR: Although the inserted value might define the material propperly " << endl
     << "ERROR: the conversion is not yet implemented... feel free to do so !!!" << endl << endl;
    exit(1);
  }
  else {
    cerr << "ERROR: In IsotropicMaterial::Set(string, double)" << endl
     << "ERROR: IsotropicMaterial has the following data" << endl
     << "ERROR: members: E, nu, rho, K, visc. " << endl
     << "ERROR: Set has been called with the following string: " << EngConst << endl << endl;
    exit(1);
  }

  //Reassign G it might have changed
  SetG();


  if(PlaneStress12_flag == up_to_date) PlaneStress12_flag = inited;
  if(PlaneStrain12_flag == up_to_date) PlaneStrain12_flag = inited;
  if(RotationSymmetric3_flag == up_to_date) RotationSymmetric3_flag = inited;
  if(ThreeDSolid_flag == up_to_date) ThreeDSolid_flag = inited;
  if(ThreeDShell_flag == up_to_date) ThreeDShell_flag = inited;
}


double IsotropicMaterial::Get(string EngConst) const{
  std::string the_const;
  if (EngConst == "E" )
    the_const = "E1";
  else if (EngConst == "nu" || EngConst == "nu12" )
    the_const = "nu12";
  else if (EngConst == "rho")
    the_const = "rho";
  else if (EngConst == "G" )
    the_const = "G12";
  else if (EngConst == "K")
    the_const = "K1";
  else if (EngConst == "visc")
    the_const = "visc";

  map<string, double>::const_iterator mapIter;
  mapIter = EngConstants.find(the_const);
  if ( mapIter == EngConstants.end() ) {
    cerr << "ERROR: In IsotropicMaterial::Get(string)" << endl
     << "ERROR: IsotropicMaterial has the following data" << endl
     << "ERROR: members: E, nu, rho, K, visc. " << endl
     << "ERROR: Get has been called with the following string: " << EngConst << endl << endl;
    exit(1);
  }
  else
    return mapIter -> second;
}

void IsotropicMaterial::SetG()
{
        if ( Get("E") != 0. && Get("nu") != 0. )
                EngConstants["G12"] = EngConstants["G31"] = EngConstants["G23"] = Get("E")/2./(1.+Get("nu"));
        else
                EngConstants["G12"] = EngConstants["G31"] = EngConstants["G23"] = 0.;
}

void IsotropicMaterial::print(ostream& stream) const
{
  stream << " Isotropic Material  : ";
  stream.setf(ios::right, ios::adjustfield);
  stream.precision(9);
  stream.width(15);     stream << Get("E");
  stream.width(15);     stream << Get("nu");
  stream.width(15);     stream << Get("rho");
  stream.width(15);     stream << Get("K");
  stream.width(15);     stream << Get("visc");
}



///////////////////////////////////////////////////////////////
/// functions member of class 'TransverseIsotropicMaterial23'
//////////////////////////////////////////////////////////////

TransverseIsotropicMaterial23::TransverseIsotropicMaterial23(){

  Set(0.0, 0.0,
        0.0, 0.0,
        0.0,
        0.0,
        0.0, 0.0,
        0.0);

  PlaneStress12_flag = not_inited;
  PlaneStrain12_flag = not_inited;
  RotationSymmetric3_flag = not_inited;
  ThreeDSolid_flag = not_inited;
  ThreeDShell_flag = not_inited;
}

TransverseIsotropicMaterial23::TransverseIsotropicMaterial23(double E1, double E2,
                                                         double G12, double G23,
                                                         double nu12,
                                                         double rho,
                                                         double K1, double K2,
                                                         double visc){

  Set(E1, E2, G12, G23, nu12, rho, K1, K2, visc);

  PlaneStress12_flag = not_inited;
  PlaneStrain12_flag = not_inited;
  RotationSymmetric3_flag = not_inited;
  ThreeDSolid_flag = not_inited;
  ThreeDShell_flag = not_inited;
}


void TransverseIsotropicMaterial23::Set(double E1, double E2,
                                        double G12, double G23,
                                        double nu12,
                                        double rho,
                                        double K1, double K2,
                                        double visc){

  Set("E1",E1);
  Set("E2",E2);
  Set("G12",G12);
  Set("G23",G23);
  Set("nu12",nu12);
  Set("rho",rho);
  Set("K1",K1);
  Set("K2",K2);
  Set("visc", visc);

  if(PlaneStress12_flag == up_to_date) PlaneStress12_flag = inited;
  if(PlaneStrain12_flag == up_to_date) PlaneStrain12_flag = inited;
  if(RotationSymmetric3_flag == up_to_date) RotationSymmetric3_flag = inited;
  if(ThreeDSolid_flag == up_to_date) ThreeDSolid_flag = inited;
  if(ThreeDShell_flag == up_to_date) ThreeDShell_flag = inited;
}

void TransverseIsotropicMaterial23::Set( string EngConst, double value){

  if ( EngConst == "E1" || EngConst == "G23" ||
        EngConst == "rho" || EngConst == "K1" || EngConst == "visc" ||
        EngConst == "Xt" || EngConst == "Xc" || EngConst == "Yt" || EngConst == "Yc" || EngConst == "S")
        EngConstants[EngConst] = value;
  else if ( EngConst == "G12" )
        EngConstants["G12"] = EngConstants["G31"] = value;
  else if ( EngConst == "K2" )
        EngConstants["K2"] = EngConstants["K3"] = value;
  else if ( EngConst == "nu12" )
        EngConstants["nu12"] = EngConstants["nu13"] = value;
  else if ( EngConst == "E2" )
        EngConstants["E2"] = EngConstants["E3"] = value;
  else if ( EngConst == "nu23" )
  {
    cerr << "ERROR: In TransverseIsotropicMaterial23::Set(string, double)" << endl
     << "ERROR: TransverseIsotropicMaterial23 has the following data" << endl
     << "ERROR: members: E1, E2, G12, G23, nu12, rho, K1, K2, visc, Xt, Xc, Yt, Yc, S." << endl
     << "ERROR: Set has been called with the following string: " << EngConst << endl
     << "ERROR: Although the inserted value might define the material propperly " << endl
     << "ERROR: the conversion is not yet implemented... feel free to do so !!!" << endl << endl;
    exit(1);
   }
  else {
    cerr << "ERROR: In TransverseIsotropicMaterial23::Set(string, double)" << endl;
    cerr << "ERROR: TransverseIsotropicMaterial23 has the following data" << endl;
    cerr << "ERROR: members: E1, E2, G12, G23, nu12, nu23, rho, K1, K2, K3, visc, Xt, Xc, Yt, Yc, S. " << endl;
    cerr << "ERROR: Set has been called with the following string: " << EngConst << endl << endl;
    exit(1);
  }

  Setnu23();

  if(PlaneStress12_flag == up_to_date) PlaneStress12_flag = inited;
  if(PlaneStrain12_flag == up_to_date) PlaneStrain12_flag = inited;
  if(RotationSymmetric3_flag == up_to_date) RotationSymmetric3_flag = inited;
  if(ThreeDSolid_flag == up_to_date) ThreeDSolid_flag = inited;
  if(ThreeDShell_flag == up_to_date) ThreeDShell_flag = inited;
}

double TransverseIsotropicMaterial23::Get(string EngConst_) const{
  map<string, double>::const_iterator mapIter;
  mapIter = EngConstants.find(EngConst_);
  if ( mapIter == EngConstants.end() ) {
    cerr << "ERROR: In TransverseIsotropicMaterial23::Get(string)" << endl;
    cerr << "ERROR: TransverseIsotropicMaterial23 has the following data" << endl;
    cerr << "ERROR: members: E1, E2, G12, G23, nu12, rho, K1, K2, K3, visc. " << endl;
    cerr << "ERROR: Get has been called with the following string: " << EngConst_ << endl << endl;
    exit(1);
  }
  else
    return mapIter -> second;
}

void TransverseIsotropicMaterial23::Setnu23()
{
        if ( Get("E2") != 0. && Get("G23") != 0. )
                EngConstants["nu23"] = (Get("E2")/2./Get("G23"))-1.;
        else
                EngConstants["nu23"] = 0.;
}

void TransverseIsotropicMaterial23::print(ostream& stream) const
{
  stream << " TransverseIsotropic23: ";
  stream.setf(ios::right, ios::adjustfield);
  stream.precision(8);

  stream.width(15); stream << Get("E1");
  stream.width(15); stream << Get("E2");

  stream.width(15); stream << Get("G12");
  stream.width(15); stream << Get("G23");

  stream.width(15); stream << Get("nu12");

  stream.width(15); stream << Get("rho");

  stream.width(15); stream << Get("K1");
  stream.width(15); stream << Get("K2");

  stream.width(15); stream << Get("visc");
  
  stream.width(15); stream << Get("Xt");
  stream.width(15); stream << Get("Xc");
  stream.width(15); stream << Get("Yt");
  stream.width(15); stream << Get("Yc");
  stream.width(15); stream << Get("S");

}
/*
mtl::Dense_Matrix TransverseIsotropicMaterial23::Get(StressStrainRelation type)
{
   switch (type)
    {
    case PlaneStrain12 :
      {
      mtl::Dense_Matrix Q(3,3);

//   std::cout << "The Transformation Matrix" << std::endl;
//   mtl::print_all_matrix(T);

      mtl::set_value(Q,0.);

      if ( Get("E1") != 0.0 && Get("E2") != 0.0 && Get("nu12") != 0.0 && Get("G12") != 0.0 ){
          double E1 = Get("E1"); double E2 = Get("E2"); 
          double nu12 = Get("nu12"); double G12 = Get("G12");
      
          Q(0,0) = E1 / (1 - nu12*nu12);
          Q(0,1) = nu12 * E2 / ( 1 - nu12*nu12 ); Q(1,0) = Q(0,1);
          Q(1,1) = E2 / (1 - nu12*nu12);
          Q(2,2) = G12;
      }
      else
      {
        cerr << "ERROR: In TransverseIsotropicMaterial23::GetStressStrainRelation" << endl;
        cerr << "ERROR: TransverseIsotropicMaterial23 has the following data" << endl;
        cerr << "ERROR: members: E1, E2, G12, G23, nu12, rho, K1, K2, K3, visc. " << endl;
        exit(1);
      }
   
      return Q;
    }
   }
}
*/

/////////////////////////////////////////////////
/// constructors of class 'WeaveMaterial'
/////////////////////////////////////////////////

WeaveMaterial::WeaveMaterial(){
  Set(0., 0., 0., 0., 0.,0.);

  SingleLayerMaterial.Set(0.,0.,0.,
                          0.,0.,0.,
                          0.,0.,0.);

  PlaneStress12_flag = not_inited;
  PlaneStrain12_flag = not_inited;
  RotationSymmetric3_flag = not_inited;
  ThreeDSolid_flag = not_inited;
  ThreeDShell_flag = not_inited;

  SetSingleLayerMaterial();
}

WeaveMaterial::WeaveMaterial(double E1, double E2,
                        double G12, double nu12,
                        double rho,
                        double E_Ratio)
{
  Set( E1, E2, G12, nu12, rho, E_Ratio);

  PlaneStress12_flag = not_inited;
  PlaneStrain12_flag = not_inited;
  RotationSymmetric3_flag = not_inited;
  ThreeDSolid_flag = not_inited;
  ThreeDShell_flag = not_inited;

  SetSingleLayerMaterial();
}

/////////////////////////////////////////////////
/// functions member of class 'WeaveMaterial'
/////////////////////////////////////////////////

void WeaveMaterial::Set(double E1, double E2,
                        double G12, double nu12,
                        double rho,
                        double E_Ratio)
{
  if (E1 == E2)
  {
        Set("E1", E1);
        Set("E2", E2);
        Set("G12", G12);
        Set("nu12", nu12);
        Set("rho", rho);
        Set("E_Ratio",E_Ratio);


        if(PlaneStress12_flag == up_to_date) PlaneStress12_flag = inited;
        if(PlaneStrain12_flag == up_to_date) PlaneStrain12_flag = inited;
        if(RotationSymmetric3_flag == up_to_date) RotationSymmetric3_flag = inited;
        if(ThreeDSolid_flag == up_to_date) ThreeDSolid_flag = inited;
        if(ThreeDShell_flag == up_to_date) ThreeDShell_flag = inited;
  }
  else
  {
     cerr << "ERROR: In WeaveMaterial::Set(double, double...)" << endl
        << "ERROR: WeaveMaterial is implemented to have E1 = E2" << endl
        << "ERROR: Please correct it..." << endl << endl;
    exit(1);
  }
  SetSingleLayerMaterial();
}

void WeaveMaterial::Set(string EngConst, double value ){

  if ( EngConst == "E1" || EngConst == "E2" )
        EngConstants["E1"] = EngConstants["E2"] = value;
  else if ( EngConst == "nu12" )
        EngConstants["nu12"] = value;
  else if ( EngConst == "G12" )
        EngConstants["G12"] =  value;
  else if ( EngConst == "rho" )
        EngConstants["rho"] = value;
  else if ( EngConst == "E_Ratio" )
        E_Ratio = value;
  else {
    cerr << "ERROR: In WeaveMaterial::Set(string, double)" << endl
     << "ERROR: WeaveMaterial has the following data" << endl
     << "ERROR: members: E1, E2, nu12, G12, E_Ratio, rho" << endl
     << "ERROR: Set has been called with the following string: " << EngConst << endl << endl;
    exit(1);
  }


  if(PlaneStress12_flag == up_to_date) PlaneStress12_flag = inited;
  if(PlaneStrain12_flag == up_to_date) PlaneStrain12_flag = inited;
  if(RotationSymmetric3_flag == up_to_date) RotationSymmetric3_flag = inited;
  if(ThreeDSolid_flag == up_to_date) ThreeDSolid_flag = inited;
  if(ThreeDShell_flag == up_to_date) ThreeDShell_flag = inited;

  SetSingleLayerMaterial();
}


double WeaveMaterial::Get(string EngConst) const{

  if (EngConst == "E_Ratio")
        return E_Ratio;
  else
  {
        std::string the_const;
        if (EngConst == "E1")
        the_const = "E1";
        else if (EngConst == "E2")
        the_const = "E2";
        else if (EngConst == "nu12")
        the_const = "nu12";
        else if (EngConst == "rho")
        the_const = "rho";
        else if (EngConst == "G12")
        the_const = "G12";
        map<string, double>::const_iterator mapIter;
        mapIter = EngConstants.find(the_const);
        if ( mapIter == EngConstants.end() ) {
        cerr << "ERROR: In WeaveMaterial::Get(string)" << endl
        << "ERROR: WeaveMaterial has the following data" << endl
        << "ERROR: members: E1, E2, nu12, G12, E_Ratio, rho" << endl
        << "ERROR: Get has been called with the following string: " << EngConst << endl << endl;
        exit(1);
        }
        else
        return mapIter -> second;
  }
}

void WeaveMaterial::SetSingleLayerMaterial()
{

        double E1 = 0., E2 = 0., nu12 = 0., G12 = 0.;

        if ( Get("E1") != 0. && Get("E2") != 0. &&
             Get("nu12") != 0. && Get("G12") != 0.)
        {

                mtl::Dense_Matrix PSW = OrthotropicMaterial::Get(PlaneStress12);

                nu12 = PSW(0,1)*(1. + E_Ratio)/2./PSW(0,0)/E_Ratio;
                E1 = PSW(0,0)*(1.-E_Ratio*nu12*nu12)/(1.+E_Ratio);
                E2 = E_Ratio*E1;
                G12 = PSW(2,2)/2.;

        }

        SingleLayerMaterial.Set(E1,E2,G12,0.,nu12,0.,0.,0.,0.);
}

mtl::Dense_Matrix WeaveMaterial::GetShearedPlaneStress12(double ShearAngle)
{
   mtl::Dense_Matrix Q(3,3), TransfQ(3,3), T(3,3);
   T = GetPlaneMaterialTransformationMatrix(ShearAngle);

//   std::cout << "The Transformation Matrix" << std::endl;
//   mtl::print_all_matrix(T);

   mtl::set_value(Q,0.);
   mtl::set_value(TransfQ,0.);

   if ( SingleLayerMaterial.Get("E1") != 0.  && SingleLayerMaterial.Get("E2") != 0. &&
        SingleLayerMaterial.Get("nu12") != 0. && SingleLayerMaterial.Get("G12") != 0. )
   {
        double div = 1.-SingleLayerMaterial.Get("E2")/SingleLayerMaterial.Get("E1")
                        *SingleLayerMaterial.Get("nu12")*SingleLayerMaterial.Get("nu12");
        Q(0,0) = SingleLayerMaterial.Get("E1")/div;
        Q(1,1) = SingleLayerMaterial.Get("E2")/div;
        Q(0,1) = Q(1,0) = SingleLayerMaterial.Get("E2")*SingleLayerMaterial.Get("nu12")/div;
        Q(2,2) = SingleLayerMaterial.Get("G12");

        mtl::copy(Q,TransfQ);
        mult_AT_B_A_add(T,Q,TransfQ);

//        std::cout << " The transformed weavematerial with " << ShearAngle << " rad angle" << std::endl;
//        mtl::print_all_matrix(TransfQ);
   }
   else
   {
        cerr << "ERROR: In WeaveMaterial::GetShearedPlaneStress" << endl
        << "ERROR: The SingleLayerMaterial is not set properly, one of the engineering constants is zero" << endl << endl;
        exit(1);
   }

   return TransfQ;
}

mtl::Dense_Matrix WeaveMaterial::GetShearedThreeDShell(double ShearAngle)
{
   mtl::Dense_Matrix ThreeD(6,6), PlaneStress12Matrix(3,3);

//   std::cout << "The unset PlaneStressMatrix " << std::endl;
//   mtl::print_all_matrix(PlaneStress12Matrix);

   PlaneStress12Matrix = GetShearedPlaneStress12(ShearAngle);

//   std::cout << "The set PlaneStressMatrix " << std::endl;
//   mtl::print_all_matrix(PlaneStress12Matrix);

   mtl::set_value(ThreeD,0.);

   ThreeD(0,0) = PlaneStress12Matrix(0,0);
   ThreeD(1,1) = PlaneStress12Matrix(1,1);

   ThreeD(0,1) = PlaneStress12Matrix(0,1);
   ThreeD(1,0) = PlaneStress12Matrix(1,0);

   ThreeD(3,3) = PlaneStress12Matrix(2,2);
   ThreeD(4,4) = PlaneStress12Matrix(2,2);
   ThreeD(5,5) = PlaneStress12Matrix(2,2);

   return ThreeD;
}

void WeaveMaterial::print(ostream& stream) const
{
  stream << " Weave Material  : ";
  stream.setf(ios::right, ios::adjustfield);
  stream.precision(9);
  stream.width(15);     stream << Get("E1");
  stream.width(15);     stream << Get("E2");
  stream.width(15);     stream << Get("nu12");
  stream.width(15);     stream << Get("G12");
  stream.width(15);     stream << Get("E_Ratio");
  stream.width(15);     stream << Get("rho");
}




//////////////////////////////////////
/// functions not member of any class
//////////////////////////////////////
ostream& fe_base::operator<< (ostream& stream, const Material& m){
  m.print(stream);

  return stream;
}

mtl::Dense_Matrix fe_base::GetPlaneMaterialTransformationMatrix( double  theta) //rotation angle about 3-direction in rads! This matrix is for stresses valid only!!!
{       //! returns the inverse AND transponed plane transformation matrix, which is used in the mult_AT_B_A_add function

        double c = cos(theta), s = sin(theta);

        mtl::Dense_Matrix T(3,3);
        // T^(-T)
        T(0,0) = c*c; T(0,1) = s*s; T(0,2) = s*c;
        T(1,0) = s*s; T(1,1) = c*c; T(1,2) = -s*c;
        T(2,0) = -2*c*s; T(2,1) = 2*c*s; T(2,2) = c*c-s*s;

        return T;
}

mtl::Dense_Matrix fe_base::GetPlaneMaterialTransformationMatrixStrains( double  theta) //rotation angle about 3-direction in rads. This matrix is for strains valid only!!!
{       //! returns the inverse AND transponed plane transformation matrix, which is used in the mult_AT_B_A_add function

        double c = cos(theta), s = sin(theta);

        mtl::Dense_Matrix T(3,3);
        // R*T^(-1)R^(-1) 
        T(0,0) = c*c; T(0,1) = s*s; T(0,2) = -s*c;
        T(1,0) = s*s; T(1,1) = c*c; T(1,2) = s*c;
        T(2,0) = 2*c*s; T(2,1) = -2*c*s; T(2,2) = c*c-s*s;

        return T;
}


mtl::Dense_Matrix fe_base::Get3DMaterialTransformationMatrix( double  theta) //rotation angle about 3-direction in rads
{       // returns the 3-dim transformation matrix

        double c = cos(theta), s = sin(theta);

        mtl::Dense_Matrix T3D(6,6);

        T3D(0,0) = c*c; T3D(0,1) = s*s; T3D(0,3) = 2*s*c;
        T3D(1,0) = s*s; T3D(1,1) = c*c; T3D(1,3) = -2*s*c;
        T3D(2,2) = 1;
        T3D(4,4) = c; T3D(4,5) = -s;
        T3D(5,4) = s; T3D(5,5) = c;
        T3D(3,0) = -c*s; T3D(3,1) = c*s; T3D(3,3) = c*c-s*s;

        return T3D;
}

