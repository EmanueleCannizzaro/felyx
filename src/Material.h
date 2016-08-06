//-----------------------------------------------------------------------------
// Material.h
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


#ifndef Material_h
#define Material_h Material_h

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "mtl/mtl_felyx_utils.h"

namespace fe_base{      //begin of namespace fe_base

mtl::Dense_Matrix fe_base::GetPlaneMaterialTransformationMatrix( double );
mtl::Dense_Matrix fe_base::GetPlaneMaterialTransformationMatrixStrains( double );

mtl::Dense_Matrix fe_base::Get3DMaterialTransformationMatrix( double );

  ////
  //// Abstract base class for materials
  ////

  class Material{
  public:

    //! Enumeration type for the stress-strain-relation.
  enum StressStrainRelation { PlaneStress12,
                                PlaneStrain12,
                                RotationSymmetric3,
                                ThreeDSolid,
                                ThreeDShell};

  //! Enumeration type for the Stress-Strain Matrix Status.
  enum StressStrain_flag{ up_to_date, inited, not_inited };

    // Constructors:
    //! Standard constructor
    Material()                  {}
    //Material(double a )               :rho(a) {}
    //! Standard destructor
    virtual ~Material()         {}

    // "Virtual Constructor" Idiom functions
    virtual Material* CreateNew()       const =0;
    virtual Material* Clone()           const =0;

    // Other
    friend std::ostream& operator<<(std::ostream&, const Material&); // overloaded << operator


    //Get material properties functions
    virtual mtl::Dense_Matrix  Get( StressStrainRelation ) = 0;
    virtual double Get( std::string ) const;
    virtual std::string ClassName() const = 0;

    //Set material properties functions
    virtual void Set(std::string, double);

  protected:
    virtual void print(std::ostream&) const =0; // virtual print-function

  }; // of class material

  ////
  //// Class for orthotropic materials, derived from "Material"
  //// serving as base class for all other linear elastic materials
  ////

  class OrthotropicMaterial:public Material{

  public:
    // Constructors:
    //! Standard constructor
    OrthotropicMaterial();

    OrthotropicMaterial(double E1, double E2, double E3,
                        double G12, double G31, double G23,
                        double nu12, double nu13, double nu23,
                        double rho,
                        double K1, double K2, double K3,
                        double visc);

    //! Standard destructor
    virtual ~OrthotropicMaterial() {}

    //! "Virtual Constructor" Idiom function.
    virtual OrthotropicMaterial* CreateNew()    const   { return new OrthotropicMaterial(); }
    //! "Virtual Constructor" Idiom function.
    virtual OrthotropicMaterial* Clone()        const   { return new OrthotropicMaterial(*this); }

    // Set functions
    void Set(double, double, double,
             double, double, double,
             double, double, double,
             double,
             double, double, double,
             double);
    virtual void Set(std::string, double); // Set single EngConstant
    
    virtual void SetUltimate(double, double, double, double, double);

    //!Get material properties function.
    virtual mtl::Dense_Matrix  Get( StressStrainRelation );
    virtual double Get(std::string) const;
    virtual std::string ClassName() const { return "OrthotropicMaterial"; };

  //protected:
    //! Print-function
    virtual void print(std::ostream&) const;

    // Data members:
    map<std::string, double> EngConstants;

    StressStrain_flag PlaneStress12_flag;
    StressStrain_flag PlaneStrain12_flag;
    StressStrain_flag RotationSymmetric3_flag;
    StressStrain_flag ThreeDSolid_flag;
    StressStrain_flag ThreeDShell_flag;

    mtl::Dense_Matrix PlaneStress12_mat;
    mtl::Dense_Matrix PlaneStrain12_mat;
    mtl::Dense_Matrix RotationSymmetric3_mat;
    mtl::Dense_Matrix ThreeDSolid_mat;
    mtl::Dense_Matrix ThreeDShell_mat;


  private:




  }; //of class Orthotropic material





  ////
  //// Class for isotropic materials, derived from "OrthotropicMaterial"
  ////
  class IsotropicMaterial:public OrthotropicMaterial{

  public:
    // Constructors:
    //! Standard constructor
    IsotropicMaterial();
    IsotropicMaterial(double E, double nu, double rho=0., double K=0., double visc=0.);
    //! Standard destructor
    virtual ~IsotropicMaterial() {}

    // "Virtual Constructor" Idiom functions
    //!Virtual constructor that creates an empty element of the type from which it is called.
    virtual IsotropicMaterial* CreateNew()      const   { return new IsotropicMaterial(); }
    //!Virtual constructor that clones the element from which it is called.
    virtual IsotropicMaterial* Clone()          const   { return new IsotropicMaterial(*this); }

    // Set functions
    void Set( double E, double nu, double rho=0., double K=0., double visc=0.); //Set ( E, nu, rho ) - rho optional
    virtual void Set( std::string, double );          //Set single EngConsts

    //A function which reassigns G to the new propper value or
    //to  0. if one of E or nu are 0
    void SetG();

    //Get material properties functions
    virtual double Get(std::string) const;
    virtual std::string ClassName() const  { return "IsotropicMaterial"; };

  protected:
    // print-function
    virtual void print(std::ostream&) const;

  }; //of class isotropicmaterial



  ////
  //// Class for transverse isotropic materials with the plane
  //// of isotropy in the 23 plane, derived from "OrthotropicMaterial"
  //// This represents a UD-Layer with the fibers in the 1 direction
  ////

  class TransverseIsotropicMaterial23:public OrthotropicMaterial{

  public:
    // Constructors:
    //! Standard constructor
    TransverseIsotropicMaterial23();

    TransverseIsotropicMaterial23(double E1, double E2,
                                  double G12, double G23,
                                  double nu12,
                                  double rho = 0.,
                                  double K1 = 0., double K2 = 0.,
                                  double visc = 0.);

    //! Standard destructor
    virtual ~TransverseIsotropicMaterial23() {}

    //! "Virtual Constructor" Idiom function.
    virtual TransverseIsotropicMaterial23* CreateNew()  const   { return new TransverseIsotropicMaterial23(); }
    //! "Virtual Constructor" Idiom function.
    virtual TransverseIsotropicMaterial23* Clone()              const   { return new TransverseIsotropicMaterial23(*this); }

    // Set functions
    void Set(double, double,
             double, double,
             double,
             double,
             double, double,
             double);
    virtual void Set(std::string, double); // Set single EngConstant

    //A function which reassigns nu23 to the new propper value or
    //to  0. if one of G23 or E2 are 0
    void Setnu23();
    
    //!Get material properties function.
    virtual double Get(std::string) const;
    virtual std::string ClassName() const { return "TransverseIsotropicMaterial23"; };

  protected:
    //! Print-function
    virtual void print(std::ostream&) const;

  }; //of class TransverseIsotropicMaterial23


  /*!   Class representing a woven material. This is to be used for problems where draping
        analysis is done in advance. Every shell element will possess a unique set of
        MaterialOrientations. The material will change its properties with the shear angle.
        When only the four engineering constants used to compute the Q-Matrix for the plane
        stress state are given an assumption for the transverse direction of one material is made
        based on the E-modul in fiber directions.
        Only the plane stress state is to be used for this material.
        The E_Ratio must be set based on the type of material (for carbon/epoxy 1/10 and
        for glass/epoxy 1/4 are good values)
        It represents the ratio of Youngs Moduli E2/E1 for the SingleLayerMaterial
  */
  class WeaveMaterial:public OrthotropicMaterial
  {
  public:
    // Constructors:
    //! Standard constructor
    WeaveMaterial();
    WeaveMaterial(double E1, double E2,
                        double G12, double nu12,
                        double rho=0.,
                        double E_Ratio = 0.);
    //! Standard destructor
    virtual ~WeaveMaterial() {}

    // "Virtual Constructor" Idiom functions
    //!Virtual constructor that creates an empty element of the type from which it is called.
    virtual WeaveMaterial* CreateNew()  const   { return new WeaveMaterial(); }
    //!Virtual constructor that clones the element from which it is called.
    virtual WeaveMaterial* Clone()              const   { return new WeaveMaterial(*this); }

    // Set functions
    void Set(double E1, double E2,
                double G12, double nu12,
                double rho=0.,
                double E_Ratio = 0.);
    virtual void Set( std::string, double );          //Set single EngConsts
    //!Function calculating the single LayerMaterial
    void SetSingleLayerMaterial();

    //Get material properties functions
    virtual double Get(std::string) const;
    virtual std::string ClassName() const  { return "WeaveMaterial"; };
    mtl::Dense_Matrix GetShearedPlaneStress12( double ); //returns the PlaneStressMatrix based on the ShearAngle in rad
    mtl::Dense_Matrix GetShearedThreeDShell( double );
    TransverseIsotropicMaterial23 GetSingleLayerMaterial() { return SingleLayerMaterial; };

  private:
  double E_Ratio;
  TransverseIsotropicMaterial23 SingleLayerMaterial;
    // print-function
    virtual void print(std::ostream&) const;
  }; //of class WeaveMaterial



  /*! A class representing the orientation of the material 1 direction and the
      shear angle of a weave material
  */
  class MaterialOrientation
  {
        public:
        //constructors
        MaterialOrientation():AngleMaterial1Direction(0.),ShearAngle(0.){};
        MaterialOrientation(double MatOneDir, double ShearAngle):AngleMaterial1Direction(MatOneDir),ShearAngle(ShearAngle) {};

        bool operator==(const MaterialOrientation& a)
        {
                return (AngleMaterial1Direction == a.AngleMaterial1Direction && ShearAngle == a.ShearAngle);
        }
        bool operator!=(const MaterialOrientation& a)
        {
                return !(AngleMaterial1Direction == a.AngleMaterial1Direction && ShearAngle == a.ShearAngle);
        }


        //Set functions
        //-------------
        void SetAngleMaterial1Direction(double value) { AngleMaterial1Direction = value; };
        void SetShearAngle(double value) { ShearAngle = value; };

        //Get functions
        //-------------
        double GetAngleMaterial1Direction() const { return AngleMaterial1Direction; };
        double GetShearAngle() const { return ShearAngle; };

        private:
        //both angles are measured in rads !!!!
        double AngleMaterial1Direction; //the angle of the material to the direction given from node 0 to node 1
        double ShearAngle; //the angle between the material 1 and the material 2 direction (pi/2 for an undistorted material)

  };//of class MaterialOrientation

}       // end of namespace fe_base

#endif
