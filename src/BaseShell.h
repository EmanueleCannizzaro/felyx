//-----------------------------------------------------------------------------
// BaseShell.h
//
// begin     : Nov 11 2002
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


#ifndef BaseShell_h
#define BaseShell_h BaseShell_h

#include "NumIntElem.h"

namespace fe_base{

  class BaseShell:public NumIntElem{

  public:
    BaseShell()            : NumIntElem(), EleCoordSysPtr(NULL) {};
    BaseShell(unsigned n_) : NumIntElem(n_), EleCoordSysPtr(NULL)  {};

    virtual ~BaseShell() {};

    virtual BaseShell* Clone() const = 0;


    // CLASS FUNCTIONS
    //----------------

    //!Sets a Matrix of the Nodalcoordinate systems
    /*!The returned matrix has the dimensions (3*NodeCount,3)
      the first three lines are the rotational matrix at the
      first node...
     */
    void GetNodalCoordSys( Dense_Matrix ) const;

    /*!Returns a rotational matrix for the point defined.
      The matrix represents a coordinate system with localZ normal
      to the shell  and localX according to the EleCoordSys.
      The point must be given in natural coordinates in the element
     */
    void GetMaterialCoordSys( const Dense_Matrix intpoints_, const unsigned ip_, Dense_Matrix CS_ );

    /*!The following function sets a coordinate system at a point on the midplane of the shell
      specified by the row i_ of matrix evalpoint_(anynumber, 3) in natural coordinates. The
      normal zeta is errected using the gradients in natural xi and eta directions. The local x direction is
      expected to lie in the global x-z-plane. If the zeta-axis is within a certain accuracy in this
      plane it will be set to the global y-z-plane.
      The return type is a rotational matrix
    */
    void GetCoordMatrix( const Dense_Matrix, const unsigned, Dense_Matrix CS_ ) const;
    void GetCoordMatrix( const double, const double, const double, Dense_Matrix CS_ ) const;

    //! Returns the normal
    /*! Returns the normal  with length 1 in the point specified by the line ..
      of the parametermatrix the point should be specified in natural coordinates
    */
    Dense_Vector EvalShellNormal( const Dense_Matrix IntPoints, unsigned ip) const;

    //!Returns the global thickness of the shell
    virtual double GetShellThickness() const = 0;

    //! Returns the global coordinates of a point in natural coords
    void Natural2GlobalCoords( const Dense_Vector, Dense_Vector );

    //! Returns the shape functions in a specified point
    void GetShapefunctions( const Dense_Vector, Dense_Vector  );

    //! This function calculates the ESM of a single layer of isotropic or anisotropic material
    /*! The angle must be given in degrees */
    Dense_Matrix SingleLayerESM( Material*, const Dense_Matrix, const MaterialOrientation&, const double );
    Dense_Matrix SingleLayerEMM( Material*, const Dense_Matrix, const MaterialOrientation& );

    //!Returns the volume of the present layer
    double GetShellVolume( const Dense_Matrix TmpDirCos_, const Dense_Matrix IntPoints,
                                vector<Dense_Matrix>&, vector<double>&  );

    //!Returns the area of the Shell
    double EvalArea() const;

    //! Evaluates the Jacobimatrix for a Shell in 3 dimensional space
    /*! the function returns a 3x3- Jacobimatrix for a specific point given by the IPth
      row of the matrix evalpoints. The shape function are beeing evaluated inside this
      function.
    */
    void EvalJacobi( Dense_Matrix, const Dense_Matrix,
                     const Dense_Matrix, const Dense_Matrix, unsigned) const;

    //!Expands the 5 DOF ESM to 6 DOF in the global coordinate system
    void ExpandESM2globalCoords(Dense_Matrix, const Dense_Matrix, const Dense_Matrix,
                                const Material*, double );

    // The following functions are here and not virtual since it is assumed that the MaterialTransformations
    // are only set once in the midplane of the shell. That is why the SetMaterialTransformationVector is set
    // from within the CalcESM function of the derived elements and not inside the SingleLayerESM function.

    //! Returns the vector of material transformation matrices as defined in SetMaterialTransformationVector;
    //vector<Dense_Matrix> GetMaterialTransformationVector() const { return MaterialTransformationVector; };
    //! Returns the vector mapping transformations to intpoints as defined in SetMaterialTransformationVector;
    //vector<unsigned> GetMaterialTransformations2IntPts() const { return MaterialTransformations2IntPts; };
    //!Set the Material transformation vector and MatTransform2IntPts for the chosen integration scheme
    //void SetMaterialTransformationVector();


    // VIRTUAL FUNCTIONS
    //------------------

    /*! Set CoordSysPtr, only if element allows element coordinate systems, otherwise
      EleCoordSysPtr is set to NULL --> therefore a virtual function */
    virtual void SetEleCoordSysPtr( CoordSys* Ptr ) { EleCoordSysPtr = Ptr; };

    //! Convertes standard coordinates to global coordinates for shells
    /*! The Jacobian is transfered to the function. Its inversion is done within the function.
      With the jacobian a new 9x9 matrix is built with the inverse Jacobian as sub-matrices
    */
    virtual void stCoord2globalCoord(const Dense_Matrix , Dense_Matrix );

    virtual void EvalDerivedDeformations( Dense_Matrix, const Dense_Matrix , const Dense_Matrix,
                                          const Dense_Matrix, unsigned) = 0;

    //! Sets the B Matrix
    virtual void SetBMatrix( const Dense_Matrix, Dense_Matrix, BMatrixTypes);
    virtual void SetNMatrix( const Dense_Matrix, Dense_Matrix, BMatrixTypes, Dense_Matrix, double, unsigned);

    //! Returns the ElementCoordinateSystemPointer
    virtual CoordSys* GetEleCoordSysPtr() const { return EleCoordSysPtr; };

    //! Returns the rotational stiffness for different shells
    virtual unsigned GetLayerCount() const = 0;

    /*! Eval strains at integration points (argument: Integration Points)
        The resulting matrix looks as follows:
        ex ey ez gxy gyz gxz  in a row for each int point
    */
    virtual Dense_Matrix EvalStrainsAtIntPoints( const Dense_Matrix );
    
    /*! Eval stresses at integration points (arguments: Integration Points, Material, and Material orientation angle)
        The resulting matrix looks as follows:
        sx sy sz txy tyz txz  in a row for each int point
    */
    virtual Dense_Matrix EvalStressesAtIntPoints( const Dense_Matrix, Material*, const float_type);
    
    /*! The original EvalStresses function throwing a warning only
    */
    virtual void EvalStresses( string );
    
    /*! Eval stresses or strains at shell element corners. 
        The function takes a matrix of stresses or strains at int points and 
        extrapolates these values to the shell element corners. The order of the 
        stresses or strains in the resulting matrix is according to the "GetlNodeCoords()"-
        function of the respective element.
    */
    virtual Dense_Matrix Extrapolate2Corners( const Dense_Matrix );

    //! Eval volume of shell element (all layers)
    virtual double EvalVolume() const;
    
    /*! Eval strains at int points in element coordinate system directions.
        The result matrix looks as folows:
        ex* ey* ez* gxy* gyz* gxz* in a row for each int point.
    */
    virtual Dense_Matrix evalElementDirectionStrains( const Dense_Matrix IntPoints );
    
    /*! Eval strains at int points in main material directions 
        The result matrix looks as follows:
        e1 e2 e3 g12 g23 g13  in a row for each int point.
    */
    virtual Dense_Matrix evalMaterialDirectionStrains( const Dense_Matrix IntPoints, const float_type );
    
    /*! Eval stresses at int points in element coordinate system directions.
        The result matrix looks as folows:
        sx* sy* sz* txy* tyz* txz* in a row for each int point.
        (arguments: Integration Points, Material, and Material orientation angle)
    */
    virtual Dense_Matrix evalElementDirectionStresses( const Dense_Matrix IntPoints, Material*, const float_type );
    
    /*! Eval stresses at in points in main material directions
      The result matrix looks as follows:
      s1 s2 s3 t12 t23 t13
      (arguments: Integration Points, Material, and Material orientation angle)
    */
    virtual Dense_Matrix evalMaterialDirectionStresses( const Dense_Matrix IntPoints, Material*, const float_type );
    
    /////
    // Failure criteria
    /////
    
    /*! Eval Maximum Stress Criteria at each int point
        The resulting matrix holds the following entries for each integration point
        s1/Xt; |s1|/Xc; s2/Yt; |s2|/Yc; |t12|/S
    */
    virtual Dense_Matrix evalMaxStressCriteria (const Dense_Matrix IntPoints, const float_type, Material* );
    
    /*! Eval Tsai-Hill Criteria at each int point
        The resulting matrix holds the Tsai-Hill value for each integration point.
        Tsai-Hill <= 1;
    */
    virtual Dense_Matrix evalTsaiHillCriteria (const Dense_Matrix IntPoints, const float_type, Material* );
    
    /*! Eval Tsai-Wu Criteria at each int point
        The resulting matrix holds the Tsai-Wu value for each integration point
        Tsai-Wu <= 1
    */
    virtual Dense_Matrix evalTsaiWuCriteria (const Dense_Matrix IntPoints, const float_type, Material* );
    
    /*! Eval Hashin Criteria at each integration point
        The resulting matrix holds all four subcriteria values at each integration points
        fibre tension failure; fibre compression failure; inter-fibre tension failure; inter-fibre compression failure;
        all values must be below 1  
    */
    virtual Dense_Matrix evalHashinCriteria (const Dense_Matrix IntPoints, const float_type, Material* );
    

  public:
    // Data members
    // -------------------
    Dense_Matrix IntPointStress;

    virtual Dense_Matrix GetlNodeCoords() = 0;

  protected:

    //! Evaluate strains of a shell for specific integration points
    Dense_Matrix evalStrains(const Dense_Matrix IntPoints );    
    vector<Dense_Matrix>     MaterialTransformationVector;
    vector<unsigned>         MaterialTransformations2IntPts;
    

  private:
    //DATA MEMBERS
    //------------
    //! Pointer to Element Coordinate system
    CoordSys* EleCoordSysPtr;


  };

    /*!This function aligns a localX direction with respect to
    a given coordinate system represented in the Matrix CS_.
    CS_ is the rotational matrix. The localZ direction is passed
    to the function as normal_. The parameter toler_ gives a tolerance
    to swap from CS_-x- to CS_-y-direction. The tolerance must be given in
    radiants. The resulting CoordSys is stored in CS_
  */
  void AlignXY( Dense_Matrix CS_, const Dense_Vector normal_, const double toler_ );

  /*!Returns the material transformation matrix as defined by cook..
     The argument is the rotation matrix
     Since not everyone is able to find cook (due to the excellent literature reference here),
     some more commenting may be useful:
     
     generally it holds:  sigma'_ij = a_ik * a_jl * sigma_kl (lovely tensors), where a_xy are the entries
     of the rotational matrix. In matrix notation something like: sigma_123 = T * sigma_xyz.
     Since the strains do not fulfill the general requirments of tensors when
     calculating a matrix T for a_ik * a_jl, the reuter matrix needs to be applied.
     epsilon_123 = R * T^(-1) * R^(-1) * epsilon_xyz = MatTransMatrix * epsilon_xyz
  */
  Dense_Matrix GetMaterialTransformationMatrix( Dense_Matrix CoordSys );
  
  /*! Returns the stress transformation matrix required to transform stresses
      from one coordinate system to another. it is based on the above defined 
      material transformation matrix by cook, whereas the reuter manipulations are 
      made applied in an inverse way to finally get the matrix T.
  */
  Dense_Matrix GetStressTransformationMatrix( Dense_Matrix CoordSys );

} //end of namespace

#endif
