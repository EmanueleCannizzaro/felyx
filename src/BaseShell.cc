//-----------------------------------------------------------------------------
// BaseShell.cc
//
// begin     : Nov 20 2001
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


#include "BaseShell.h"

using namespace fe_base;

//Returns the ESM of a single layer in global coords
Dense_Matrix BaseShell::SingleLayerESM( Material* MatPtr_, const Dense_Matrix IntPoints_,
                                        const MaterialOrientation& LayerAngle_, const double LayerThickness_){

  unsigned NodeCount = GetNodeCount(), nSize = NodeCount*(GetDofSet().count()-1);

  double PI = 4*atan(1.);

  Dense_Matrix
    SingleLayerESM(nSize+NodeCount, nSize+NodeCount),
    J(3,3),
    B(6,nSize),
    K(nSize,nSize),
    Result_at_current_intpoint(nSize,nSize),
    TransformedMaterial(6,6),
    MatTransformationMatrix(6,6),
    TmpCoordSys(3,3),
    DerivedDeformations(9, nSize),
    PlaneStressMatrix(3,3),
    RotatedPlaneStressMatrix(3,3),
    D(6,6);

  vector<Dense_Matrix> J_vec(IntPoints_.nrows());
  vector<double> detJ_vec(IntPoints_.nrows());

  //Compute the nodal coordinate systems
  Dense_Matrix NodalCoordSystems(3*NodeCount,3);
  GetNodalCoordSys( NodalCoordSystems );

  //Get the material data (the stress-strain-relation)
  //The diagonal elements 4 and 5 have to be scaled by a parameter
  //depending on the thickness and the area of the element
  //see next comment
  //mtl::copy((MatPtr_->Get(Material::ThreeDShell)), D);


  if ( MatPtr_ -> ClassName() == "WeaveMaterial" && (fabs(LayerAngle_.GetShearAngle()-PI/2.)) > 0.01 )// if it is a weave material and it is distorted
  {
     mtl::copy(((WeaveMaterial*)MatPtr_)->GetShearedThreeDShell(LayerAngle_.GetShearAngle()),D);
  }
  else
  {
    mtl::copy((MatPtr_->Get(Material::ThreeDShell)), D);
  }



  //To account for the thickness-direction variation of transverse shear-strain which is nearly parabolic
  //and not constant like assumed by the Mindlin theory. The different scaling is taken from ANSYS.
  //The cost of the function calculating the area could be intense... to be checked
  //Here one could calculate the shell area, instead of the volume...
  //The vector of Jacobi Matrices is set as well inside...
  double Volume = GetShellVolume( NodalCoordSystems, IntPoints_, J_vec, detJ_vec );
//  double thickness = GetShellThickness();

  double f = 1.0 + 0.2*Volume/25.0/LayerThickness_/LayerThickness_/LayerThickness_;
  if ( f < 1.2 ) f = 1.2;

  D(4,4) /= f;
  D(5,5) /= f;

  //Rotate the Material about the zeta-Axis
if ( (MatPtr_ -> ClassName() == "TransverseIsotropicMaterial23" || MatPtr_ -> ClassName() == "OrthotropicMaterial" ||
      MatPtr_->ClassName() == "WeaveMaterial" ) &&
        fabs(LayerAngle_.GetAngleMaterial1Direction()) > 0.01 ){
    Dense_Matrix rot(3,3), copyD(6,6), TransformationMatrix(6,6);
    mtl::copy(CoordSys::Vec2Mat(0.0,0.0,LayerAngle_.GetAngleMaterial1Direction()), rot);
    mtl::copy(D,copyD);
    mtl::set_value(TransformationMatrix,0.);
    TransformationMatrix = GetMaterialTransformationMatrix( rot );
    mtl::set_value(D,0.);
    mult_AT_B_A(TransformationMatrix,copyD,D); //Materialproperties rotated to the appropriate angle
  }

/*
 ///SetMaterialTransformationVector();
  std::vector<Dense_Matrix> TransformedMaterialVector(GetMaterialTransformationVector().size());
  for (unsigned k = 0 ; k < TransformedMaterialVector.size() ; ++k )
  {
        //apply the transformation
        TransformedMaterialVector[k] = mtl::Dense_Matrix(6,6);
        mtl::set_value(TransformedMaterialVector[k],0.0);
        mult_AT_B_A(GetMaterialTransformationVector()[k], D , TransformedMaterialVector[k]);
  }
*/

 //Initialization of all variables
  double fdummy;

  for (unsigned i = 0 ; i < IntPoints_.nrows(); ++i){


    //Derived shape functions N with respect to  xi and eta coordinates
    //computed at the point represented by the row i of the matrix
    //IntPoints. They are stored in shapefunc, which will be used in different functions
    //The third column represents the shapefunctions (not derived) which will be used
    //for shell elements
    //Compute the derived deformations with respect to natural
    //coordinates
    //uses the matrix shapefunc
    EvalDerivedDeformations(DerivedDeformations, NodalCoordSystems, EvalDerivedShapeFunc( IntPoints_[i] ), IntPoints_, i);
    //Transforme the shapefunctions to global x, y, z
    //coordinates

    stCoord2globalCoord(J_vec[i], DerivedDeformations);
    //The integration
    SetBMatrix(DerivedDeformations, B, ShellIntegration);
    //IntPoints(i,3) is the weight of the ith point.
    //It is included to maintain the ability to change the
    //degree of integration easely

    fdummy = detJ_vec[i]*IntPoints_(i,IntPoints_.ncols()-1);

    //Materialtransformation at the present gausspoint
    GetMaterialCoordSys(IntPoints_, i,TmpCoordSys);

    MatTransformationMatrix = fe_base::GetMaterialTransformationMatrix( TmpCoordSys );

//    AlignMaterial2CoordSys(TmpCoordSys,D,TransformedMaterial);
    mult_AT_B_A(MatTransformationMatrix,D, TransformedMaterial);

    scale(TransformedMaterial, fdummy);
/*
    //Materialtransformation at the present gausspoint
    mtl::copy(TransformedMaterialVector[GetMaterialTransformations2IntPts()[i]],TransformedMaterial);

    fdummy = detJ_vec[i]*IntPoints_(i,IntPoints_.ncols()-1);
    scale(TransformedMaterial, fdummy);
*/
    //add the stiffness of the current gausspoint to the stiffnessmatrix K
    mult_AT_B_A_add(B,TransformedMaterial,K);


  } //of the integration loop
  ExpandESM2globalCoords(SingleLayerESM, K,NodalCoordSystems, MatPtr_, Volume);

  return SingleLayerESM;
}

Dense_Matrix BaseShell::SingleLayerEMM( Material* MatPtr_, const Dense_Matrix IntPoints_,
                                        const MaterialOrientation& LayerAngle_ ){

  unsigned NodeCount = GetNodeCount(), nSize = NodeCount*(GetDofSet().count()-3);

  Dense_Matrix
    SingleLayerESM(nSize+3*NodeCount, nSize+3*NodeCount),
    J(3,3),
    N(3,nSize),
    M(nSize,nSize),
    NodalCoordSystems(3*NodeCount,3);

  GetNodalCoordSys( NodalCoordSystems );

  double thickness = GetShellThickness();

  double detJ, fdummy;

  for (unsigned i = 0 ; i < IntPoints_.nrows(); ++i){
    Dense_Matrix derivedshapefunc = EvalDerivedShapeFunc(IntPoints_[i]);

    EvalJacobi(J, derivedshapefunc, NodalCoordSystems, IntPoints_, i);

    SetNMatrix(derivedshapefunc, N, ShellIntegration, IntPoints_, thickness, i);

    detJ = determinant(J);

    fdummy = detJ*IntPoints_(i,IntPoints_.ncols()-1)*MatPtr_->Get("rho");

    mult(trans(N),scaled(N,fdummy),M);

  } //of the integration loop

  unsigned u, v, i, j;

  for ( j = 0 ; j < NodeCount ; ++j ){
    for ( i = 0 ; i < NodeCount ; ++i ){
      for ( u = 0 ; u < 3 ; ++u ){
        for ( v = 0 ; v < 3 ; ++v ){
          SingleLayerESM(6*j+u,6*i+v) = M(3*j+u,3*i+v);
        }
      }
    }
  }

  return SingleLayerESM;
}

void BaseShell::SetBMatrix( const Dense_Matrix D, Dense_Matrix B, BMatrixTypes type ){

  mtl::set_value(B,0.0);

  switch (type)
    {
    case ShellIntegration :
      {
        if (GetDofSet().IsElement3DShell()){

          Dense_Matrix H(6,9);
          mtl::set_value(H,0.0);
          H(0,0) = 1;
          H(1,4) = 1;
          H(2,8) = 1;
          H(3,1) = 1;     H(3,3) = 1;
          H(4,5) = 1;     H(4,7) = 1;
          H(5,2) = 1;     H(5,6) = 1;

          //B_ = H * DerivedDeformations_;
          mult(H, D ,B);
        }

        else
          {
            cerr << endl << "####         Error in NumIntElem::ShellIntegration           ####"
                 << endl << "#### The Element calling ShellIntegration is not a Shell! ####"
                 << endl << "#### It is a " << GetName() <<  " number " << GetId() << " ##### " << endl;
          }
        break;
      }
    default:
      {
        cerr << endl << "####         Error in BaseShell::SetBMatrix          ####"
             << endl << "#### The selected integration scheme is not implemented! ####" << endl;
        break;
      }
    }

}


void BaseShell::SetNMatrix( const Dense_Matrix S, Dense_Matrix N, BMatrixTypes type, Dense_Matrix intpoints, double t, unsigned i ){

  mtl::set_value(N,0.0);

  switch (type)
    {
    case ShellIntegration :
      {
        if (GetDofSet().IsElement3DShell()){
          for (unsigned k = 0 ; k < GetNodeCount() ; k++)
            {
              N(0,(k+1)*3-3) = S(k,2);
              N(1,(k+1)*3-2) = S(k,2);
              N(2,(k+1)*3-1) = S(k,2);
              //(1,(k+1)*5-2) = -S(k,2)*intpoints(i,zeta)*t/2.0;
              //0,(k+1)*5-1) = S(k,2)*intpoints(i,zeta)*t/2.0;
            }
        }
        else
          {
            cerr << endl << "####         Error in NumIntElem::ShellIntegration           ####"
                 << endl << "#### The Element calling ShellIntegration is not a Shell! ####"
                 << endl << "#### It is a " << GetName() <<  " number " << GetId() << " ##### " << endl;
          }
        break;
      }
    default:
      {
        cerr << endl << "####         Error in BaseShell::SetBMatrix          ####"
             << endl << "#### The selected integration scheme is not implemented! ####" << endl;
        break;
      }
    }

}

void BaseShell::EvalJacobi( Dense_Matrix J_, const Dense_Matrix DerivedShapeFunc_,
                            const Dense_Matrix TmpDirCos, const Dense_Matrix EvalPoints_,
                            unsigned IP) const{

  //initialize the matrix to zero
  mtl::set_value(J_,0.0);

  double t = GetShellThickness(), zetaCoord;
  unsigned zeta = EvalPoints_.ncols()-2;

  for (unsigned k = 0 ; k < GetNodeCount() ; ++k){
    //areacoords have one dimension more,
    //this is the  natural zetacoordinate
    zetaCoord = EvalPoints_(IP,zeta);

    J_(0,0) += DerivedShapeFunc_(k,0)*(NodeVec[k]->Cx + zetaCoord*t*TmpDirCos(2+k*3,0)*0.5);
    J_(1,0) += DerivedShapeFunc_(k,1)*(NodeVec[k]->Cx + zetaCoord*t*TmpDirCos(2+k*3,0)*0.5);
    J_(2,0) += DerivedShapeFunc_(k,2)*t*TmpDirCos(2+k*3,0)*0.5;

    J_(0,1) += DerivedShapeFunc_(k,0)*(NodeVec[k]->Cy + zetaCoord*t*TmpDirCos(2+k*3,1)*0.5);
    J_(1,1) += DerivedShapeFunc_(k,1)*(NodeVec[k]->Cy + zetaCoord*t*TmpDirCos(2+k*3,1)*0.5);
    J_(2,1) += DerivedShapeFunc_(k,2)*t*TmpDirCos(2+k*3,1)*0.5;

    J_(0,2) += DerivedShapeFunc_(k,0)*(NodeVec[k]->Cz + zetaCoord*t*TmpDirCos(2+k*3,2)*0.5);
    J_(1,2) += DerivedShapeFunc_(k,1)*(NodeVec[k]->Cz + zetaCoord*t*TmpDirCos(2+k*3,2)*0.5);
    J_(2,2) += DerivedShapeFunc_(k,2)*t*TmpDirCos(2+k*3,2)*0.5;
 }

}


void BaseShell::stCoord2globalCoord(const Dense_Matrix Jacobi_, Dense_Matrix DerivedDeformations_){
  if (!GetDofSet().IsElement3DShell()){
    cerr << endl << "#####################################################";
    cerr << endl << "ERROR: in fe_base::Element::stCoord2globalCoordShell: ";
    cerr << endl << "ERROR: The Element is not a Shell.";
    cerr << endl << "#####################################################";
  }
  else{
    unsigned size = Jacobi_.nrows();
    Dense_Matrix invJacobi(size, size);
    inversion(Jacobi_, invJacobi);

    Dense_Matrix tripleJacobi(9,9), dummyMatrix(DerivedDeformations_.nrows(), DerivedDeformations_.ncols() );
    copy(DerivedDeformations_,dummyMatrix);
    for ( unsigned i = 0 ; i < 3 ; ++i ){
      for ( unsigned k = 0 ; k < 3 ; ++k ){
        tripleJacobi(i,k)     = invJacobi(i,k);
        tripleJacobi(3+i,3+k) = invJacobi(i,k);
        tripleJacobi(6+i,6+k) = invJacobi(i,k);
      }
    }

    mtl::set_value(DerivedDeformations_,0.0);
    mult(tripleJacobi,dummyMatrix,DerivedDeformations_);

  }
}

void BaseShell::GetNodalCoordSys( Dense_Matrix NodalCoordSys_ ) const{

  unsigned nc = GetNodeCount();

  // Check if NodalCoordSys_ has appropriate size
  MTL_ASSERT( NodalCoordSys_.nrows() == 3*nc && NodalCoordSys_.ncols() == 3, "felyx::BaseShell::GetNodalCoordSys()");

  // Eval the matrix of node point locations
  Dense_Matrix  NPL(nc,4);
  mtl::set_value(NPL,0.0);
  if ( nc == 8 ) {  //if it is a quadrilateral
    NPL(0,0) = -1;   NPL(0,1) = -1;    NPL(0,2) = 0;
    NPL(1,0) =  1;   NPL(1,1) = -1;    NPL(1,2) = 0;
    NPL(2,0) =  1;   NPL(2,1) =  1;    NPL(2,2) = 0;
    NPL(3,0) = -1;   NPL(3,1) =  1;    NPL(3,2) = 0;
    NPL(4,0) =  0;   NPL(4,1) = -1;    NPL(4,2) = 0;
    NPL(5,0) =  1;   NPL(5,1) =  0;    NPL(5,2) = 0;
    NPL(6,0) =  0;   NPL(6,1) =  1;    NPL(6,2) = 0;
    NPL(7,0) = -1;   NPL(7,1) =  0;    NPL(7,2) = 0;
  }
  else if ( nc == 6 ) { // if it is a triangle
    //The coordinates are given in area coordinates
    //Zienkiewicz p. 180 ff.
    NPL(0,0) = 1.0;  NPL(0,1) = 0.0; NPL(0,2) = 0.0; NPL(0,3) = 0.0;
    NPL(1,0) = 0.0;  NPL(1,1) = 1.0; NPL(1,2) = 0.0; NPL(1,3) = 0.0;
    NPL(2,0) = 0.0;  NPL(2,1) = 1.0; NPL(2,2) = 1.0; NPL(2,3) = 0.0;
    NPL(3,0) = 0.5;  NPL(3,1) = 0.5; NPL(3,2) = 0.0; NPL(3,3) = 0.0;
    NPL(4,0) = 0.0;  NPL(4,1) = 0.5; NPL(4,2) = 0.5; NPL(4,3) = 0.0;
    NPL(5,0) = 0.5;  NPL(5,1) = 0.0; NPL(5,2) = 0.5; NPL(5,3) = 0.0;
  }

  mtl::set_value(NodalCoordSys_,0.0);

  Dense_Matrix TmpCoordSys(3,3);
  for ( unsigned i = 0 ; i < nc ; ++i ){

    GetCoordMatrix( NPL, i, TmpCoordSys );

    //copy the row vectors
    for ( unsigned k = 0 ; k < 3 ; ++k )
      copy(TmpCoordSys[k],NodalCoordSys_[3*i+k]);
  }
}

void BaseShell::GetCoordMatrix( const Dense_Matrix pointlocation_, const unsigned ip_, Dense_Matrix CS_ ) const{

  // Set vals
  mtl::set_value(CS_,0.0);
  mtl::set_diagonal(CS_,1.0);

  Dense_Vector normal = EvalShellNormal( pointlocation_, ip_ );
  //mtl::copy(GetShellNormal(pointlocation_, ip_),normal);
  //GetShellNormal(pointlocation_, ip_, normal);

  AlignXY(CS_, normal, 0.1);
}

//Check Shell93::GetCoordMatrix(Dense_Matrix, unsigned) for more info
void BaseShell::GetCoordMatrix( const double xi_, const double eta_, const double zeta_, Dense_Matrix CS_ ) const {
  Dense_Matrix Pointlocation(1,3);
  Pointlocation(0,0) = xi_;
  Pointlocation(0,1) = eta_;
  Pointlocation(0,2) = zeta_;

  GetCoordMatrix( Pointlocation, 0, CS_ );
}

Dense_Vector BaseShell::EvalShellNormal( const Dense_Matrix IntPoints, unsigned ip) const{

  // Eval derived shape functions
  Dense_Matrix DerivedShapeFuncs = EvalDerivedShapeFunc( IntPoints[ip] );

  // Build the gradients in xi and eta directions
  // one could use the jacobimatrix for this step but in this
  // shell implementation the shapefunctions derived with respect
  // to zeta coordinates are not stored...
  Dense_Vector localX(3,0.0), localY(3,0.0), normal(3,0.0);
  for (  unsigned i = 0 ; i < GetNodeCount() ; ++i ){
    localX[0] += DerivedShapeFuncs(i,0)*NodeVec[i]->Cx;
    localX[1] += DerivedShapeFuncs(i,0)*NodeVec[i]->Cy;
    localX[2] += DerivedShapeFuncs(i,0)*NodeVec[i]->Cz;

    localY[0] += DerivedShapeFuncs(i,1)*NodeVec[i]->Cx;
    localY[1] += DerivedShapeFuncs(i,1)*NodeVec[i]->Cy;
    localY[2] += DerivedShapeFuncs(i,1)*NodeVec[i]->Cz;
  }

  //the jacobian in the z-direction is not defined in a shell
  cross_prod_3d(localX, localY, normal);

  float_type norm2 = 1.0/ mtl::two_norm( normal );
  mtl::scale(normal,norm2);

  return normal;
}


//This function returns the global coords of an arbitrary
//point p_ in natural coordinates. p_ can be given in 3 cartesian or
//4 area coords (actually 3 areacoords and 1 cartesian zeta coord.
//The size of p_ is used to determine the coordinate system!!!!
void BaseShell::Natural2GlobalCoords( const Dense_Vector p_, Dense_Vector global){

  unsigned nc = GetNodeCount();

  Dense_Vector normal(3,0.0), shapef(nc,0.0);

  //Get the vector of shape functions
  GetShapefunctions(p_,shapef);

  Dense_Matrix NCS(3*nc,3);
  GetNodalCoordSys( NCS );

  //The xi and eta coords are already included in the shapefuncs
  double zeta = 0.0, l,m,n, t = 0.8;

    //  = GetShellThickness();

  if (p_.size() == 3) zeta = p_[2];
  if (p_.size() == 4) zeta = p_[3];

  for (unsigned i = 0 ; i < nc ; ++i){
    mtl::copy(NCS[2+i*3], normal);
    l = normal[0]; m = normal[1]; n = normal[2];

    global[0] += shapef[i]*(NodeVec[i]->Cx + zeta*t*l/2.0);
    global[1] += shapef[i]*(NodeVec[i]->Cy + zeta*t*m/2.0);
    global[2] += shapef[i]*(NodeVec[i]->Cz + zeta*t*n/2.0);
  }
}

void BaseShell::GetMaterialCoordSys( const Dense_Matrix pointlocation_, const unsigned ip_, Dense_Matrix CS_ ) {

  // Reset vals in CS_
  mtl::set_value(CS_,0.0);

  unsigned size = pointlocation_.ncols() - 1; //only the coordinates are important, not the weight
  //initialize the precision
  double tolerance = 0.0;

  //Set the normal
  Dense_Vector normal=EvalShellNormal( pointlocation_, ip_ );

  //if no special coordinate system is set the local x  is aligned parallel
  //to the connection from node 0 to node 1
  if ( EleCoordSysPtr == NULL ) {
    //they are used as projection lines, they dont have to be perpendicular
    mtl::add(NodeVec[1]->Get() , scaled( NodeVec[0]->Get(), -1.0) , CS_[0] );
    mtl::add(NodeVec[3]->Get() , scaled( NodeVec[0]->Get(), -1.0) , CS_[1] );
  }

  //if a special coordinate system is set it depends on the type
  //of EleCoordSys. The function GetLocalCS(..) makes the distinction
  else if ( EleCoordSysPtr->GetCSType() == CoordSys::cartesian ||   EleCoordSysPtr->GetCSType() == CoordSys::cylindrical ){

    //tolerance is 45 deg this is as ANSYS does it ....
    tolerance =  1.0/sqrt(2.0);

    //get the point (in natural coords!)
    Dense_Vector localX(size,0.0);
    mtl::copy(pointlocation_[ip_](0,3), localX);

    //transforme it to global coords
    Dense_Vector globalX(3,0.0);
    Natural2GlobalCoords(localX, globalX);

    //get the CS at the specified point
    mtl::copy(EleCoordSysPtr->GetLocalCS( globalX ), CS_ );
  }
  else {
    //std::cout << "WARNING: Element coordinate system type is not implemented: " << EleCoordSysPtr->GetType() << std::endl;
    FELYX_RUNTIME_THROW("BaseShell::GetMaterialCoordSys(): unknown coord sys type");
  }

  AlignXY(CS_, normal, tolerance);
}


//Takes the location point_ in natural coords and returns a
//Vector containing the shapefunctions for that point the point
//can be given in 3 cartesian coords or in 4 area coords
void BaseShell::GetShapefunctions( const Dense_Vector p_, Dense_Vector shapef ) {

  FELYX_LOGIC_ASSERT( abs(p_[0]) <= 1.0 && abs(p_[1]) <= 1.0 && abs(p_[2]) <= 1.0, "BaseShell::GetShapefunctions()");
  FELYX_LOGIC_ASSERT( p_.size() == 3 || p_.size() ==4 , "BaseShell::GetShapefunctions()" );

  if (p_.size() == 3){
    double xi = p_[0], eta = p_[1]; //zeta is not used for this shell

    shapef[0] = 0.25*(1-xi)*(1-eta)*(-xi-eta-1);
    shapef[1] = 0.25*(1+xi)*(1-eta)*(xi-eta-1);
    shapef[2] = 0.25*(1+xi)*(1+eta)*(xi+eta-1);
    shapef[3] = 0.25*(1-xi)*(1+eta)*(-xi+eta-1);
    shapef[4] = 0.5*(1-xi*xi)*(1-eta);
    shapef[5] = 0.5*(1+xi)*(1-eta*eta);
    shapef[6] = 0.5*(1-xi*xi)*(1+eta);
    shapef[7] = 0.5*(1-xi)*(1-eta*eta);
  }
  else{
    double L1 =p_[0], L2 =p_[1], L3 =p_[2];

    shapef[0] = L1 * (2*L1 - 1.0);
    shapef[1] = L2 * (2*L2 - 1.0);
    shapef[2] = L3 * (2*L3 - 1.0);
    shapef[3] = 4*L1*L2;
    shapef[4] = 4*L2*L3;
    shapef[5] = 4*L3*L1;
  }

}



double BaseShell::GetShellVolume( const Dense_Matrix TmpDirCos_, const Dense_Matrix IntPoints, vector<Dense_Matrix>& J_vec, vector<double>& detJ_vec )  {

  Dense_Matrix J(3,3);// IntPoints = GetIntPoints();

  double volume = 0;

  for (unsigned i = 0 ; i < IntPoints.nrows() ; ++i ){

    EvalJacobi(J, EvalDerivedShapeFunc( IntPoints[i] ), TmpDirCos_, IntPoints, i);
    J_vec[i].resize(3,3);
    mtl::copy(J,J_vec[i]);

    detJ_vec[i] = determinant(J_vec[i]);

    volume += (detJ_vec[i]*IntPoints(i,IntPoints.ncols()-1));
  }

  return volume;
}

double BaseShell::EvalVolume() const{

  //Compute the nodal coordinate systems
  Dense_Matrix NodalCoordSystems(3*GetNodeCount(),3);
  GetNodalCoordSys( NodalCoordSystems );

  Dense_Matrix J(3,3);
  Dense_Matrix IntPoints = GetIntPoints();
  float_type volume = 0.0, detJ=0.0;

  for (unsigned i = 0 ; i < IntPoints.nrows() ; ++i ){
    EvalJacobi(J, EvalDerivedShapeFunc( IntPoints[i] ), NodalCoordSystems, IntPoints, i);

    detJ = determinant(J);
    volume += (detJ*IntPoints(i,IntPoints.ncols()-1));
  }
  return volume;
  //return volume/GetShellThickness();
}

double BaseShell::EvalArea() const{
  return EvalVolume()/GetShellThickness();
}


//Expand the ESM to meet global coords dimensions
//from 5 to 6 DOF per node
//The drilling dof is given a small stiffness value (check the function)
//Returned is a Elementstiffnessmatrix with all global DOF active for all nodes
void BaseShell::ExpandESM2globalCoords( Dense_Matrix K, const Dense_Matrix K_,
                                        const Dense_Matrix NodalCoordSys_,
                                        const Material* MatPtr, double Volume){

  unsigned u, v, l, m, nc = GetNodeCount(), size = nc*GetDofSet().count();

  Dense_Matrix T(size,size),final_K(size,size), dummy(3,3), dummyT(3,3), dummy2(3,3), dummyT2(3,3), dummyTT(3,3);

  mtl::set_value(T,0.0);
  mtl::set_value(final_K,0.0);

  //The next two factors are derived from explicit data within ANSYS
  //The ansys manual sais, a small value is added to the drilling DOF
  //but it is not mentioned what the value looks like
  //The following factors give good results for isotropic material
  //They will have to be checked again for orthotropic material
  //could be, that these factors are dependent on shape functions...

  //For triangular shells the factor for node 3 (Ansys counting) is wrong
  //within ANSYS. It is different from all the others!!!!

  double thewonderfactor = 0.0;

  if( MatPtr->ClassName() == "IsotropicMaterial")
    thewonderfactor = Volume * MatPtr->Get("E") * 1e-5 / (static_cast<double>(GetLayerCount()));

  else if( MatPtr->ClassName() == "TransverseIsotropicMaterial23"  || MatPtr->ClassName() == "WeaveMaterial" )

    thewonderfactor = (Volume * ((2*(MatPtr->Get("E1")) + (MatPtr->Get("E2")))/3 * 1e-5)) / (static_cast<double>(GetLayerCount()));

  else
    cerr << endl << " The requested material is not implemented " << endl;

  double thesmallerwonderfactor = -thewonderfactor / 7.0;

  // upper left block
  for ( u = 0 ; u < nc  ; ++u ){
    for ( v = 0 ; v <= u ; ++v ){
      for ( l = 0 ; l < 3 ; ++l ){
        for ( m = 0 ; m < 3 ; ++m ){
          K(6*u+l,6*v+m) = K_(5*u+l,5*v+m);
        }
      }
    }
  }

  // upper right block
  for ( u = 0 ; u < nc  ; ++u ){
    for ( l = 0 ; l < 3 ; ++l ){
      for ( m = 0 ; m < 3 ; ++m ){
        dummyT(l,m) = NodalCoordSys_(u*3+l,m);
      }
    }
    for ( v = u ; v < nc ; ++v ){
      for ( l = 0 ; l < 3 ; ++l ){
        for ( m = 0 ; m < 2 ; ++m ){
          dummy(l,m) = K_(5*v+l,5*u+m+3);
        }
      }
      mtl::set(dummy2,0.0);
      mult(dummy,dummyT,dummy2);
      for ( l = 0 ; l < 3 ; ++l ){
        for ( m = 0 ; m < 3 ; ++m ){
          K(6*v+l,6*u+m+3) = dummy2(l,m);
        }
      }
    }
  }

  // lower left block
  for ( u = 0 ; u < nc  ; ++u ){
    for ( l = 0 ; l < 3 ; ++l ){
      for ( m = 0 ; m < 3 ; ++m ){
        dummyT(l,m) = NodalCoordSys_(u*3+l,m);
      }
    }
    mtl::set(dummyTT,0.0);
    transpose(dummyT,dummyTT);
    mtl::set(dummy,0.0);
    for ( v = 0 ; v <= u ; ++v ){
      for ( l = 0 ; l < 2 ; ++l ){
        for ( m = 0 ; m < 3 ; ++m ){
          dummy(l,m) = K_(5*u+l+3,5*v+m);
        }
      }
      mtl::set(dummy2,0.0);
      mult(dummyTT,dummy,dummy2);
      for ( l = 0 ; l < 3 ; ++l ){
        for ( m = 0 ; m < 3 ; ++m ){
          K(6*u+l+3,6*v+m) = dummy2(l,m);
        }
      }
    }
  }

  // lower right block
  for ( u = 0 ; u < nc  ; ++u ){
    for ( l = 0 ; l < 3 ; ++l ){
      for ( m = 0 ; m < 3 ; ++m ){
        dummyT(l,m) = NodalCoordSys_(u*3+l,m);
      }
    }
    mtl::set(dummyTT,0.0);
    transpose(dummyT,dummyTT);
    for ( v = 0 ; v <= u ; ++v ){
      for ( l = 0 ; l < 3 ; ++l ){
        for ( m = 0 ; m < 3 ; ++m ){
          dummyT2(l,m) = NodalCoordSys_(v*3+l,m);
        }
      }
      dummy(0,0) = K_(5*u+3,5*v+3);
      dummy(0,1) = K_(5*u+3,5*v+4);
      dummy(1,0) = K_(5*u+4,5*v+3);
      dummy(1,1) = K_(5*u+4,5*v+4);
      if(v==u) dummy(2,2)=thewonderfactor;
      else dummy(2,2)=thesmallerwonderfactor;
      dummy(0,2)=0.0; dummy(1,2)=0.0;
      dummy(2,0)=0.0; dummy(2,1)=0.0;
      mtl::set(dummy2,0.0);
      mult(dummyTT,dummy,dummy2);
      mtl::set(dummy,0.0);
      mult(dummy2,dummyT2,dummy);
      for ( l = 0 ; l < 3 ; ++l ){
        for ( m = 0 ; m < 3 ; ++m ){
          K(6*u+l+3,6*v+m+3) = dummy(l,m);
        }
      }
    }
  }

  // fill symmetricly
  for ( u = 0 ; u < size  ; ++u ){
    for ( v = u ; v < size  ; ++v ){
          K(u,v)=K(v,u);
        }
      }

} //of ExpandESM2globalCoords

///
// Evaluate strains of a shell for specific integration points
///
Dense_Matrix BaseShell::evalStrains( const Dense_Matrix IntPoints_ ){

  unsigned ndofNode = GetDofSet().count()-1; // all global dofs are required
  unsigned ndofElem = GetNodeCount() * ndofNode; 

  // Resize and initalize matrix to store strains
  Dense_Matrix Strains( IntPoints_.nrows(), 6 );

  // Write deformations of all nodes into single vector, i.e. the global deformations at nodes
  // [x1 y1 z1 a1 b1 c1 x2 y2 z2 a2 b2 c2 ... xn yn zn an bn cn] 
  Dense_Vector u(ndofElem);
  for (unsigned node = 0 ; node< GetNodeCount() ; ++node){
   for (unsigned dof = 0 ; dof < ndofNode ; ++dof){
      u[ node*ndofNode + dof ] =  NodeVec[node]->GetDeformations()[dof];
    }
  }

  // Eval nodal coordinate systems required for evaluating the local deformations
  Dense_Matrix NodalCoordSystems(3*GetNodeCount(),3);
  GetNodalCoordSys( NodalCoordSystems );

  // Evaluate strains at all integration points
  Dense_Matrix DerivedDeformations(9, ndofElem), Jacobi(3,3), B(6,ndofElem);;
  for (unsigned i = 0 ; i < IntPoints_.nrows(); ++i) {

    EvalDerivedDeformations(DerivedDeformations, NodalCoordSystems, EvalDerivedShapeFunc( IntPoints_[i] ), IntPoints_, i);

    EvalJacobi(Jacobi, EvalDerivedShapeFunc( IntPoints_[i] ), NodalCoordSystems, IntPoints_, i);

    stCoord2globalCoord( Jacobi , DerivedDeformations);

    SetBMatrix(DerivedDeformations, B, ShellIntegration);
    
    mult( B,u,Strains[i] ); 
  }

  return Strains;
}

Dense_Matrix BaseShell::EvalStrainsAtIntPoints( const Dense_Matrix IntPoints_ ){
  return BaseShell::evalStrains(IntPoints_);
}

Dense_Matrix BaseShell::EvalStressesAtIntPoints( const Dense_Matrix IntPoints_, Material* MatPtr_, const float_type angle_ ){
  // get the element coordinate system stresses at the integration points.
  Dense_Matrix EDStresses = BaseShell::evalElementDirectionStresses( IntPoints_, MatPtr_, angle_);
  
  // initialization of the result matrix
  Dense_Matrix IntPointStresses(EDStresses.nrows(),EDStresses.ncols());
  mtl::set_value(IntPointStresses,0.);
 
  // iterate over all int points and determine the global stresses
  for (unsigned i = 0; i < IntPoints_.nrows(); ++i){
    // get the local coordinate system at the respective int point
    Dense_Matrix CSi(3,3);
    GetMaterialCoordSys( IntPoints_, i, CSi ) ; 
 
    // the material transformation matrix is corrected to make it usable for 
    // stress transformation. One could also implement this matrix
    Dense_Matrix M(6,6); //, T(6,6), Tinv(6,6), Tmp(6,6);
    M = GetStressTransformationMatrix( CSi );  
    //mult(Rinv,M,Tmp); 
    //mult(Tmp,R,Tinv);
    //lu_inversion(Tinv,T);
  
    // the global strains in a "vector"
    Dense_Vector elvec(6,0.);
    for ( unsigned j = 0; j < elvec.size(); ++j) elvec[j] = EDStresses(i,j);
    
    // the resulting local results in main material direction
    Dense_Vector globres(6,0.);
    mult(M,elvec,globres);
    
    // write the results into the mainmaterialdirectionstrains matrix
    for ( unsigned j = 0; j < globres.size(); ++j) IntPointStresses(i,j) = globres[j];  
  }
  return IntPointStresses;
}

Dense_Matrix BaseShell::Extrapolate2Corners( const Dense_Matrix IntPtValues_ ){
    
  Dense_Matrix CornerVals(IntPtValues_.nrows(),IntPtValues_.ncols()); 
  
  // extrapolation to the nodes of the element. 
  unsigned gnc = GetNodeCount();
  unsigned dim = GetDofSet().ElementDimension();
  Dense_Matrix IntPoints = GetIntPoints(); 
  unsigned nip = IntPoints.nrows();
  Dense_Vector p(nip);
  Dense_Matrix dummyM1(nip,nip), inv_dummyM1(nip,nip), dummyM2(gnc,nip);
  
  if (nip==8) {
    for (unsigned i = 0 ; i < nip; ++i) {
      dummyM1(i,0) = 1;
      dummyM1(i,1) = IntPoints(i,0);
      dummyM1(i,2) = IntPoints(i,1);
      dummyM1(i,3) = IntPoints(i,2);
      dummyM1(i,4) = IntPoints(i,0) * IntPoints(i,1);
      dummyM1(i,5) = IntPoints(i,1) * IntPoints(i,2);
      dummyM1(i,6) = IntPoints(i,0) * IntPoints(i,2);
      dummyM1(i,7) = IntPoints(i,0) * IntPoints(i,1) * IntPoints(i,2);
    }
  }
  else if (nip==4) {
    for (unsigned i = 0 ; i < nip; ++i) {
      dummyM1(i,0) = 1;
      dummyM1(i,1) = IntPoints(i,0);
      dummyM1(i,2) = IntPoints(i,1);
      dummyM1(i,3) = IntPoints(i,0) * IntPoints(i,1);
    }
  }
  lu_inversion(dummyM1, inv_dummyM1);

  if (nip==8) {
    for (unsigned i = 0 ; i < gnc; ++i) {
      dummyM2(i,0) = 1;
      dummyM2(i,1) = GetlNodeCoords()(i,0);
      dummyM2(i,2) = GetlNodeCoords()(i,1);
      dummyM2(i,3) = GetlNodeCoords()(i,2);
      dummyM2(i,4) = GetlNodeCoords()(i,0) * GetlNodeCoords()(i,1);
      dummyM2(i,5) = GetlNodeCoords()(i,1) * GetlNodeCoords()(i,2);
      dummyM2(i,6) = GetlNodeCoords()(i,0) * GetlNodeCoords()(i,2);
      dummyM2(i,7) = GetlNodeCoords()(i,0) * GetlNodeCoords()(i,1) * GetlNodeCoords()(i,2);
    }
  }
  else if (nip==4) {
    for (unsigned i = 0 ; i < gnc; ++i) {
      dummyM2(i,0) = 1;
      dummyM2(i,1) = GetlNodeCoords()(i,0);
      dummyM2(i,2) = GetlNodeCoords()(i,1);
      dummyM2(i,3) = GetlNodeCoords()(i,0) * GetlNodeCoords()(i,1);
    }
  }

  for (unsigned i=0; i<(dim-1)*3; i++)
    {
      mtl::set(p,0.0);
      mult(inv_dummyM1, columns(IntPtValues_)[i], p);

      mult(dummyM2, p, columns(CornerVals)[i]);

    }  

  return CornerVals;
}


void BaseShell::EvalStresses(string StressType) {
  std::cerr << GetName() << "::EvalStresses() no more implemented for shell elements." << '\n' 
            << "Shell stresses cannot be evaluated for the nodes in a reasonable way" << '\n'
            << "so that a single stress can be evaluated for each node. There are"    << '\n'
            << "some specialized functions available for the evaluation of stresses in" << '\n'
            << "any coordinate system (global, element, material direction). Check BaseShell.h !" 
            << std::endl; 
  
/*
std::cout << "Strains at int points" << std::endl;
print_all_matrix(EvalStrainsAtIntPoints(GetIntPoints()));

std::cout << "Stresses at int points" << std::endl;
print_all_matrix(EvalStressesAtIntPoints(GetIntPoints(), GetMaterialPtr(), GetPropertiesPtr()->GetDouble("Theta")*PI/180.));

std::cout << "element strains" << std::endl;
print_all_matrix(evalElementDirectionStrains(GetIntPoints()));

std::cout << "element stresses" << std::endl;
print_all_matrix(evalElementDirectionStresses(GetIntPoints(), GetMaterialPtr(), GetPropertiesPtr()->GetDouble("Theta")*PI/180.));

std::cout << "material strains" << std::endl;
print_all_matrix(evalMaterialDirectionStrains(GetIntPoints(), GetPropertiesPtr()->GetDouble("Theta")*PI/180.));

std::cout << "material stresses" << std::endl;
print_all_matrix(evalMaterialDirectionStresses(GetIntPoints(), GetMaterialPtr(), GetPropertiesPtr()->GetDouble("Theta")*PI/180.));
*/

}

// function determinig the material direction strains, i.e. e1, e2, ...
Dense_Matrix BaseShell::evalElementDirectionStrains( const Dense_Matrix IntPoints_ ){
  // get the global strains at the integration points.
  Dense_Matrix GlobalStrains = BaseShell::evalStrains( IntPoints_ );
  
  // the matrix containing the transformed strains
  Dense_Matrix ElementDirectionStrains(GlobalStrains.nrows(),6);
  
  // iterate over all int points and determine the local strains
  for (unsigned i = 0; i < IntPoints_.nrows(); ++i){
    // get the local coordinate system at the respective int point
    Dense_Matrix CSi(3,3);
    GetMaterialCoordSys( IntPoints_, i, CSi ) ;
 
    // the material transformation matrix inlcuding the reuter correction as defined by cook
    Dense_Matrix mattr(6,6);
    mattr = GetMaterialTransformationMatrix( CSi ); 
  
    // the global strains in a "vector"
    Dense_Vector globvec(6,0.);
    for ( unsigned j = 0; j < globvec.size(); ++j) globvec[j] = GlobalStrains(i,j);
    
    // the resulting local results in main material direction
    Dense_Vector locres(6,0.);
    mult(mattr,globvec,locres);
    
    // write the results into the mainmaterialdirectionstrains matrix
    for ( unsigned j = 0; j < locres.size(); ++j) ElementDirectionStrains(i,j) = locres[j];  
  }
  
  return ElementDirectionStrains;
}

Dense_Matrix BaseShell::evalMaterialDirectionStrains( const Dense_Matrix IntPoints_ , const float_type angle_){
  
  // get the element coordinate strains at the integration points.
  Dense_Matrix ElementDirectionStrains = BaseShell::evalElementDirectionStrains( IntPoints_ );
  
  // the matrix containing the material direction strains
  Dense_Matrix MaterialDirectionStrains(ElementDirectionStrains.nrows(),6);
  
  // the element direction strains need to be transformed to the material direction strains
  // this requires a transformation about the local zeta direction.
  for (unsigned i = 0; i < IntPoints_.nrows(); ++i){
    Dense_Matrix rot(3,3);
    mtl::copy(CoordSys::Vec2Mat(0.0,0.0,angle_), rot);  

    // the material transformation matrix inlcuding the reuter correction as defined by cook
    Dense_Matrix mattr(6,6);
    mattr = GetMaterialTransformationMatrix( rot ); 
  
    // the global strains in a "vector"
    Dense_Vector globvec(6,0.);
    for ( unsigned j = 0; j < globvec.size(); ++j) globvec[j] = ElementDirectionStrains(i,j);
    
    // the resulting local results in main material direction
    Dense_Vector locres(6,0.);
    mult(mattr,globvec,locres);
    
    // write the results into the mainmaterialdirectionstrains matrix
    for ( unsigned j = 0; j < locres.size(); ++j) MaterialDirectionStrains(i,j) = locres[j];  
 }
  return MaterialDirectionStrains;
}

Dense_Matrix BaseShell::evalElementDirectionStresses( const Dense_Matrix IntPoints_, Material* MatPtr_, const float_type angle_ ){
  // get the element coordinate strains at the integration points.
  Dense_Matrix MDStresses = BaseShell::evalMaterialDirectionStresses( IntPoints_, MatPtr_, angle_);

  // the matrix containing the material direction strains
  Dense_Matrix ElementDirectionStresses(MDStresses.nrows(),6);
  mtl::set_value(ElementDirectionStresses,0.);
  
  // get the transformation matrix
  Dense_Matrix transform(6,6);
  transform = Get3DMaterialTransformationMatrix(-angle_);
  
  // transpose the material direction stresses
  Dense_Matrix TransposedMDStresses(MDStresses.ncols(), MDStresses.nrows());
  transpose(MDStresses,TransposedMDStresses);
  
  Dense_Matrix TransposedEDStresses(ElementDirectionStresses.ncols(),ElementDirectionStresses.nrows());
  // multiply
  mult(transform,TransposedMDStresses,TransposedEDStresses);
  transpose(TransposedEDStresses,ElementDirectionStresses);
  
  return ElementDirectionStresses;
}

Dense_Matrix BaseShell::evalMaterialDirectionStresses( const Dense_Matrix IntPoints_, Material* MatPtr_ , const float_type angle_){
  // get the material direction strains for all int points
  Dense_Matrix MDStrains = evalMaterialDirectionStrains( IntPoints_, angle_ );

  /* this is a 6 x X matrix holding all the material direction strains for each int point
     int point 1: e1* e2* e3* g12* g23* g13*
     ...
     int point N: ...
  */

  // initialize result matrix
  Dense_Matrix MaterialDirectionStresses(MDStrains.nrows(),MDStrains.ncols());
  
  // the transposed matrix of MDStrains
  Dense_Matrix TransposedMDS(MDStrains.ncols(),MDStrains.nrows());
  transpose(MDStrains,TransposedMDS);
  
  Dense_Matrix D(6,6);
  if (MatPtr_->ClassName() == "TransverseIsotropicMaterial23" ){
    mtl::copy((MatPtr_->Get(Material::ThreeDShell)), D);
  }
  else if( MatPtr_->ClassName() == "IsotropicMaterial" ){
    throw std::logic_error("BaseShell::evalElementDirectionStresses: no implementation for isotropic materials yet, use transverse isotropic material instead!");
  }
  // if another material type is required.
  else {
    throw std::logic_error("BaseShell::evalElementDirectionStresses: not implemented for this material type");
  }

  // some corrections using a scaling factor
  Dense_Matrix NodalCoordSystems(3*IntPoints_.nrows(),3);
  GetNodalCoordSys( NodalCoordSystems );
  vector<Dense_Matrix> J_vec(IntPoints_.nrows());
  vector<double> detJ_vec(IntPoints_.nrows());

  double Volume = GetShellVolume( NodalCoordSystems, IntPoints_, J_vec, detJ_vec );
  double thickness = GetShellThickness();

  double f = 1.0 + 0.2*Volume/25.0/thickness/thickness/thickness;
  if ( f < 1.2 ) f = 1.2;

  D(4,4) /= f;
  D(5,5) /= f;

  // calculate the stresses by matrix multiplication. 
  Dense_Matrix TransposedMDStresses(TransposedMDS.nrows(), TransposedMDS.ncols());
  mult(D,TransposedMDS,TransposedMDStresses);
  transpose(TransposedMDStresses,MaterialDirectionStresses);
  
  return MaterialDirectionStresses;
}

Dense_Matrix BaseShell::evalMaxStressCriteria(const Dense_Matrix IntPoints_, const float_type angle_, Material* MatPtr_){
  
  // get the material direction stresses
  Dense_Matrix MDStresses = evalMaterialDirectionStresses(IntPoints_, MatPtr_, angle_);

  // initialize the resulting matrix
  Dense_Matrix Criteria(MDStresses.nrows(), 5);
  mtl::set_value(Criteria,0.0);
  
  // check whether the failure criteria values are properly defined
  if ( MatPtr_->Get("Xt") == 0. || MatPtr_->Get("Xc") == 0. ||
       MatPtr_->Get("Yt") == 0. || MatPtr_->Get("Yc") == 0. || 
       MatPtr_->Get("S") == 0.){
       cerr << "BaseShell::evalMaxStressCriteria: A failure criteria value" << endl; 
       cerr << "(ultimate tension or compression) is set to zero. Check the" << endl;
       cerr << "material definition! The current values are" << endl;
       cerr << "Xt " << MatPtr_->Get("Xt") << endl;
       cerr << "Xc " << MatPtr_->Get("Xc") << endl;
       cerr << "Yt " << MatPtr_->Get("Yt") << endl;
       cerr << "Yc " << MatPtr_->Get("Yc") << endl;
       cerr << "S " << MatPtr_->Get("S") << endl;
       exit(1);
       }
  
  // iterate over all integration points and evaluate maximum stress criteria
  for (unsigned i = 0; i < MDStresses.nrows(); ++i){
    if (MDStresses(i,0) >= 0.){ // in case of tension in 1 direction
      Criteria(i,0) = MDStresses(i,0)  / MatPtr_->Get("Xt");
    }
    else{ // in case of compression in 1 direction
      Criteria(i,1) = abs(MDStresses(i,0)) / MatPtr_->Get("Xc");
    }
    if (MDStresses(i,1) >= 0.){ // in case of tension in 2 direction
      Criteria(i,2) = MDStresses(i,1) / MatPtr_->Get("Yt");
    }
    else{ // in case of compression in 2 direction
      Criteria(i,3) = abs(MDStresses(i,1)) / MatPtr_->Get("Yc");
    }
    Criteria(i,4) = abs(MDStresses(i,3)) / MatPtr_->Get("S");  
  }
  
  return Criteria;
}

Dense_Matrix BaseShell::evalTsaiHillCriteria(const Dense_Matrix IntPoints_, const float_type angle_, Material* MatPtr_ ){
  // get the material direction stresses
  Dense_Matrix MDS = evalMaterialDirectionStresses(IntPoints_, MatPtr_, angle_);
  
  // initialize the resulting matrix
  Dense_Matrix Criteria(MDS.nrows(), 1);
  mtl::set_value(Criteria,0.0);
  
  // check whether the failure criteria values are properly defined
  if ( MatPtr_->Get("Xt") == 0. || MatPtr_->Get("Yt") == 0. || MatPtr_->Get("S") == 0.){
       cerr << "BaseShell::evalMaxStressCriteria: A failure criteria value" << endl; 
       cerr << "(ultimate tension or compression) is set to zero. Check the" << endl;
       cerr << "material definition! The current values are" << endl;
       cerr << "Xt " << MatPtr_->Get("Xt") << endl;
       cerr << "Xc " << MatPtr_->Get("Xc") << endl;
       cerr << "Yt " << MatPtr_->Get("Yt") << endl;
       cerr << "Yc " << MatPtr_->Get("Yc") << endl;
       cerr << "S " << MatPtr_->Get("S") << endl;
       exit(1);
       }
       
  // this criteria does not distinguish between tension and compression, thus
  // these values are set equal, i.e. Xt = Xc and Yt = Yc
  
  double X = MatPtr_->Get("Xt"); double Y = MatPtr_->Get("Yt"); double S = MatPtr_->Get("S");
  
  // iterate over all integration points and evaluate maximum stress criteria
  for (unsigned i = 0; i < MDS.nrows(); ++i){
    Criteria(i,0) = pow(MDS(i,0)/X,2) - MDS(i,0)*MDS(i,1)/pow(X,2) + pow(MDS(i,1)/Y,2) + pow(MDS(i,3)/S,2);
  }
  
  return Criteria;     
}

Dense_Matrix BaseShell::evalTsaiWuCriteria(const Dense_Matrix IntPoints_, const float_type angle_, Material* MatPtr_ ){
  // get the material direction stresses
  Dense_Matrix MDS = evalMaterialDirectionStresses(IntPoints_, MatPtr_, angle_);
  
  // initialize the resulting matrix
  Dense_Matrix Criteria(MDS.nrows(), 1);
  mtl::set_value(Criteria,0.0);
  
  // check whether the failure criteria values are properly defined
  if ( MatPtr_->Get("Xt") == 0. || MatPtr_->Get("Xc") == 0. ||
       MatPtr_->Get("Yt") == 0. || MatPtr_->Get("Yc") == 0. ||
       MatPtr_->Get("S") == 0.){
       cerr << "BaseShell::evalMaxStressCriteria: A failure criteria value" << endl; 
       cerr << "(ultimate tension or compression) is set to zero. Check the" << endl;
       cerr << "material definition! The current values are" << endl;
       cerr << "Xt " << MatPtr_->Get("Xt") << endl;
       cerr << "Xc " << MatPtr_->Get("Xc") << endl;
       cerr << "Yt " << MatPtr_->Get("Yt") << endl;
       cerr << "Yc " << MatPtr_->Get("Yc") << endl;
       cerr << "S " << MatPtr_->Get("S") << endl;
       exit(1);
       }
       
  double Xt = MatPtr_->Get("Xt"); double Xc = MatPtr_->Get("Xc");
  double Yt = MatPtr_->Get("Yt"); double Yc = MatPtr_->Get("Yc");
  double S = MatPtr_->Get("S");
  // evaluation of the constants
  double F1 = 1/Xt - 1/Xc; double F2 = 1/Yt - 1/Yc; double F11 = 1/(Xt*Xc);
  double F22 = 1/(Yt*Yc); double F66 = 1/pow(S,2);
  
  // The constant F12 is omitted, might lead to approx. 10% difference
  
  // iterate over all integration points and evaluate maximum stress criteria
  for (unsigned i = 0; i < MDS.nrows(); ++i){
    Criteria(i,0) = F1*MDS(i,0) + F2*MDS(i,1) + F11*pow(MDS(i,0),2) + F22*pow(MDS(i,1),2) + F66*pow(MDS(i,3),2);
  }
  
  return Criteria;     
}

Dense_Matrix BaseShell::evalHashinCriteria(const Dense_Matrix IntPoints_, const float_type angle_, Material* MatPtr_){
  // get the material direction stresses
  Dense_Matrix MDS = evalMaterialDirectionStresses(IntPoints_, MatPtr_, angle_);
  
  // initialize the resulting matrix
  Dense_Matrix Criteria(MDS.nrows(), 4);
  mtl::set_value(Criteria,0.0);
  
  // check whether the failure criteria values are properly defined
  if ( MatPtr_->Get("Xt") == 0. || MatPtr_->Get("Xc") == 0. ||
       MatPtr_->Get("Yt") == 0. || MatPtr_->Get("Yc") == 0. ||
       MatPtr_->Get("S") == 0.){
       cerr << "BaseShell::evalMaxStressCriteria: A failure criteria value" << endl; 
       cerr << "(ultimate tension or compression) is set to zero. Check the" << endl;
       cerr << "material definition! The current values are" << endl;
       cerr << "Xt " << MatPtr_->Get("Xt") << endl;
       cerr << "Xc " << MatPtr_->Get("Xc") << endl;
       cerr << "Yt " << MatPtr_->Get("Yt") << endl;
       cerr << "Yc " << MatPtr_->Get("Yc") << endl;
       cerr << "S " << MatPtr_->Get("S") << endl;
       exit(1);
       }
       
  double Xt = MatPtr_->Get("Xt"); double Xc = MatPtr_->Get("Xc");
  double Yt = MatPtr_->Get("Yt"); double Yc = MatPtr_->Get("Yc");
  // the original values Sq and Sl of the original Hashin criterion are assumed to be equal.
  double S = MatPtr_->Get("S");
    
  // iterate over all integration points and evaluate maximum stress criteria
  for (unsigned i = 0; i < MDS.nrows(); ++i){
    if (MDS(i,0) >= 0) Criteria(i,0) = pow(MDS(i,0)/Xt,2) + pow(MDS(i,3)/S,2);
    if (MDS(i,0) <  0) Criteria(i,1) = abs(MDS(i,0))/Xc;
    if (MDS(i,1) >= 0) Criteria(i,2) = pow(MDS(i,1)/Yt,2) + pow(MDS(i,3)/S,2);
    if (MDS(i,1) <  0) Criteria(i,3) = pow(MDS(i,1)/(2*S),2) + MDS(i,1)/Yc*(pow(Yc/(2*S),2)-1) - pow(MDS(i,3)/S,2);
  }
  
  return Criteria;     
}

/*
void BaseShell::SetMaterialTransformationVector()
{
        if ( GetIntPoints().nrows() == 8 ) //if it is GaussQuadratic3D
        {
                Dense_Matrix tmpCoordSys(3,3);
                MaterialTransformationVector.resize(4);
                MaterialTransformations2IntPts.resize(GetIntPoints().nrows());

                //Get the points in the plane of the nodes...
                Dense_Matrix tmpIntPoints(GetIntPoints().nrows(),GetIntPoints().ncols());
                mtl::copy(GetIntPoints(),tmpIntPoints);
                for (unsigned k = 0 ; k < tmpIntPoints.nrows() ; ++k )
                        tmpIntPoints(k,tmpIntPoints.ncols()-2) = 0.;

                for (unsigned i = 0 ; i < 4 ;  ++i)
                {
                        GetMaterialCoordSys( tmpIntPoints, i, tmpCoordSys );
                        MaterialTransformationVector[i] = mtl::Dense_Matrix(6,6);
                        mtl::copy(GetMaterialTransformationMatrix(tmpCoordSys),MaterialTransformationVector[i]);
                        MaterialTransformations2IntPts[i] = MaterialTransformations2IntPts[i+4] = i;
                }
        }
        else if ( GetIntPoints().nrows() == 12 ) //if it is Gauss2x2x3
        {
                Dense_Matrix tmpCoordSys(3,3);
                MaterialTransformationVector.resize(4);
                MaterialTransformations2IntPts.resize(GetIntPoints().nrows());

                for (unsigned i = 4 ; i < 8 ;  ++i)
                {
                        GetMaterialCoordSys( GetIntPoints(), i, tmpCoordSys );
                        MaterialTransformationVector[i-4] = mtl::Dense_Matrix(6,6);
                        mtl::copy(GetMaterialTransformationMatrix(tmpCoordSys),MaterialTransformationVector[i-4]);
                        MaterialTransformations2IntPts[i] =
                                MaterialTransformations2IntPts[i+4] =
                                MaterialTransformations2IntPts[i-4] = i-4;
                }
        }
        else if ( GetIntPoints().nrows() == 4 ) //if it is Gauss2x2x1
        {
                Dense_Matrix tmpCoordSys(3,3);
                MaterialTransformationVector.resize(4);
                MaterialTransformations2IntPts.resize(GetIntPoints().nrows());

                for (unsigned i = 0 ; i < 4 ;  ++i)
                {
                        GetMaterialCoordSys( GetIntPoints(), i, tmpCoordSys );
                        MaterialTransformationVector[i] = mtl::Dense_Matrix(6,6);
                        mtl::copy(GetMaterialTransformationMatrix(tmpCoordSys),MaterialTransformationVector[i]);
                        MaterialTransformations2IntPts[i] = i;
                }
        }
        else if ( GetIntPoints().nrows() == 6 ) //if it is Gauss3x2
        {
                Dense_Matrix tmpCoordSys(3,3);
                MaterialTransformationVector.resize(3);
                MaterialTransformations2IntPts.resize(GetIntPoints().nrows());

                //Get the points in the plane of the nodes...
                Dense_Matrix tmpIntPoints(GetIntPoints().nrows(),GetIntPoints().ncols());
                mtl::copy(GetIntPoints(),tmpIntPoints);
                for (unsigned k = 0 ; k < tmpIntPoints.nrows() ; ++k )
                        tmpIntPoints(k,tmpIntPoints.ncols()-2) = 0.;

                for (unsigned i = 0 ; i < 3 ;  ++i)
                {
                        GetMaterialCoordSys( tmpIntPoints, i, tmpCoordSys );
                        MaterialTransformationVector[i] = mtl::Dense_Matrix(6,6);
                        mtl::copy(GetMaterialTransformationMatrix(tmpCoordSys),MaterialTransformationVector[i]);
                        MaterialTransformations2IntPts[i] = MaterialTransformations2IntPts[i+3] = i;
                }
        }
        else if ( GetIntPoints().nrows() == 9 ) //if it is Gauss3x3
        {
                Dense_Matrix tmpCoordSys(3,3);
                MaterialTransformationVector.resize(3);
                MaterialTransformations2IntPts.resize(GetIntPoints().nrows());

                for (unsigned i = 3 ; i < 6 ;  ++i)
                {
                        GetMaterialCoordSys( GetIntPoints(), i, tmpCoordSys );
                        MaterialTransformationVector[i-3] = mtl::Dense_Matrix(6,6);
                        mtl::copy(GetMaterialTransformationMatrix(tmpCoordSys),MaterialTransformationVector[i-3]);
                        MaterialTransformations2IntPts[i] =
                                 MaterialTransformations2IntPts[i+3] =
                                 MaterialTransformations2IntPts[i-3] = i-3;
                }
        }
        else if ( GetIntPoints().nrows() == 3 ) //if it is AreaCoords3Points
        {
                Dense_Matrix tmpCoordSys(3,3);
                MaterialTransformationVector.resize(3);
                MaterialTransformations2IntPts.resize(GetIntPoints().nrows());

                for (unsigned i = 0 ; i < 3 ;  ++i)
                {
                        GetMaterialCoordSys( GetIntPoints(), i, tmpCoordSys );
                        MaterialTransformationVector[i] = mtl::Dense_Matrix(6,6);
                        mtl::copy(GetMaterialTransformationMatrix(tmpCoordSys),MaterialTransformationVector[i]);
                        MaterialTransformations2IntPts[i] = i;
                }
        }
        else
        {
                std::cout << "ERROR: In Shell93::SetMaterialTransformationVector() " << endl
                                << "ERROR: The requested Integration scheme is not implemented." << endl;
                exit(1);
        }
}
*/

//Functions not belonging to class BaseShell
//------------------------------------------

//Align x and y direction based on an element coordsys equal to ANSYS
void fe_base::AlignXY( Dense_Matrix CS_, const Dense_Vector normal_, const double toler_ ){

  // Check for dims of input args
  MTL_ASSERT( CS_.nrows() == 3 && CS_.ncols() == 3 && normal_.size() ==3, "felyx::AlignXY()");

  // Define local data containers
  Dense_Vector localX(3,0.0), localY(3,0.0);

  // Build the local y-direction
  mtl::cross_prod_3d(normal_, CS_[0], localY );
  if (mtl::two_norm( localY ) < toler_)
    mtl::cross_prod_3d(normal_, CS_[1], localY );

  double norm2 = 1.0 / two_norm( localY ) ;
  mtl::scale(localY, norm2  );

  // orient localX perpendicular to the other two
  mtl::cross_prod_3d( localY , normal_, localX );

  mtl::copy(localX , CS_[0]);
  mtl::copy(localY , CS_[1]);
  mtl::copy(normal_, CS_[2]);

}

Dense_Matrix fe_base::GetMaterialTransformationMatrix( Dense_Matrix CoordSys_ )
{
  Dense_Matrix MatTransform(6,6);
  mtl::set_value(MatTransform,0.0);
  unsigned n,p,i,k;

  //Check Cook for more details...
  for ( i = 0 ; i < 3 ; ++i ){
    for ( k = 0 ; k < 3 ; ++k ){
      p = k+1; n = i+1;
      if ( p == 3 ) p = 0; //trick to enable permutations
      if ( n == 3 ) n = 0; //trick to enable permutations

      MatTransform(i,k)     = CoordSys_(i,k)*CoordSys_(i,k); //Set the submatrix T11
      MatTransform(i,k+3)   = CoordSys_(i,k)*CoordSys_(i,p); //Set the submatrix T12
      MatTransform(i+3,k)   = 2*CoordSys_(i,k)*CoordSys_(n,k); //Set the submatrix T21
      MatTransform(i+3,k+3) = CoordSys_(i,k)*CoordSys_(n,p) +  CoordSys_(n,k)*CoordSys_(i,p); //Set T22
    }
  }

  return MatTransform;
}
                      
Dense_Matrix fe_base::GetStressTransformationMatrix( Dense_Matrix CoordSys_ )
{
  Dense_Matrix MatTransform(6,6);
  mtl::set_value(MatTransform,0.0);
  
  MatTransform = GetMaterialTransformationMatrix( CoordSys_ );
  
    // reuter matrix and its invers
  Dense_Matrix R(6,6), Rinv(6,6); 
  mtl::set_value(R,0.);
  R(0,0) = R(1,1) = R(2,2) = 1.;
  R(3,3) = R(4,4) = R(5,5) = 2.;
  lu_inversion(R,Rinv);
 
  // the material transformation matrix is corrected to make it usable for 
  // stress transformation.
  Dense_Matrix M(6,6), T(6,6), Tinv(6,6), Tmp(6,6);
  mtl::set_value(T,0.); mtl::set_value(Tinv,0.); mtl::set_value(Tmp,0.);
  M = GetMaterialTransformationMatrix( CoordSys_ );  
  mult(Rinv,MatTransform,Tmp); 
  mult(Tmp,R,Tinv);
  lu_inversion(Tinv,T);
  
  return T;
  
}
