//-----------------------------------------------------------------------------
// Darcy3D4.cc
//
// begin     : Sept 2 2003
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
#include "Darcy3D4.h"

using namespace fe_base;

////
//// Initialize static data members of class "Darcy3D4"
////
const int    Darcy3D4::NodeCount     	= 4;
const int    Darcy3D4::Id	       	= 470;
const string Darcy3D4::Name	       	= "LCM 3D Tetrahedron";
LcmDofSet Darcy3D4::ElementDofSet 	= LcmDofSet("100000");


// EvalDerivedShapeFunc is a Matrix of the values of the derived shape functions with respect to the standard
// coords L1, L2, L3 and L4 in a specific point.
// Its rows correspond to a specific shape function. The cols to the direction of derivation 
// EvalPoints is a Matrix of points, where the shape functions are to be evaluated. IP determines the specific point,
// means row of Evalpoints. The EvalPoints are given in volumetric coords
//----------------------------------------------------------------------------------------------------------------------

void Darcy3D4::EvalDerivedShapeFunc(Dense_Matrix shapefunc)
{
      shapefunc(0,0) = -1.0;
      shapefunc(1,0) = 1.0;
      shapefunc(2,0) = 0.0;
      shapefunc(3,0) = 0.0;

      shapefunc(0,1) = -1.0;
      shapefunc(1,1) = 0.0;
      shapefunc(2,1) = 1.0;
      shapefunc(3,1) = 0.0;

      shapefunc(0,2) = -1.0;
      shapefunc(1,2) = 0.0;
      shapefunc(2,2) = 0.0;
      shapefunc(3,2) = 1.0;
}

Dense_Matrix Darcy3D4::CalcEM() {

  double detJ;
  unsigned nSize = GetEMSize();
  unsigned dim = 3;

  //Initialization of some Matrices 
  Dense_Matrix 
    J(dim,dim),                     /* Jacobimatrix */
    K(nSize,nSize),                 /* "Stiffness" Matrix */
    shapefunc(GetNodeCount(), 3),   /* derived - shapefuncs */
    D,                              /* Permeabilty Matrix */
    dummy(4,3);

  //Get the material data
  D.resize(3,3);
  mtl::set_value(D,0.0);
  
  if ( MaterialPtr -> ClassName() == "TransverseIsotropicMaterial23" ) {
    mtl::copy(K_Mat,D);
  }
  else {
    D(0,0) = MaterialPtr->Get("K");
    D(1,1) = MaterialPtr->Get("K");
    D(2,2) = MaterialPtr->Get("K");
  }
  
  EvalDerivedShapeFunc(shapefunc);
  
  //Compute the Jacobian matrix
  EvalJacobi(J, shapefunc);
  
  //Transforme the shapefunctions to global x, y, z coordinates
  stCoord2globalCoord(J, shapefunc);
  
  detJ = determinant(J)/6.0;
  
  scale(D,detJ);
  
  mult(shapefunc,D,dummy);
  mult(dummy,trans(shapefunc),K);
  
  return K;
}


void Darcy3D4::EvalJacobi(Dense_Matrix J_, Dense_Matrix shapefunc_){

  //initialize the matrix to zero   
  mtl::set_value(J_,0.0);

    for (unsigned k = 0 ; k < GetNodeCount() ; ++k){
      J_(0,0) += shapefunc_(k,0)*NodeVec[k]->Cx;
      J_(0,1) += shapefunc_(k,0)*NodeVec[k]->Cy;
      J_(0,2) += shapefunc_(k,0)*NodeVec[k]->Cz;
      J_(1,0) += shapefunc_(k,1)*NodeVec[k]->Cx;

      J_(1,1) += shapefunc_(k,1)*NodeVec[k]->Cy;
      J_(1,2) += shapefunc_(k,1)*NodeVec[k]->Cz;
      J_(2,0) += shapefunc_(k,2)*NodeVec[k]->Cx;

      J_(2,1) += shapefunc_(k,2)*NodeVec[k]->Cy;
      J_(2,2) += shapefunc_(k,2)*NodeVec[k]->Cz;
    }
}



void Darcy3D4::stCoord2globalCoord(const Dense_Matrix Jacobi_,  Dense_Matrix shapefunc_){

  //inversion of the Jacobian
  unsigned size = Jacobi_.nrows();
  Dense_Matrix invJacobi(size, size);

  inversion(Jacobi_, invJacobi);

    Dense_Vector dummyvec1(size), dummyvec2(size);
    for (unsigned k = 0 ; k < GetNodeCount() ; k++)
      {
	dummyvec1[0] = shapefunc_(k,0);
	dummyvec1[1] = shapefunc_(k,1);
	dummyvec1[2] = shapefunc_(k,2);
	
	mtl::set_value(dummyvec2,0.0);
	mult(invJacobi, dummyvec1, dummyvec2);
	
	shapefunc_(k,0) = dummyvec2[0];
	shapefunc_(k,1) = dummyvec2[1];
	shapefunc_(k,2) = dummyvec2[2];
      }

}


void Darcy3D4::CompEleVelocity(){
  
  vx=0; vy=0; vz=0;

  Dense_Vector ab(3,0.0), ac(3,0.0), ad(3,0.0);
  ab[0] = GetNodeIter(1)->Cx - GetNodeIter(0)->Cx;
  ab[1] = GetNodeIter(1)->Cy - GetNodeIter(0)->Cy;
  ab[2] = GetNodeIter(1)->Cz - GetNodeIter(0)->Cz;
  ac[0] = GetNodeIter(2)->Cx - GetNodeIter(0)->Cx;
  ac[1] = GetNodeIter(2)->Cy - GetNodeIter(0)->Cy;
  ac[2] = GetNodeIter(2)->Cz - GetNodeIter(0)->Cz;
  ad[0] = GetNodeIter(3)->Cx - GetNodeIter(0)->Cx;
  ad[1] = GetNodeIter(3)->Cy - GetNodeIter(0)->Cy;
  ad[2] = GetNodeIter(3)->Cz - GetNodeIter(0)->Cz;
  
  double B, C, D, E;
  E = (ab[0] * ac[1] * ad[2] + ab[1] * ac[2] * ad[0] + ab[2] * ac[0] * ad[1]
       - ab[2] * ac[1] * ad[0] - ab[0] * ac[2] * ad[1] - ab[1] * ac[0] * ad[2]);

  Dense_Matrix temp(3,3), temp2(3,3) ;

  for (unsigned i=0; i<GetNodeCount(); ++i) {

    temp(0,0) = GetNodeIter((i+1)%4)->Cx;
    temp(1,0) = GetNodeIter((i+2)%4)->Cx;
    temp(2,0) = GetNodeIter((i+3)%4)->Cx;
    temp(0,1) = GetNodeIter((i+1)%4)->Cy;
    temp(1,1) = GetNodeIter((i+2)%4)->Cy;
    temp(2,1) = GetNodeIter((i+3)%4)->Cy;
    temp(0,2) = GetNodeIter((i+1)%4)->Cz;
    temp(1,2) = GetNodeIter((i+2)%4)->Cz;
    temp(2,2) = GetNodeIter((i+3)%4)->Cz;
    
    mtl::copy(temp,temp2);
    temp2(0,0)=1.0;
    temp2(1,0)=1.0;
    temp2(2,0)=1.0;
    
    B = determinant(temp2); 

    mtl::copy(temp,temp2);
    temp2(0,1)=1.0;
    temp2(1,1)=1.0;
    temp2(2,1)=1.0;
    
    C = determinant(temp2);

    mtl::copy(temp,temp2);
    temp2(0,2)=1.0;
    temp2(1,2)=1.0;
    temp2(2,2)=1.0;
    
    D = determinant(temp2);

    if (i%2==1) {
      B=-B; C=-C; D=-D;
    }

    vx += GetNodeIter(i)->Pressure / E * B;
    vy += GetNodeIter(i)->Pressure / E * C;
    vz += GetNodeIter(i)->Pressure / E * D;
  }
  if ( MaterialPtr -> ClassName() == "TransverseIsotropicMaterial23" ){
    if ( EleCoordSysPtr == NULL){
      vx *= MaterialPtr->Get("K1")/MaterialPtr->Get("visc");
      vy *= MaterialPtr->Get("K2")/MaterialPtr->Get("visc");
      vz *= MaterialPtr->Get("K3")/MaterialPtr->Get("visc");
    }
    else {
      double visc = MaterialPtr->Get("visc");
      double A = vx;
      double B = vy;
      double C = vz;
      vx = (K_Mat(0,0)*A + K_Mat(0,1)*B + K_Mat(0,2)*C)/visc;
      vy = (K_Mat(1,0)*A + K_Mat(1,1)*B + K_Mat(1,2)*C)/visc;
      vz = (K_Mat(2,0)*A + K_Mat(2,1)*B + K_Mat(2,2)*C)/visc;
    }
  }
  else {
    vx *= MaterialPtr->Get("K")/MaterialPtr->Get("visc");
    vy *= MaterialPtr->Get("K")/MaterialPtr->Get("visc");
    vz *= MaterialPtr->Get("K")/MaterialPtr->Get("visc");
  }

  // cout << "vx = " << vx << ", vy = " << vy << ", vz = " << vz << endl;

}


void Darcy3D4::CompVolFlux() {
  Dense_Vector norm(3,0.0), v(3,0.0);
  v[0]=vx; v[1]=vy; v[2]=vz;

  Flux(0,1,v,ssab,ssc);
  Flux(0,2,v,ssac,ssd);
  Flux(0,3,v,ssad,ssb);
  Flux(1,2,v,ssbc,ssa);
  Flux(1,3,v,ssbd,ssc);
  Flux(2,3,v,sscd,ssa);
}

void Darcy3D4::Flux(unsigned n1, unsigned n2, Dense_Vector v, Dense_Vector s1, Dense_Vector s2){

  Dense_Vector norm(3,0.0);
  cross_prod_3d(s1,s2,norm);
  double q_ = scalar_prod_3d(norm,v);
  if (q_>0) q_ *= GetNodeIter(n2)->Fill;
  else      q_ *= GetNodeIter(n1)->Fill;
  GetNodeIter(n1)->q += q_;
  GetNodeIter(n2)->q -= q_;
}

void Darcy3D4::EvalCentroid(){
  unsigned n = GetNodeCount();
  for (unsigned i=0; i<n; ++i) {
    sx += GetNodeIter(i)->Cx;
    sy += GetNodeIter(i)->Cy;
    sz += GetNodeIter(i)->Cz;
  }
  sx /= n; sy /= n; sz /= n;

  Dense_Vector s(3,0.0), sa(3,0.0), sb(3,0.0), sc(3,0.0), sd(3,0.0),
    sab(3,0.0), sac(3,0.0), sad(3,0.0), sbc(3,0.0), sbd(3,0.0), scd(3,0.0);

  s[0]=sx; s[1]=sy; s[2]=sz; 

  sab[0]= (GetNodeIter(0)->Cx + GetNodeIter(1)->Cx)/2.0;
  sab[1]= (GetNodeIter(0)->Cy + GetNodeIter(1)->Cy)/2.0;
  sab[2]= (GetNodeIter(0)->Cz + GetNodeIter(1)->Cz)/2.0;

  sac[0]= (GetNodeIter(0)->Cx + GetNodeIter(2)->Cx)/2.0;
  sac[1]= (GetNodeIter(0)->Cy + GetNodeIter(2)->Cy)/2.0;
  sac[2]= (GetNodeIter(0)->Cz + GetNodeIter(2)->Cz)/2.0;

  sad[0]= (GetNodeIter(0)->Cx + GetNodeIter(3)->Cx)/2.0;
  sad[1]= (GetNodeIter(0)->Cy + GetNodeIter(3)->Cy)/2.0;
  sad[2]= (GetNodeIter(0)->Cz + GetNodeIter(3)->Cz)/2.0;

  sbc[0]= (GetNodeIter(1)->Cx + GetNodeIter(2)->Cx)/2.0;
  sbc[1]= (GetNodeIter(1)->Cy + GetNodeIter(2)->Cy)/2.0;
  sbc[2]= (GetNodeIter(1)->Cz + GetNodeIter(2)->Cz)/2.0;

  sbd[0]= (GetNodeIter(1)->Cx + GetNodeIter(3)->Cx)/2.0;
  sbd[1]= (GetNodeIter(1)->Cy + GetNodeIter(3)->Cy)/2.0;
  sbd[2]= (GetNodeIter(1)->Cz + GetNodeIter(3)->Cz)/2.0;

  scd[0]= (GetNodeIter(2)->Cx + GetNodeIter(3)->Cx)/2.0;
  scd[1]= (GetNodeIter(2)->Cy + GetNodeIter(3)->Cy)/2.0;
  scd[2]= (GetNodeIter(2)->Cz + GetNodeIter(3)->Cz)/2.0;

  sa[0]= (GetNodeIter(1)->Cx + GetNodeIter(2)->Cx + GetNodeIter(3)->Cx)/3.0;
  sa[1]= (GetNodeIter(1)->Cy + GetNodeIter(2)->Cy + GetNodeIter(3)->Cy)/3.0;
  sa[2]= (GetNodeIter(1)->Cz + GetNodeIter(2)->Cz + GetNodeIter(3)->Cz)/3.0;

  sb[0]= (GetNodeIter(0)->Cx + GetNodeIter(2)->Cx + GetNodeIter(3)->Cx)/3.0;
  sb[1]= (GetNodeIter(0)->Cy + GetNodeIter(2)->Cy + GetNodeIter(3)->Cy)/3.0;
  sb[2]= (GetNodeIter(0)->Cz + GetNodeIter(2)->Cz + GetNodeIter(3)->Cz)/3.0;

  sc[0]= (GetNodeIter(0)->Cx + GetNodeIter(1)->Cx + GetNodeIter(3)->Cx)/3.0;
  sc[1]= (GetNodeIter(0)->Cy + GetNodeIter(1)->Cy + GetNodeIter(3)->Cy)/3.0;
  sc[2]= (GetNodeIter(0)->Cz + GetNodeIter(1)->Cz + GetNodeIter(3)->Cz)/3.0;

  sd[0]= (GetNodeIter(0)->Cx + GetNodeIter(1)->Cx + GetNodeIter(2)->Cx)/3.0;
  sd[1]= (GetNodeIter(0)->Cy + GetNodeIter(1)->Cy + GetNodeIter(2)->Cy)/3.0;
  sd[2]= (GetNodeIter(0)->Cz + GetNodeIter(1)->Cz + GetNodeIter(2)->Cz)/3.0;

  ssa.resize(3); mtl::copy(sa,ssa); add(scaled(s,-1.0),ssa);
  ssb.resize(3); mtl::copy(sb,ssb); add(scaled(s,-1.0),ssb);
  ssc.resize(3); mtl::copy(sc,ssc); add(scaled(s,-1.0),ssc);
  ssd.resize(3); mtl::copy(sd,ssd); add(scaled(s,-1.0),ssd);

  ssab.resize(3); mtl::copy(sab,ssab); add(scaled(s,-1.0),ssab);
  ssac.resize(3); mtl::copy(sac,ssac); add(scaled(s,-1.0),ssac);
  ssad.resize(3); mtl::copy(sad,ssad); add(scaled(s,-1.0),ssad);
  ssbc.resize(3); mtl::copy(sbc,ssbc); add(scaled(s,-1.0),ssbc);
  ssbd.resize(3); mtl::copy(sbd,ssbd); add(scaled(s,-1.0),ssbd);
  sscd.resize(3); mtl::copy(scd,sscd); add(scaled(s,-1.0),sscd);

   if ( MaterialPtr -> ClassName() == "TransverseIsotropicMaterial23" ) {
     K_Mat.resize(3,3); 
     Dense_Matrix D(3,3);
 
    D(0,0) = MaterialPtr->Get("K1");
    D(1,1) = MaterialPtr->Get("K2");
    D(2,2) = MaterialPtr->Get("K3");
    if (EleCoordSysPtr==NULL){
       mtl::copy(D,K_Mat);}
    else{
       mult_AT_B_A_add(EleCoordSysPtr->GetRotMat(),D,K_Mat);}
  }
 }


void Darcy3D4::EvalCV() {
  Dense_Vector ab(3,0.0), ac(3,0.0), ad(3,0.0), temp(3,0.0);
  
  ab[0] = GetNodeIter(1)->Cx - GetNodeIter(0)->Cx;
  ab[1] = GetNodeIter(1)->Cy - GetNodeIter(0)->Cy;
  ab[2] = GetNodeIter(1)->Cz - GetNodeIter(0)->Cz;
  ac[0] = GetNodeIter(2)->Cx - GetNodeIter(0)->Cx;
  ac[1] = GetNodeIter(2)->Cy - GetNodeIter(0)->Cy;
  ac[2] = GetNodeIter(2)->Cz - GetNodeIter(0)->Cz;
  ad[0] = GetNodeIter(3)->Cx - GetNodeIter(0)->Cx;
  ad[1] = GetNodeIter(3)->Cy - GetNodeIter(0)->Cy;
  ad[2] = GetNodeIter(3)->Cz - GetNodeIter(0)->Cz;
  
  double d = (ab[0] * ac[1] * ad[2] + ab[1] * ac[2] * ad[0] + ab[2] * ac[0] * ad[1]
	      - ab[2] * ac[1] * ad[0] - ab[0] * ac[2] * ad[1] - ab[1] * ac[0] * ad[2]) / 24.0;
  
  for (unsigned i=0; i<GetNodeCount(); ++i) {
    GetNodeIter(i)->Volume += d;
  }
  
}
