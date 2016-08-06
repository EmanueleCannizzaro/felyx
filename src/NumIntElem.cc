//-----------------------------------------------------------------------------
// NumIntElem.cc
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
   

#include "NumIntElem.h"


using namespace fe_base;

void NumIntElem::EvalShapeFunc(Dense_Matrix A, const Dense_Matrix B, int n){
  std::string error = GetName() +"::EvalShapeFunc(...) is not implemented, aborting!";
  FELYX_RUNTIME_THROW( error.c_str() );
}

Dense_Matrix NumIntElem::GetIntegrationPoints(IntegrationScheme type)
{
  switch (type)
    {
    case AreaCoords3Points :
      {
	//---------------------------------------------------------
	//Triangular 2D integration according to Ansys.
	//The cols are the areacoords L1, L2, L3 of the specific
	//integration points. The areal weights are stored,
	//in the fourth  col. The cordinates given, are those used
	//by ansys for triangular plane Elements
	//---------------------------------------------------------
	Dense_Matrix IP(3,4);
	
	double coord1 = 2.0/3.0;
	double coord2 = 0.16666666666666666666666666;
	//the weight is divided by 6 because the determinant gives the area of a
	//rectangle and a triangle  has a 2 times smaller volume (1.0/3.0/2.0 = 0.166666)
	double weight = 0.1666666666666666666666666;
	
	IP(0,0) = coord1; IP(0,1) = coord2; IP(0,2) = coord2; IP(0,3) = weight;
	IP(1,0) = coord2; IP(1,1) = coord1; IP(1,2) = coord2; IP(1,3) = weight;
	IP(2,0) = coord2; IP(2,1) = coord2; IP(2,2) = coord1; IP(2,3) = weight;
	
	return IP;
      	break;
      }
    case Gauss3x3 :
      {
	//---------------------------------------------------------
	//Triangular 2D integration according to Ansys.
	//The cols are the areacoords L1, L2, L3 of the specific
	//integration points. The areal weights are stored,
	//in the fourth  col. The cordinates given, are those used
	//by ansys for triangular plane Elements
	//---------------------------------------------------------
	Dense_Matrix IP(9,5);

	double thCoord = sqrt(0.6);
	double weight1 = 8.0/9.0;
	double weight2 = 5.0/9.0;
	

	double coord1 = 2.0/3.0;
	double coord2 = 0.16666666666666666666666666;
	//the weight is divided by 6 because the determinant gives the area of a
	//rectangle and a triangle  has a 2 times smaller area (1.0/3.0/2.0 = 0.166666)
	double weight = 0.1666666666666666666666666;
	
	IP(0,0) = coord1; IP(0,1) = coord2; IP(0,2) = coord2; IP(0,3) = -thCoord; IP(0,4) = weight*weight2;
	IP(1,0) = coord2; IP(1,1) = coord1; IP(1,2) = coord2; IP(1,3) = -thCoord; IP(1,4) = weight*weight2;
	IP(2,0) = coord2; IP(2,1) = coord2; IP(2,2) = coord1; IP(2,3) = -thCoord; IP(2,4) = weight*weight2;

	IP(3,0) = coord1; IP(3,1) = coord2; IP(3,2) = coord2; IP(3,3) = 0       ; IP(3,4) = weight*weight1;
	IP(4,0) = coord2; IP(4,1) = coord1; IP(4,2) = coord2; IP(4,3) = 0       ; IP(4,4) = weight*weight1;
	IP(5,0) = coord2; IP(5,1) = coord2; IP(5,2) = coord1; IP(5,3) = 0       ; IP(5,4) = weight*weight1;

	IP(6,0) = coord1; IP(6,1) = coord2; IP(6,2) = coord2; IP(6,3) =  thCoord; IP(6,4) = weight*weight2;
	IP(7,0) = coord2; IP(7,1) = coord1; IP(7,2) = coord2; IP(7,3) =  thCoord; IP(7,4) = weight*weight2;
	IP(8,0) = coord2; IP(8,1) = coord2; IP(8,2) = coord1; IP(8,3) =  thCoord; IP(8,4) = weight*weight2;

	
	return IP;
      	break;
      }    case Gauss3x2 :
      {
	//---------------------------------------------------------
	//Triangular 2D integration according to Ansys.
	//The cols are the areacoords L1, L2, L3 of the specific
	//integration points. The areal weights are stored,
	//in the fourth  col. The cordinates given, are those used
	//by ansys for triangular plane Elements
	//---------------------------------------------------------
	Dense_Matrix IP(6,5);
	//thickness weight is 1.0
	double thCoord = 1.0 / sqrt(3.0);

	double coord1 = 2.0/3.0;
	double coord2 = 0.16666666666666666666666666;
	//the weight is divided by 6 because the determinant gives the area of a
	//rectangle and a triangle  has a 2 times smaller area (1.0/3.0/2.0 = 0.166666)
	double weight = 0.1666666666666666666666666;
	
	IP(0,0) = coord1; IP(0,1) = coord2; IP(0,2) = coord2; IP(0,3) = -thCoord; IP(0,4) = weight;
	IP(1,0) = coord2; IP(1,1) = coord1; IP(1,2) = coord2; IP(1,3) = -thCoord; IP(1,4) = weight;
	IP(2,0) = coord2; IP(2,1) = coord2; IP(2,2) = coord1; IP(2,3) = -thCoord; IP(2,4) = weight;

	IP(3,0) = coord1; IP(3,1) = coord2; IP(3,2) = coord2; IP(3,3) =  thCoord; IP(3,4) = weight;
	IP(4,0) = coord2; IP(4,1) = coord1; IP(4,2) = coord2; IP(4,3) =  thCoord; IP(4,4) = weight;
	IP(5,0) = coord2; IP(5,1) = coord2; IP(5,2) = coord1; IP(5,3) =  thCoord; IP(5,4) = weight;

	
	return IP;
      	break;
      }

    case VolumeCoords4Points :
      {
	//----------------------------------------------------------
	//Tetraedral 3D integration.
	//The cols are the volumecoords L1, L2, L3, L4 of the
	//integration points. The volumial weights are stored in the
	//fifth col.
	//-----------------------------------------------------------
	Dense_Matrix IP(4,5);
	
	double coord1 = 0.585410196624968;
	double coord2 = 0.138196601125010;
	//the weight is divided by 6 because the determinant gives the volume of a
	//qube and a tetraedron has a 6 times smaller volume (0.25/6.0 = 0.04166666) 
	double weight = 0.0416666666666666;
	
	IP(0,0) = coord1; IP(0,1) = coord2; IP(0,2) = coord2; IP(0,3) = coord2; IP(0,4) = weight;
	IP(1,0) = coord2; IP(1,1) = coord1; IP(1,2) = coord2; IP(1,3) = coord2; IP(1,4) = weight;
	IP(2,0) = coord2; IP(2,1) = coord2; IP(2,2) = coord1; IP(2,3) = coord2; IP(2,4) = weight;
	IP(3,0) = coord2; IP(3,1) = coord2; IP(3,2) = coord2; IP(3,3) = coord1; IP(3,4) = weight;
	
	return IP;  
	break;
      }

    case GaussQuadratic3D :
      {
	double coord = 1.0 / sqrt(3.0);
	double weight = 1.0;

	//----------------------------------------------------------
	//3D - Gauss-Integration with 2 points per dimension (quadratic)
	//The cols are the rectangular xi, eta, zeta coords of the points
	//The fourth col stores the three dimensional weighting of the points
	//----------------------------------------------------------
	Dense_Matrix IP(8,4); 
	
	IP(0,0) = -coord; IP(0,1) = -coord; IP(0,2) = -coord; IP(0,3) = weight;
	IP(1,0) =  coord; IP(1,1) = -coord; IP(1,2) = -coord; IP(1,3) = weight;
	IP(2,0) =  coord; IP(2,1) =  coord; IP(2,2) = -coord; IP(2,3) = weight;
	IP(3,0) = -coord; IP(3,1) =  coord; IP(3,2) = -coord; IP(3,3) = weight;
	IP(4,0) = -coord; IP(4,1) = -coord; IP(4,2) =  coord; IP(4,3) = weight;
	IP(5,0) =  coord; IP(5,1) = -coord; IP(5,2) =  coord; IP(5,3) = weight;
	IP(6,0) =  coord; IP(6,1) =  coord; IP(6,2) =  coord; IP(6,3) = weight;
	IP(7,0) = -coord; IP(7,1) =  coord; IP(7,2) =  coord; IP(7,3) = weight;
	
	return IP;
	break;
      }
    case GaussLinear3D :
      {
	Dense_Matrix IP(1,4);
	mtl::set_value(IP, 0.0);
	IP(0,3) = 8.0;
	return IP;

	break;
      }
    case GaussQuadratic2D :	
      {

	double coord = 0.577350269189626;
	double weight = 1.0;

	//----------------------------------------------------------
	//2D - Gauss-Integration with 2 points per dimension (quadratic)
	//The cols are the rectangular xi, eta coords of the points
	//The third col stores the two dimensional weighting of the points
	//----------------------------------------------------------
	Dense_Matrix IP(4,3);
	IP(0,0) = -coord; IP(0,1) = -coord; IP(0,2) = weight;
	IP(1,0) =  coord; IP(1,1) = -coord; IP(1,2) = weight;
	IP(2,0) =  coord; IP(2,1) =  coord; IP(2,2) = weight;
	IP(3,0) = -coord; IP(3,1) =  coord; IP(3,2) = weight;
	
	return IP;
	break;
      }

    case GaussQubic2D :
      {
	double coord = 0.774596669241483; 
	double weight1 = 0.888888888888889; 
	double weight2 = 0.555555555555556;
	
	//----------------------------------------------------------
	//2D - Gauss-Integration with 3 points per dimension (qubic)
	//The cols are the rectangular xi and  eta coords of the points
	//The fourth col stores the two  dimensional weighting of the points
	//----------------------------------------------------------

	Dense_Matrix IP(9,3);
	
	IP(0,0) = -coord;    IP(0,1) = -coord;   IP(0,2) = weight2*weight2;
	IP(1,0) = 0;         IP(1,1) = -coord;   IP(1,2) = weight1*weight2;
	IP(2,0) = coord;     IP(2,1) = -coord;   IP(2,2) = weight2*weight2;
	IP(3,0) = -coord;    IP(3,1) = 0;        IP(3,2) = weight1*weight2;
	IP(4,0) = 0;         IP(4,1) = 0;        IP(4,2) = weight1*weight1;
	IP(5,0) = coord;     IP(5,1) = 0;        IP(5,2) = weight1*weight2;
	IP(6,0) = -coord;    IP(6,1) = coord;    IP(6,2) = weight2*weight2;
	IP(7,0) = 0;         IP(7,1) = coord;    IP(7,2) = weight1*weight2;
	IP(8,0) = coord;     IP(8,1) = coord;    IP(8,2) = weight2*weight2;
	
	return IP;
	break;
      }
    case Gauss2x2x3 :
      {
	const double thCoord = sqrt(0.6);
	const double aCoord = 1.0/(sqrt(3.0));
	const double weight1 = 8.0/9.0;
	const double weight2 = 5.0/9.0;
	
	//----------------------------------------------------------
	//3D - Gauss-Integration with 2x2 points in plane and 3 in thickness direction
	//The cols are the rectangular xi and  eta coords of the points
	//The fourth col stores the two  dimensional weighting of the points
	//----------------------------------------------------------
	
	Dense_Matrix IP(12,4);
	
	IP(0,0)=-aCoord; IP(0,1)=-aCoord; IP(0,2)=-thCoord; IP(0,3)=weight2;
	IP(1,0)= aCoord; IP(1,1)=-aCoord; IP(1,2)=-thCoord; IP(1,3)=weight2;
	IP(2,0)= aCoord; IP(2,1)= aCoord; IP(2,2)=-thCoord; IP(2,3)=weight2;
	IP(3,0)=-aCoord; IP(3,1)= aCoord; IP(3,2)=-thCoord; IP(3,3)=weight2;

	IP(4,0)=-aCoord; IP(4,1)=-aCoord; IP(4,2)= 0.0; IP(4,3)=weight1;
	IP(5,0)= aCoord; IP(5,1)=-aCoord; IP(5,2)= 0.0; IP(5,3)=weight1;
	IP(6,0)= aCoord; IP(6,1)= aCoord; IP(6,2)= 0.0; IP(6,3)=weight1;
	IP(7,0)=-aCoord; IP(7,1)= aCoord; IP(7,2)= 0.0; IP(7,3)=weight1;

	IP( 8,0)=-aCoord; IP( 8,1)=-aCoord; IP( 8,2)= thCoord; IP( 8,3)=weight2;
	IP( 9,0)= aCoord; IP( 9,1)=-aCoord; IP( 9,2)= thCoord; IP( 9,3)=weight2;
	IP(10,0)= aCoord; IP(10,1)= aCoord; IP(10,2)= thCoord; IP(10,3)=weight2;
	IP(11,0)=-aCoord; IP(11,1)= aCoord; IP(11,2)= thCoord; IP(11,3)=weight2;



	return IP;
	break;
      }
    case Gauss2x2x1 :
      {
	const double aCoord = 1.0/(sqrt(3.0));
	const double weight = 2.0;
	
	//----------------------------------------------------------
	//3D - Gauss-Integration with 2x2 points in plane and 3 in thickness direction
	//The cols are the rectangular xi and  eta coords of the points
	//The fourth col stores the two  dimensional weighting of the points
	//----------------------------------------------------------
	
	Dense_Matrix IP(4,4);
	
	IP(0,0)=-aCoord; IP(0,1)=-aCoord; IP(0,2)= 0.0; IP(0,3)=weight;
	IP(1,0)= aCoord; IP(1,1)=-aCoord; IP(1,2)= 0.0; IP(1,3)=weight;
	IP(2,0)= aCoord; IP(2,1)= aCoord; IP(2,2)= 0.0; IP(2,3)=weight;
	IP(3,0)=-aCoord; IP(3,1)= aCoord; IP(3,2)= 0.0; IP(3,3)=weight;

	return IP;
	break;
      }
    default :
      {
	Dense_Matrix dummy;
	return dummy;
      }

    } // of switch....

}  //of GetIntegrationPoints....

