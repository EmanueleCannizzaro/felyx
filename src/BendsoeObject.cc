//-----------------------------------------------------------------------------
// StructObject.cc
//
// begin     : Jan 2 2003
// copyright : (c) 2001 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {okoenig, wintermantel, nzehnder}@imes.mavt.ethz.ch
// www       : www.imes.ethz.ch/st
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

#include "BendsoeObject.h"
#include "PreProcessing.h"

using namespace felyx;



////
// Print global status of actual FEM evaluation
////
void BendsoeObject::PrintGlobalStatus(){
  OutStream << "# Print Global Status: " << endl;
  OutStream << "\t - Model loaded from        : " << InterfacePtr->GetLoadPath() << endl; 
  OutStream << "\t - # Elements               : " << Elements.size() << endl;
  OutStream << "\t - # Nodes                  : " << Nodes.size() << endl;
  OutStream << "\t - # Materials              : " << Materials.size() << endl;
  OutStream << "\t - # Active DOF's           : " << DofCount << endl;
  OutStream << "\t - # Memory needs GSM [MB]  : " << (int)(Profile * sizeof(float_type) / (1024*1024) ) << endl;
#ifdef USE_FLOATS
  OutStream << "\t - Precision of val type    : float" << endl;
#else
  OutStream << "\t - Precision of val type    : double" << endl;
#endif
}


////
// Function which calculates the Bendsoe
////

  
int BendsoeObject::CalcBend(){

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Link node DOF's to GSM index" << endl; 
  DofCount = LinkNodes2Gsm(Nodes, Elements);


  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Eval envelope and profile of GSM : ";

  dense1D<unsigned> Envelope(DofCount);
  Profile = EvalEnvelope(Envelope, Elements);

  if (noise > 0)
    OutStream << Profile
	      << " -> Memory needs of GSM: "
	      << (int)( Profile * sizeof(float_type) / (1024*1024)) << " MB " << endl;

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Initialize GSM of size : " << DofCount << endl;
  // typedef mtl::matrix< float_type , symmetric<lower>, envelope<>, row_major >::type EnvelopeMatrix;
  EnvelopeMatrix GSM(Envelope, Profile );

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Evaluate ESM's and assemble them to GSM" << endl;
  AssembleGM(Elements, GSM);

  if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Compute Bendsoe" << endl;
  int solved = topopt(GSM,Envelope);

   if (noise > 0) OutStream << "--< " << my_timer.elapsed() << " >-- " << "Solution done - status : " << solved << endl;

  return solved;

}

//-----------------------Calculate Deformations u--------------------------------------------------------------------
int BendsoeObject::FE(Vector x, dense1D<unsigned> Envelope, EnvelopeMatrix &K, Vector DofSolution, Vector Eold, Vector nuold)
{
 Vector F(Envelope.size());

 if (loop==1)
   {
     for (unsigned i = 0; i < Esize; i++)
       {
	 IsotropicMaterial* mat;                     // mat is a pointer
	 mat = new IsotropicMaterial;                // create pointer
	 Elements[i]->SetMaterialPtr(mat);
	 Elements[i]->GetMaterialPtr()->Set("nu",nuold[i]);  // set nu	  
       }
   }

     for (unsigned i = 0; i < Esize; i++)
       {	 
	 double Enew=Eold[i]*pow(x[i],3);                    // calculate new E-modul
	 Elements[i]->GetMaterialPtr()->Set("E",Enew);       // set new E-Modul
       }
	
      mtl::copy(DofSolution,F);                         
      mtl::set(K,0.0);                                // set K (=GSM) to 0 
      AssembleGM(Elements,K);                         // assemble new GSM
      skyline_solve(K, F, Envelope);                  // solve system of equation Ku=F 

	  for (unsigned i = 0; i < Nodes.size() ; i++)  
	    {
	      Nodes[i].SetDeformations(F);	      // assign new calculated deformations of nodes 
	    }
     return 0;
 }
//---------------------Optimality Criteria Update---------------------------------------------------------------------
int  BendsoeObject::calc_xnew( Vector &x, Vector dc, unsigned num, vector<unsigned>  x_fix, double total_x)
{
     double l1=0, l2=100000;
     double lmid;
     float_type  min_x=0.01;              // min. value of density, not smaller, because system is getting singular
     float_type  move=0.2;                // move-limit for density update
     Vector xnew(Esize);                  // vector for new densities

     while (l2-l1 > 0.0001)
      {
       lmid= 0.5*(l2+l1);
 
	       for (unsigned k = 0; k < Esize ; k++)  
		 {
		   if ( x[k]+move < x[k]*sqrt(-dc[k]/lmid)) 
		     xnew[k] = x[k]+move;
		   else
		     xnew[k] = x[k]*sqrt(-dc[k]/lmid);
		 }
	       for (unsigned k = 0; k < Esize ; k++)  
		 {
		   if ( 1.0 < xnew[k]  )               // max. value of density equal 1	
		     xnew[k]= 1.0 ;
		 } 
	       for (unsigned k = 0; k < Esize ; k++)  
		 {
		   if (x[k]-move > xnew[k]  )     
		     xnew[k] = x[k]-move ;
		 }
	       for (unsigned k = 0; k < Esize ; k++)  
		 {
		   if (min_x > xnew[k] )              // min. value of density equal min_x	
		     xnew[k] = min_x ;
		 }

	       for (unsigned k = 0; k < num ; k++)
		 {
		   xnew[x_fix[k]]=1;                  // set densities of fixed elements to 1
		 }
	       
	       float_type sum=0; 
	       for (unsigned k = 0; k < Esize ; k++)  
		 {
		   sum=sum + xnew[k];                 // sum of all new densities
		 }
	   if ( sum - total_x > 0)                    // difference between old volume and new volume
	     l1 = lmid;
	   else
	     l2 = lmid;
      }
     	mtl::copy(xnew, x);
	   return 0;
}
//=================================Filter==================================================
int  BendsoeObject::check(Vector x, Vector &dc, vector<vector<unsigned> > &neighbours, vector<unsigned> ele_neighbours)
{
double fac, sum, sum1;
  Vector dcn(Esize); 
  Vector mdc(Esize);                  // copy of dc

  Dense_Vector Coord1(dim); 
  Dense_Vector Coord2(dim); 
  Dense_Vector Coordvector(dim);

    for (unsigned i=0; i<Esize; i++)
      {
        mdc[i]=dc[i];
        dcn[i] =0;
      }

     if (loop==1)
       {
        //-----evaluate the lenght of a element (only once)--(coordinates 2th node minus coord. 1th node)-------
        mtl::copy(Elements[0]->GetNodeIter(0)->GetCoords(),Coord1);
        mtl::copy(Elements[0]->GetNodeIter(1)->GetCoords(),Coord2);
        rmin=1.5*(Coord2[0]-Coord1[0])  ;                   //rmin=radius of filtering=filtering in a fixed neighborhood 
       //----------------------------------------------------------------------------------------------------
       for (unsigned i=0; i<Esize; i++)
        {
	 for (unsigned k=0; k<Esize; k++)
	    {    
        	mtl::copy(Elements[i]->GetNodeIter(0)->GetCoords(),Coord1);
        	mtl::copy(Elements[k]->GetNodeIter(0)->GetCoords(),Coord2);
        	for (unsigned d=0; d<dim; d++)
	           {
		     Coordvector[d]=Coord2[d]-Coord1[d];
		   }
	       fac = rmin - two_norm(Coordvector);
	       if(fac>0)
		 {
		 ele_neighbours.push_back(k);  //fill vector ele_neighbours with the number of the neighbourood-elementsy
		 }
	    }
	 neighbours[i]=ele_neighbours;        //copy vector ele_neighbours in neighbours[i]
	 ele_neighbours.clear();
	}
       }

   for (unsigned i=0; i<Esize; i++)
    {
	sum=0.0;
	for (unsigned k=0; k < neighbours[i].size(); k++)
	  {    
        	mtl::copy(Elements[i]->GetNodeIter(0)->GetCoords(),Coord1);
        	mtl::copy(Elements[neighbours[i][k]]->GetNodeIter(0)->GetCoords(),Coord2);
        	for (unsigned d=0; d<dim; d++)
	           {
		     Coordvector[d]=Coord2[d]-Coord1[d];
		   }
	       fac = rmin - two_norm(Coordvector);
	       sum1=max(0.0,fac);
	       sum= sum +sum1 ;
	       dcn[i]=dcn[i] +sum1*x[neighbours[i][k]]*mdc[neighbours[i][k]];	       
	  }
	  dc[i]=dcn[i]/(x[i]*sum);
    } 
    return 0;    
}

//===============================Print file to plot====================================================
void  BendsoeObject::PrintFile( Vector x, unsigned Nnr)
{
  vector<StructNode>::iterator nodeit;
  PtrVector<StructElement*>::iterator eleit;
  Vector Nodenumber(Elements.size());       // node number of the 1.node of the 1.elment ...

  unsigned i=0;   
     for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit,i++ ) {
	     nodeit = Nodes.begin();
	     while (nodeit != (*eleit)->GetNodeIter(0)) {
	       ++nodeit;
	     }
	   Nodenumber[i]= distance(Nodes.begin(),nodeit);
	}

      char buffer[80];      
      cout << "Write output file ";
      sprintf(buffer, "%04d.dat", loop);
      puts(buffer);     
      string fname = buffer;
      fstream FS(fname.c_str(), ios::out | ios::trunc);     
      FS << "title = \"Bendsoe-Optimierung\" "<< endl;

      FS << "variables = \"x\", \"y\" ";
      if (dim==3) {
        FS <<",\"z\", \"Dichte\"  " << endl;
        }
      else {
        FS <<", \"Dichte\"  " << endl;
        }

      FS << "zone n=" << Nodes.size() << ", e=" << Elements.size() << ", f=fepoint " ;

      if (dim==3 && Nnr !=10) { 
	  FS <<", et=BRICK " << endl;
        }
      if (dim==3 && Nnr ==10) { 
	  FS <<", et=TETRAHEDRON " << endl;
        }
      if (dim==2 && Nnr==6)  { 
           FS << ", et=TRIANGLE " << endl;
      }
      if (dim==2 && Nnr!=6)  { 
           FS << ", et=QUADRILATERAL " << endl;
      }

      unsigned j=0;
      for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit, j++ )
	{
	  FS << nodeit->Cx << "\t" << nodeit->Cy << "\t";
      if (dim==3) { 
	  FS  << nodeit->Cz << "\t"; 
        }
      else  { 
           FS  << "\t" ;
      }
            unsigned i=0;
	    while(Nodenumber[i] != j  && i < Elements.size()){
		i++;
	      }
	    if(i < Elements.size()){
		FS << x[i] << endl;
	      }	  
	    if(i == Elements.size()){	    
		FS << 0.01 << endl;   
	      }
	}
      FS << endl;
      unsigned k=0;
      unsigned dim2=4;
      if (dim==3 && Nnr!=10){
	dim2=8;                 // element is a brick, 8 nodes for geometry
      }
      if (dim==3 && Nnr==10){
	dim2=4;                 // element (Solid187) is a tetrahedron, 4 nodes for geometry
      }
      if (dim==2 && Nnr!=6) {
	dim2=4;                 // element is quadrilateral, 4  nodes for geometry
      }
      if (dim==2 && Nnr==6) {
	dim2=3;                 // element (Plane2) is triangle, 3 nodes for geometry
      }

      for ( eleit = Elements.begin(); eleit != Elements.end(); ++eleit, k++ ){
	      for (unsigned i=0; i<dim2 ; i++) {                              // brick i < 8 or quadrilateral i<4
		  nodeit = Nodes.begin();
		  while (nodeit != (*eleit)->GetNodeIter(i)){
		      ++nodeit;
		    }
		  FS << distance(Nodes.begin(),nodeit)+1 << "\t" ;
		}
	      FS << endl;
	}
      FS << endl;
      FS.close();     
}
   
//====================================Mainprogram Topopt===============================================
int  BendsoeObject::topopt(EnvelopeMatrix &K, dense1D<unsigned> Envelope)
{
  Esize=Elements.size();
  unsigned Nnr=Elements[0]->GetNodeCount();    //number of nodes of 0.element (all elements are the same type)
  int Ndof=Elements[1]->GetEMSize();           //number of DOF of 0.element
  int Ndofn= Ndof/int(Nnr);                    //number of DOF of a node
  dim = Elements[0]->GetDofSet().ElementDimension();  //dimension of the optimization
  unsigned filter;
  double volfrac;
  double limit;

     //---------------------------------------Inputs--------------------------------------------------------
      cout << endl;
      cout << " with filtering the sensitivies write 1, without write 0: ";   
      cin >> filter;
      while ( filter != 0 && filter != 1)
        {
          cout << " number invalid, write 0 or 1: ";
          cin >> filter;
        }
      cout << endl;

      cout << " volumenfraction, value between 0 and 1 (example: 0.5): ";   
      cin >> volfrac;
      while ( volfrac <= 0 || volfrac >= 1)
        {
          cout << " value invalid (example: 0.3): ";
          cin >> volfrac;
        }
      cout << endl;

      cout << " stop criterion (max. density change), value between 0 and 0.2 (usual: 0.05): ";   
      cin >> limit;
      while ( limit <= 0 || limit > 0.2)
        {
          cout << " value invalid (example: 0.05): ";
          cin >> limit;
        }
      cout << endl;
     //--------------------------------------------Timer-------------------------------------------------------
     //! Object to evaluate runtimes of different parts    
        boost::timer my_timer1;
	boost::timer my_timer2;  
	vector<double> times(5);
  	my_timer2.restart();
     //-------------------------------------------------------------------------------------------------------

  Dense_Matrix  Ke_fix(Ndof,Ndof);        

  Vector DofSolution(Envelope.size());
  Vector Eold(Esize);
  Vector nuold(Esize);
  Vector rhoold(Esize);

  unsigned num=0;
     for (unsigned i = 0; i < Esize; i++)
       {
	 Eold[i]=Elements[i]->GetMaterialPtr()->Get("E");        //copy E-modul in Eold
	 nuold[i]=Elements[i]->GetMaterialPtr()->Get("nu");      //copy nu in nuold
	 rhoold[i]=Elements[i]->GetMaterialPtr()->Get("rho");    //copy real densities of element in rhoold
	 if (rhoold[i] == 1)
	   {
	     num=num+1;                //num=number of fixed densities
	   }
       }

     ApplyLoads(Nodes, K, DofSolution, Envelope);    //vector of forces F=DofSolution


     penal=3;                         // penalty-factor
     Vector x(Esize);                 // vector of the densities
     Vector xold(Esize);
     Vector dc(Esize);                // vector of the objective function
     double total_x=0.0;              // total_x=sum of all element-densities at beginning

     vector<unsigned> x_fix(num);
     unsigned j=0;
  
    for (unsigned k=0; k < Esize; k++)
      {
	x[k] = volfrac;    
	xold[k] = 0;
	dc[k] = 0;
	if(rhoold[k] == 1)
	  {
	    x[k]=1;                   //all fixed elements have density=1
	    x_fix[j]=k;               //in vector x_fix are the elementnumbers of the fixed elements
	    j++;
	  }
	total_x = total_x + x[k];     //volumen at beginning (=sum of all element-densities
      }

    vector<vector<unsigned> > neighbours(Esize);
    vector<unsigned> ele_neighbours;

    loop=0;                           //start interations counter
    double change = 1.0;              //initialisize stop criterion

    while (change > limit)
      {
	loop = loop + 1;
	mtl::copy(x, xold);

	my_timer1.restart();
	FE( x, Envelope, K, DofSolution, Eold, nuold);
	times[0]+=my_timer1.elapsed();

	double c = 0.0;                //objective function
	double s = 0.0;                //skalar product   
	Vector Ue(Ndof);               //vector of the deformations for one element
	Vector T(Ndof);
	Dense_Vector B(6);

	for(unsigned i=0; i< Esize; i++)
	  {
	    for (int k=0; k<int(Nnr); k++)                    //go through all nodes of an element
	      {
        	mtl::copy(Elements[i]->GetNodeIter(k)->GetDeformations(),B); 
		for (int j=0; j<Ndofn; j++)
		  {
		    Ue[Ndofn*k+j] = B[j];                     //write deformations of the element in Ue 
		  }
 	      } 
	   Elements[i]->GetMaterialPtr()->Set("E",Eold[i]);   //set old E-modul again (to calculate EM)
	   Ke_fix=Elements[i]->CalcEM();                      //calculate Element Matrix
	   mult(Ke_fix,Ue,T);
           s=mtl::dot(T,Ue);
	   c = c + pow(x[i],penal)*s;                         //calculate objective function
	   dc[i] = -penal* pow(x[i],(penal-1)) * s ;          //derivative of objective function
	  }

	if (filter==1)
	  {
          my_timer1.restart();
	  check( x, dc, neighbours, ele_neighbours);         //call function check (filtering)
	  times[3]+=my_timer1.elapsed();
	  }

	my_timer1.restart();
        calc_xnew( x, dc,num, x_fix, total_x);              //call function calc_xnew (update densities)
	times[1]+=my_timer1.elapsed();

	    double v=0, vmax=0, sum1=0;
	    for (unsigned k=0; k < Esize; k++)
	      {
		v = abs(x[k] - xold[k]);                     //change of density (new - old)
		sum1 = sum1 + x[k];                          //calculate sum of all new densities
		  if (v > vmax)
		    {
		      vmax = v;                              //max. change of density
		    }
	      }	   
	    change = vmax;
	    cout <<"iterations: "<< loop << endl;            //plot number of iterations
	    cout <<"Vol.:  "<< sum1/(Esize) << endl;         //plot percentage of volume
	    cout <<"Obj.:  "<< c << endl;                    //plot objective function
	    cout <<"change:  "<< change << endl;             //plot max. change of density

            my_timer1.restart();
	    PrintFile( x, Nnr);  
	    times[4]+=my_timer1.elapsed();
	    cout << endl;  
       }
    times[2]=my_timer2.elapsed();
    cout<<"time Topopt: "<<times[2]<<endl;
    cout<<"time FE: "<<times[0]<<endl;
    cout<<"time filter: "<<times[3]<<endl;
    cout<<"time xnew: "<<times[1]<<endl;
    cout<<"time PrintFile: "<<times[4]<<endl;
    cout<<endl;
    return 0;
}
