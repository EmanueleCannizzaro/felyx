/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// PreProcessing.h
//
// begin     : Novembre 8 2002
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

#ifndef Lanczos_h
#define Lanczos_h Lanczos_h

#include <float.h>

using namespace std;

namespace mtl{		// Put classes into namespace fe_base

typedef matrix<float_type, symmetric<lower>, envelope<>, row_major >::type EnvelopeMatrix;
typedef matrix<float_type, rectangle<>, dense<>, column_major>::type Dense_Matrix_Col;
typedef dense1D<float_type> Vector;

// declaration of functions
int ewerte(const Vector&, const Vector&, Vector& , const unsigned&, const unsigned& ) ;
int ev_sort(Vector &, Vector &, const unsigned &, const double &);
int eigenvectors(const Vector &, const unsigned &, const Dense_Matrix_Col &, const Vector &, const Vector &, Dense_Matrix_Col &);
void get_envelope(EnvelopeMatrix , dense1D<unsigned> );
void print_ev(Dense_Matrix_Col );
int lanczos(EnvelopeMatrix&, EnvelopeMatrix&, dense1D<unsigned> , unsigned , bool );

inline int ewerte(const Vector &alpha, const Vector &beta, Vector &ev, const unsigned &nev, const unsigned &n) 
{
  double a,b,an,bn,mu,c1(0.00001),eps;
  unsigned vf=0;
  const unsigned kmax=50;
  while (1.0+c1!=1.0)
    {
      eps=c1;
      c1=c1/2.0;
    }
  eps = DBL_EPSILON ;
  Vector temp(n),q(kmax);
  temp[0]=alpha[0]+beta[1];
  temp[n-1]=alpha[n-1]+beta[n-1];
  for (unsigned k=1;k<n-1;k++)
    temp[k]=alpha[k]+beta[k]+beta[k+1];

  b=mtl::max(temp);
  a=-b;

  for (unsigned nrev=0;nrev<nev;nrev++)
    {
      an=a; bn=b;
      for (unsigned k=1;k<kmax;k++)
	{
	  mu=(an+bn)/2.0;
	  q[0]=mu-alpha[0];
	  if (q[0]==0) {q[0]=-eps; cout << "alarm" << endl;}
	  for (unsigned i=1;i<n;i++)
	    {
	      q[i]=mu-alpha[i]-beta[i]*beta[i]/q[i-1];
	      if (q[i]==0) {q[i]=-eps; cout << "alarm" << endl;}
	    }
	  vf=0;
	  for (unsigned i=0;i<n;i++)
	    {
	      if (q[i]<0) vf++;
	    }
        if (vf>=nrev+1) an=mu; 
        else bn=mu;
	}
    ev[nrev]=(an+bn)/2;
    b=bn;
    }

  return 0;
}


int ev_sort(Vector &ev, Vector &ef, const unsigned &n, const double &mu)
{
  double aux(0);
  unsigned j(0), iev(0);
  do
    {
      aux=-1.0e50;
      for (unsigned i=0; i<n; i++)
	{
	  if (ev[i] < 0 && ev[i]>aux)
	    {
	      aux=ev[i];
	      iev=i;
	    }
	}
      if (aux != -1.0e50)
	{
	  ef[j]=aux;
	  ev[iev]=0.0;
	  j++;
	}
    }
  while (aux != -1.0e50);
  do
    {
      aux=0.0;
      for (unsigned i=0; i<n; i++)
	{
	  if (ev[i] > aux)
	    {
	      aux=ev[i];
	      iev=i; 
	    }
	}
      if (aux != 0)
	{
	  ef[j]=aux;
	  ev[iev]=0.0;
	  j++;
	}
    }
  while (aux != 0);
  mtl::copy(ef,ev);
  for (unsigned i=0; i<n; i++)
    {
      ef[i]=sqrt(mu+1.0/ev[i])/(2.0*3.14159265358979);
    }
  return 0;
}


int eigenvectors(const Vector &ev, const unsigned &nEV, const Dense_Matrix_Col &Q, const Vector &alpha, const Vector &beta, Dense_Matrix_Col &X)
{
  unsigned m=Q.nrows(), n=nEV;
  unsigned jmax=alpha.size(), iter;
  Vector al(jmax,0), be(jmax,0), ga(jmax,0), de(jmax,0), vert(jmax,0);
  Vector v(jmax,0), v1(jmax,0);
  double aux, s, s1, s2, u, u1, u2, dif, vn;

  for (unsigned j=0; j<n; j++)
    {
      for (unsigned i=0; i<jmax-1; i++)
	{
	  al[i]=alpha[i]-ev[j];
	  be[i]=beta[i+1];
	  ga[i]=beta[i+1];
	}
      al[jmax-1]=alpha[jmax-1]-ev[j];
      be[jmax-1]=0;
      for (unsigned i=0; i<jmax-1; i++)
	{
        s1=abs(al[i])+abs(be[i]);
        s2=abs(ga[i])+abs(al[i+1])+abs(be[i+1]);
        if (abs(al[i])/s1 >= abs(ga[i])/s2)
	  {
            de[i]=0.0;
            u=ga[i];
            u1=al[i+1];
            u2=be[i+1];
            vert[i]=0;
	  }
        else
	  {
            u=al[i];
	    u1=be[i];
	    u2=0;
	    al[i]=ga[i];
            be[i]=al[i+1];
	    de[i]=be[i+1];
	    vert[i]=1;
	  }
        ga[i]=u/al[i];
	al[i+1]=u1-ga[i]*be[i];
        be[i+1]=u2-ga[i]*de[i];
	}
      set_value(v,0.0);
      set_value(v1,1.0);
      iter=0;
      dif=1.0;
      while (dif >= 1.0e-8)
	{
	  v1[jmax-1]=-v1[jmax-1]/al[jmax-1];
	  s=v1[jmax-1];
	  v1[jmax-2]=-(v1[jmax-2]+be[jmax-2]*v1[jmax-1])/al[jmax-2];
	  if (abs(v1[jmax-2]) > abs(s)) s=v1[jmax-2];
	  for (int i=jmax-3; i>=0; i--)
	    {
	      v1[i]=-(v1[i]+be[i]*v1[i+1]+de[i]*v1[i+2])/al[i];
	      if(abs(v1[i]) > abs(s)) s=v1[i];
	    }
	  dif=0; vn=0;
	  for (unsigned i=0; i<jmax; i++)
	    {
	      aux=v[i];
	      v[i]=v1[i]/s;
	      v1[i]=v[i];
	      vn=vn+v[i]*v[i];
	      dif=std::max(dif,abs(aux-v[i]));
	    }
	  if (dif >= 1.0e-8)
	    {
	      if (iter > 4) cout << "Fehler! Keine Konvergenz." << endl; 
	      for (unsigned i=0; i<jmax-1; i++)
		{
		  if (vert[i] == 0) u=v1[i+1];
		  else 
		    {
		      u=v1[i];
		      v1[i]=v1[i+1];
		    }
		  v1[i+1]=u-ga[i]*v1[i];
		}
	      iter++;
	    }
	}
      vn=sqrt(vn);
      scale(v,1.0/vn);
        for (unsigned i=0; i<m; i++)
	  {
	    s=0;
	    for (unsigned l=0; l<jmax; l++)
	      {
		s=s+v[l]*Q(i,l);
	      }
	    X(i,j)=s;
	  }
    }

  return 0;
}

// void max_norm(const Vector &alpha, const Vector &beta, Vector &norm, const unsigned k, double temp)
// {
//   if (k==0) norm[k]=abs(alpha[k]);
//   if (k==1) 
//     {
//       norm[k]=norm[k-1]+abs(beta[1]);
//       temp=abs(beta[1])+abs(alpha[1]);
//       norm[k]=std::max(norm[k],temp);
//     }
//   if (k>=2)
//     {
//       temp += abs(beta[k]);
//       norm[k]=std::max(norm[k],temp);
//       temp=abs(beta[k])+abs(alpha[k]);
//       norm[k]=std::max(norm[k],temp);
//     }
// }

void get_envelope(EnvelopeMatrix K, dense1D<unsigned> Envelope)
{
  EnvelopeMatrix::iterator 		matrixit;
  EnvelopeMatrix::OneD::iterator 	oneDit;
  unsigned j=0;  
  for(matrixit = K.begin(); matrixit != K.end(); ++matrixit, ++j)
    { 
      oneDit = (*matrixit).begin();
      Envelope[j] = j-oneDit.index()+1;
    }
}

  void print_ev(Dense_Matrix_Col ev){
   string fname = "eigenvectors.txt";
   fstream FS(fname.c_str(), ios::out | ios::trunc);
   for (unsigned j=0;j<ev.nrows();j++)
     {
       FS << "[" ;
       for (unsigned i=0;i<ev.ncols();i++)
	 {
	   if (i>0)
	     FS << ", " << ev(j,i) ;
	   else
	     FS << ev(j,i) ;
	 }
       FS << "];" ;
     }
   FS << endl;
   FS.close();
   cout << "----> wrote eigenvectors to eigenvectors.txt" << endl;
  }

inline int lanczos(EnvelopeMatrix &K, EnvelopeMatrix &M, dense1D<unsigned> Envelope, unsigned nEV, bool evectors, std::vector<double> &EfRes)
{
  //symmetric_tag detlef;
  //twod_copy(M,M2,detlef);

  unsigned jmax=30;
  unsigned m=K.nrows();
  unsigned nconv(0);
  double mu=0.0;
  double tol=1.0e-8;
  // double temp_norm(0.0);
  int err_msg = 0;

  //  cout << " ---> doing Lanczos... " << endl;

  Vector r(m,0);               // vectors of size m
  for (unsigned i=0; i<m; i++)
    r[i]=1;
  Vector u(m,0);
  Vector temp(m,0);
  Vector h(m,0);

  Vector h1(jmax,0);           // vectors of size jmax(+1)
  // Vector norm(jmax);
  Vector alpha(jmax);
  Vector beta(jmax+1);

  Vector ev(nEV);              // vectors of size nEV
  Vector ev_old(nEV);
  Vector ef(nEV);

  Dense_Matrix_Col  Q(m,jmax); // dense matrices
  Dense_Matrix_Col  X(m,nEV);

  row_tag heini;               // tags
  column_tag rosamunde;
  dense_tag udo;

  mult(M,r,h);
  beta[0] = sqrt(dot(h,r));
  //  if (mu!=0.0) 
  //  {
  //   scale(M2,-mu);
  //   add(M2,F);
  // }

  unsigned j=0;
  while (j<jmax && nconv<nEV)
    {
      mtl::copy(r,Q[j]);
      scale(Q[j],1.0/beta[j]);
      scale(h,1.0/beta[j]);
      mtl::copy(h,u);
      if (j==0) sky_decomposition(K,Envelope);
      sky_backsubstitution(K,u,Envelope);
      mtl::copy(u,r);
      if (j>0) 
	{copy(Q[j-1], temp);
	scale(temp,-beta[j]);
	add(temp,r);    }
      alpha[j] = dot(h,r);
      copy(Q[j], temp);
      scale(temp,-alpha[j]);
      add(temp,r);
 
      if (j>1)                 // reorthogonalize
	{
	  mtl::set_value(u,0.0);
	  mult(M,r,u);
	  mtl::set_value(h1,0.0);
	  rect_mult(trans(Q),u,h1,heini,udo);
	  mtl::set_value(temp,0.0);   
	  rect_mult(Q,h1,temp,rosamunde,udo);
	  scale(temp, -1.0);
	  add(temp,r); 
	}

      mtl::set_value(h,0.0);
      mult(M, r, h);
      beta[j+1] = sqrt(dot(h,r));

      //  max_norm(alpha, beta, norm, j, temp_norm);
      if (j>=nEV-1 && j>1)  
      	{
	  mtl::copy(ev,ev_old);
	  ewerte(alpha,beta,ev,nEV,j);
	  nconv=0;
	  for (unsigned i=0; i<nEV; i++)
	    {
	      if ((ev[i]-ev_old[i])/ev[i] < tol) nconv++;
	    }
	}
      j++;
      // cout << " ---> lanczos " << j << " done."  << endl;
    }

  //  cout << "     ------>" << j << " Lanczos steps done." << endl;
  ev_sort(ev,ef,nEV,mu);
  if (evectors) eigenvectors(ev, nEV, Q, alpha, beta, X);
    for (unsigned i=0; i<nEV; i++)
      {
        EfRes[i] = ef[i];
	cout << "---> EF " << i+1 << " = " << ef[i] << " Hz" << endl; 
      }
    if (evectors) print_ev(X);
  return err_msg;
}

} // end of namespace fe_base

#endif
