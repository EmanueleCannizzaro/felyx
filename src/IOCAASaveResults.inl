//-----------------------------------------------------------------------------
// IOCAASaveResults.cc
//
// begin     : Mar 18 2001
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


////
//// Implementation of the Member function SaveModel of class IOAnsysFormat
////
template<class analysis_type>
string IOCAAFormat<analysis_type>::SaveResults ( string Fname_, string DataDir_ ) {

  //// -------------
  //// File handling
  //// -------------

  // Modify File name settings, if necessary
  SetFileNames( Fname_, DataDir_ );

  // Create path of file to be openend
  string path = GetResPath();

  ofstream OUT(path.c_str());

  //// -----------------
  //// Print File Header
  //// -----------------

  OUT<<"! FELyX FE Results "<<endl;
  OUT<<"! File: " << path << endl;
  OUT<<"! ------------------------------------------------------"<<endl;
  OUT<<"! "<<endl<<endl;


  //// -----------------------
  //// Printing Deformations
  //// -----------------------

  OUT << endl;
  OUT << "! Deformations" << endl;
  OUT << "! ----------------------------------" << endl;
  OUT << "! nodenumber, ux, uy, uz, rotx, roty, rotz" << endl;
  OUT << "!   " << Nodes.size() << endl;

  Dense_Vector values(6,0.0);
  Dense_Vector maxvalues(6,0.0);
  dense1D<int> maxvaluenodes(6,0);

  int i = 1;
  for ( typename vector<typename analysis_type::node_type>::const_iterator nit = Nodes.begin(); nit != Nodes.end(); ++nit){
    values = nit->GetDeformations();

    // Print node index
    OUT << setw(6) <<  i;

    // Print deformation values
    for (int j=0; j<6; ++j){
      OUT << setw(12) << values[j];

      // Store max vals
      if ( abs( values[j] ) > abs( maxvalues[j]) ){
        maxvalues[j] = values[j];
        maxvaluenodes[j] = i;
      }
    }
    ++i;
    OUT << endl;
  }

  // Print Max values
  OUT << endl;
  OUT << "! Max Values:" << endl;
  OUT << "! Nodes:  ";
  for (int i=0; i < 6; i++)
    OUT << setw(12) << maxvaluenodes[i];
  OUT << endl;
  OUT << "! Values: ";
  for (int i=0; i < 6; i++)
    OUT << setw(12) << maxvalues[i];
  OUT << endl;


  //// -----------------------
  //// Printing element stresses
  //// -----------------------
  OUT << endl << endl;
  OUT << "! Element stresses (if implemented!)" << endl;
  OUT << "! ----------------------------------" << endl;
  OUT << "! elementnumber, elementtype" << endl;
  OUT << "! nodenumber, sx, sy, sz, sxy, syz, sxz (in general)" << endl;
  OUT << "! OR:" << endl;
  OUT << "! elementnumber, elementtype" << endl;
  OUT << "! nodenumber, s_ax, s_bending_y_max, s_bending_z_max, s_phi-x_max (for beam4 elements)" << endl;
  OUT << "!   " << Elements.size() << endl;

  typename PtrVector<typename analysis_type::element_type*>::const_iterator eleit_begin, eleit;
  eleit_begin = Elements.begin();
  eleit = eleit_begin;
  for (eleit = eleit_begin; eleit !=Elements.end(); ++eleit ){

    if ( (*eleit)->GetStatus() && (*eleit)->Stresses.nrows() > 0 ){

      OUT << setw(6) << distance( eleit_begin, eleit)+1
          << setw(6) << (*eleit)->GetId() << endl;

      for (unsigned n=0; n< (*eleit)->Stresses.nrows(); ++n){
        OUT << setw(6) << n+1;

        for (unsigned v=0; v < (*eleit)->Stresses.ncols(); ++v){
          OUT << setw(15) << (*eleit)->Stresses(n,v);

        }
        OUT << endl;

      }
    }
  }

  //// ---------------------------------------------------------------
  //// Close File and finish Load Function
  //// ---------------------------------------------------------------
  OUT.close();

  return path;

}


