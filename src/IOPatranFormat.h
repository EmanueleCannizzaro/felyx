//-----------------------------------------------------------------------------
// IOPatranFormat.h
//
// begin     : Jan 20 2004
// copyright : (c) 2001 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {okoenig, wintermantel, nzehnder}@imes.mavt.ethz.ch
// www       : www.imes.ethz.ch/st
//
// $Id: IOPatranFormat.h,v 1.1.1.1 2004/12/13 17:55:16 okoenigl Exp $
/*
   This file is part of FELyX (Finite Element Library eXperiment).

la   FELyX is free software; you can redistribute it and/or modify
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

#ifndef IOPatranFormat_h
#define IOPatranFormat_h IOPatranFormat_h

#include <string>
#include <iomanip>
#include "IOFelyx.h"

using namespace std;

namespace fe_base{

  //! Class to store FE-Data in ANSYS-Format
  template<class analysis_type> class IOPatranFormat : public IOFelyx<analysis_type>
  {

  public:

    //! Constructor, initializing references to all FE-Data objects.
    /*! Additionally, filename and datadirectory can optionally be specified */
    IOPatranFormat( vector<typename analysis_type::node_type>&          _Nodes,
                   vector<CoordSys>&            NodeCoordSysList_,
                   PtrVector<typename analysis_type::element_type*>&   Elements_,
                   PtrVector<Material*>&        _Materials,
                   vector<PropertySet>&         _Properties,
                   vector<typename analysis_type::bc_type>&     _BoundaryConditions,
                   vector<Layer>&               Layers_,
                   vector<Laminate>&            Laminates_,
                   map<string, double>&         _PatranParameters,
                   string                       _Fname   = "felyx-3d-mixed-718",
                   string                       _DataDir = "../fedata" );

    //! Virtual Destructor
    ~IOPatranFormat() {cout << " ------------------- iopatranformat destructor" << endl;};

    //! Load FE Model from file
    /*! Return type:        Path of file loaded
        Optional Arguments: Filename, DataDirectory */
    virtual string LoadModel ( string Fname_ ="", string DataDir_ ="" );

  private:

    //! Additional Patran Specific Data format to be handled by reference
    map<string, double>& PatranParameters;

    void readOutFile(string path);
    void readLayupFile(string path);

    class LineReader {
      private:
        string str_;
        int pos_;
      public:
        LineReader(string str) : str_(str) {
          pos_ = 0;
        }

        int Int(int len) {
          int res = atoi(str_.substr(pos_, len).c_str());
          pos_ += len;
          return res;
        }

        double Double(int len) {
          double res = atof(str_.substr(pos_, len).c_str());
          pos_ += len;
          return res;
        }

        void skip(int len) {
          pos_ += len;
        }
    };

    class CoordStore {
      private:
        vector<CoordSys> *coordList;
      public:
        CoordStore(vector<CoordSys> *list) : coordList(list) {
          coordList->clear();
        }

        int add(double th1, double th2, double th3) {
          CoordSys sys(CoordSys::cartesian, CoordSys::deg, CoordSys::Euler312, th1, th2, th3);
          vector<CoordSys>::iterator it = find(coordList->begin(), coordList->end(), sys);
          if(it == coordList->end()) {
            coordList->push_back(sys);
            return coordList->size();
          } else {
            return distance(coordList->begin(), it);
          }
        }
    };

  };


#include "IOPatranFormat.inl"
#include "IOPatranLoadModel.inl"


} // of namespace
#endif

