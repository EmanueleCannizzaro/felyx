/** -*- mode: c++; c-indent-level: 2; c++-member-init-indent: 8; comment-column: 35; -*-  */
//-----------------------------------------------------------------------------
// LcmObject.h
//
// begin     : Jan 2 2003
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

#ifndef LcmObject_h
#define LcmObject_h LcmObject_h


#include "FelyxObject.h"
#include "LcmAnalysisType.h"
#include "Entrapment.h"
#include "IOLcmBc.h"
#include "CG.h"


using namespace std;
using namespace fe_base;
using namespace mtl;

namespace felyx{

  // ! Lcm Object class 
  class LcmObject:public FelyxObject<LcmAnalysisType>{

  public:
    
    // Constructors
    // ------------    
    //! Default constructor with optional noise argument
    LcmObject( unsigned noise_ = 0 ) :
      FelyxObject<LcmAnalysisType>( noise_ ),
      LcmBcFileInterface(Nodes, NodeCoordSysList, Elements, Materials, Properties, BoundaryConditions, Layers, Laminates, "", "") {
      solved = 0;
      time = 0;
      it = 0;
      conjunctionlength = 0;
      voidBoundConPtr = new LcmBoundCon;
      voidBoundConPtr->set(P,0);
    }
    
    //! File constructor, reading in FE model 
    //! Args: filenam, datadir, noiselevel (optional)
    LcmObject( string filenam_, string LcmBCfname_, string datadir_, unsigned noise_ = 0, double output_ = 0, double entrapment_ = 0, double overfill_ = 1):
      FelyxObject<LcmAnalysisType>( filenam_, datadir_, noise_ ), LcmBcFileInterface(Nodes, NodeCoordSysList, Elements, Materials, Properties, BoundaryConditions, Layers, Laminates, LcmBCfname_, datadir_ ) {
      if (LcmBCfname_ != "") LoadLcmBcFile(LcmBCfname_, datadir_);
      output = output_;
      entrapment = entrapment_;
      overfill = overfill_;
      solved = 0;
      time = 0;
      it = 0;
      conjunctionlength = 0;
      voidBoundConPtr = new LcmBoundCon;
      voidBoundConPtr->set(P,0);
    }
    ~LcmObject() {
    };


    string LoadLcmBcFile( string, string);

    void NodeNumbering() {
      unsigned i=1;
      for ( nodeit = Nodes.begin(); nodeit != Nodes.end(); ++nodeit, ++i ){
	nodeit->number = i;
      }
    }
    
    void PreProcessing();

    //! Direct solver for pressure distribution.
    void DirectSolver();

    //! Iterative solver for pressure distribution.
    void IterativeSolver();

    //! Detecting nodes that belong to one vent.
    void DetectVents();

    //! Fill nodes where pressure bc is not 0. Set temporary bc where fill is 0.
    void SetNewBC();

    void SetTdepBC();
    
    //! Evaluate neighbour nodes of a node.
    /*! Write neighbour node numbers into a vector.
      f.e. the vector neighbours[4] contains the numbers of the neighbour nodes of node 4. */
    void EvalNeighbours();

    //! Compute timestep and update fill state of each node
    double Comp_dt();

    //! Prints a .dat file containing fill, pressure and flux state at each node.
    void PrintFile();

    //! Constructs new entrapments, tracks entrapments and computes their
    /*! new volume and pressure.
    */
    void UpdateEntrapments();

    void CompareEntrapments(vector<Entrapment*>&, vector<Entrapment*>&);

    //! Sets Nodes[i].entrapment to
    /*! 0 if there is a free way from the node to a vent or
        1, 2, ... if node belongs to a certain air entrapment. Nodes that 
        belong to the same entrapment have same numbers.
        \param neighbours Neighbours list to each node
    */
    int EvalEntrapment();

    //! Called by EvalEntrapment.
    int EntrapCheck(LcmNode*);
    int SetFreeWay(LcmNode*);

    //! Closes a vent if the flow front has completely reached it.
    bool CloseVent();

    double GetFillingTime() {return FillingTime;};
    double GetTotalVolume() {return TotalVolume;};
    double GetFilledVolume() {return FilledVolume;};
    double GetConjunctionLength() {return conjunctionlength;};

    void PostProcessing();

    void PrintGlobalStatus();

  public:
    //DATA MEMBERS
    //------------
    IOLcmBc LcmBcFileInterface;

    unsigned it;
    double entrapment, overfill, output, dt;

    bool solved;
    double time;

    double conjunctionlength;

    //! Filling Time 
    double FillingTime;
    //! Total Volume 
    double TotalVolume;
    //! Filled Volume Fraction
    double FilledVolume;

    LcmBoundCon *voidBoundConPtr;

    vector<LcmNode>::iterator nodeit;
    vector<LcmElement*>::iterator eleit;

    vector<vector<unsigned> > neighbours;
    vector<vector<LcmNode*> > vents;
    vector<Entrapment> Entraps;
    vector<TdepBC> FileBC;

    // permutation vector
    // vector index: original numbering
    // vector entries: reordered numbering
    static vector<unsigned> Perm;
  };
  
}

#endif
