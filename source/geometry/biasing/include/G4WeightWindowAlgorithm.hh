//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4WeightWindowAlgorithm.hh,v 1.3 2002-08-29 15:30:50 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4WeightWindowAlgorithm
// 
// Class description:
// see also G4VWeightWindowAlgorithm.
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4WeightWindowAlgorithm_hh
#define G4WeightWindowAlgorithm_hh G4WeightWindowAlgorithm_hh 

#include "G4Nsplit_Weight.hh"
#include "G4VWeightWindowAlgorithm.hh"

class G4WeightWindowAlgorithm : public G4VWeightWindowAlgorithm
{

public:  // with description
  
  G4WeightWindowAlgorithm();
  
  ~G4WeightWindowAlgorithm(){}

  void SetUpperLimit(G4double Upper);
    // set upper limiting factor for window
    // - Upper is the maximum factor by which the weight of a particle may
    //   be higher than it should be according to the importance
  
  void SetLowerLimit(G4double Lower);
    // set lower limiting facotr for window
    // - Lower is the minimal factor by which the wight may be
    //   lower than it should be according to the importance



  G4Nsplit_Weight Calculate(G4double init_w, 
			    G4double importance) const;
    // calculate the number of tracks and their weight according 
    // to the upper and lower limmiting factors of the window
    // and the initial weight
    // - init_w is the initial weight of the particle
    // - importance is the importance of the cell

private:
  G4double fUpper;
  G4double fLower;
};

#endif
