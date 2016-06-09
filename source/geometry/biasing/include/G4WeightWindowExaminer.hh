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
// $Id: G4WeightWindowExaminer.hh,v 1.2 2003/08/19 15:44:57 dressel Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// ----------------------------------------------------------------------
// Class G4WeightWindowExaminer
//
// Class description:
// 
// Implementation of an weight window examiner according to
// G4VWeightWindowExaminer interface. 
// See also description in G4VWeightWindowExaminer.hh

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4WeightWindowExaminer_hh
#define G4WeightWindowExaminer_hh G4WeightWindowExaminer_hh

#include "G4VWeightWindowExaminer.hh"

class G4VWeightWindowAlgorithm;
class G4VParallelStepper;
class G4VWeightWindowStore;

class G4WeightWindowExaminer: public G4VWeightWindowExaminer
{

public:  // with description

  G4WeightWindowExaminer(const G4VWeightWindowAlgorithm &aWWalg,
			 const G4VParallelStepper &astepper,
			 const G4VWeightWindowStore &wwstore);
    // constructor 

  virtual ~G4WeightWindowExaminer();
    // destructor

  virtual G4Nsplit_Weight Examine(G4double w, G4double energy) const; 
    // evaluate the number of tracks and their weight according
    // to the initial track weight and energy

  

private:

  G4WeightWindowExaminer(const G4WeightWindowExaminer &);
  G4WeightWindowExaminer &operator=(const G4WeightWindowExaminer &);

private:

  const G4VWeightWindowAlgorithm &fWWalgorithm;
  const G4VParallelStepper &fPStepper;
  const G4VWeightWindowStore &fWWStore;
};


#endif
