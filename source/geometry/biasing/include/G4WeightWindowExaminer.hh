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
// $Id: G4WeightWindowExaminer.hh,v 1.1 2003-08-19 15:17:26 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4WeightWindowExaminer
//
// Class description:
//

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
    // initialisation and construct G4ImportanceFinder

  virtual ~G4WeightWindowExaminer();
    // delete G4ImportanceFinder

  virtual G4Nsplit_Weight Examine(G4double w, G4double energy) const; 
    // Get  G4Nsplit_Weight for a given mother track weight.
  

private:

  G4WeightWindowExaminer(const G4WeightWindowExaminer &);
  G4WeightWindowExaminer &operator=(const G4WeightWindowExaminer &);

private:

  const G4VWeightWindowAlgorithm &fWWalgorithm;
  const G4VParallelStepper &fPStepper;
  const G4VWeightWindowStore &fWWStore;
};


#endif
