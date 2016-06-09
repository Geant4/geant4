//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4WeightWindowExaminer.hh,v 1.3 2006/06/29 18:17:02 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
