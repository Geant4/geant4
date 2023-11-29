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
//
// ----------------------------------------------------------------------
// Class G4CellScoreValues
//
// Class description:
// This struct holds sums of scores for standard scoring.
//
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4CellScoreValues_hh
#define G4CellScoreValues_hh G4CellScoreValues_hh

#include "globals.hh"

class G4CellScoreValues
{
 public:
  G4double fSumSL = 0;
  G4double fSumSLW = 0;
  G4double fSumSLW_v = 0;
  G4double fSumSLWE = 0;
  G4double fSumSLWE_v = 0;
  G4int fSumTracksEntering = 0;
  G4int fSumPopulation = 0;
  G4int fSumCollisions = 0;
  G4double fSumCollisionsWeight = 0;
  G4double fNumberWeightedEnergy = 0;
  G4double fFluxWeightedEnergy = 0;
  G4double fAverageTrackWeight = 0;
  G4double fImportance = 0;
};

#endif
