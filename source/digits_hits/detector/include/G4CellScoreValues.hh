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
// $Id: G4CellScoreValues.hh 67992 2013-03-13 10:59:57Z gcosmo $
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

class G4CellScoreValues {
public:
  G4CellScoreValues() :
    fSumSL(0),
    fSumSLW(0),
    fSumSLW_v(0),
    fSumSLWE(0),
    fSumSLWE_v(0),
    fSumTracksEntering(0),
    fSumPopulation(0),
    fSumCollisions(0),
    fSumCollisionsWeight(0),
    fNumberWeightedEnergy(0),
    fFluxWeightedEnergy(0),
    fAverageTrackWeight(0),
    fImportance(0)
  {}
    
  G4double fSumSL;
  G4double fSumSLW;
  G4double fSumSLW_v;
  G4double fSumSLWE;
  G4double fSumSLWE_v;
  G4int fSumTracksEntering;
  G4int fSumPopulation;
  G4int fSumCollisions;
  G4double fSumCollisionsWeight;
  G4double fNumberWeightedEnergy;
  G4double fFluxWeightedEnergy;
  G4double fAverageTrackWeight;
  G4double fImportance;
};

#endif
