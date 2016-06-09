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
// $Id: G4CellScoreValues.hh,v 1.1 2003/10/03 10:07:07 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
