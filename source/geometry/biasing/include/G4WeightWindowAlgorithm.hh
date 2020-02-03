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
// G4WeightWindowAlgorithm
// 
// Class description:
//
// Implementation of a weight window algorithm. 
// The arguments in the constructor configure the algorithm:
//   - upperLimitFactor: the factor defining the upper weight limit 
//      W_u = upperLimitFactor * W_l (W_l lower weight bound)
//   - survivalFactor: used in calculating the survival weight
//      W_s = survivalFactor * W_l
//   - maxNumberOfSplits: the maximal number of splits allowed
//     to be created in one go, and the reciprocal of the minimal
//     survival probability in case of Russian roulette
//
// In case of upperLimitFactor=survivalFactor=1 the algorithm becomes
// the expected weight algorithm of the importance sampling technique.
//
// See also G4VWeightWindowAlgorithm.

// Author: Michael Dressel (CERN), 2003
// ----------------------------------------------------------------------
#ifndef G4WEIGHTWINDOWALGORITHM_HH
#define G4WEIGHTWINDOWALGORITHM_HH 1

#include "G4Nsplit_Weight.hh"
#include "G4VWeightWindowAlgorithm.hh"

class G4WeightWindowAlgorithm : public G4VWeightWindowAlgorithm
{
  public:  // with description
  
    G4WeightWindowAlgorithm(G4double upperLimitFactor = 5,
                            G4double survivalFactor = 3,
                            G4int maxNumberOfSplits = 5);
  
    virtual ~G4WeightWindowAlgorithm();

    virtual G4Nsplit_Weight Calculate(G4double init_w,
                                      G4double lowerWeightBound) const;
      // calculate number of tracks and their weight according
      // to the initial track weight and the lower energy bound

  private:

    G4double fUpperLimitFactor = 0.0;
    G4double fSurvivalFactor = 0.0;
    G4int fMaxNumberOfSplits = 0;
};

#endif
