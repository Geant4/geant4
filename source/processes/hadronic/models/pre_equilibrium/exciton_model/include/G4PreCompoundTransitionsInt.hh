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
// GEANT4 Class header file
//
// File name:  G4PreCompoundTransitionInt
//
// Author:  V.Ivantchenko, 25 January 2025
//
// Class Description:
// Model implementation for pre-equilibrium transition inside a nucleus.
// It is an alternative to the default model.
//

#ifndef G4PreCompoundTransitionsInt_h
#define G4PreCompoundTransitionsInt_h 1

// Compute transition probailities:
// TransitionProb1 => probability of transition with  \Delta N = +1
//                    number of excitons will be increased on 2
// TransitionProb2 => probability of transition with  \Delta N = -1
//                    number of excitons will be decreased on 2
// TransitionProb3 => probability of transition with  \Delta N = 0
//                    number of excitons will be the same

#include "G4VPreCompoundTransitions.hh"
#include "globals.hh"

class G4ParticleDefinition;
class G4Fragment;
class G4NuclearLevelData;

class G4PreCompoundTransitionsInt : public G4VPreCompoundTransitions
{
public:

  G4PreCompoundTransitionsInt(G4int verb);

  ~G4PreCompoundTransitionsInt() override = default;

  G4double CalculateProbability(const G4Fragment & aFragment) override;
  
  void PerformTransition(G4Fragment & aFragment) override;
  
  G4PreCompoundTransitionsInt(const G4PreCompoundTransitionsInt&) = delete;
  const G4PreCompoundTransitionsInt& operator=
  (const G4PreCompoundTransitionsInt& right) = delete;
  G4bool operator==(const G4PreCompoundTransitionsInt& right) const = delete;
  G4bool operator!=(const G4PreCompoundTransitionsInt& right) const = delete;

private:

  const G4ParticleDefinition* proton;
  G4NuclearLevelData* fNuclData;

  G4double FermiEnergy;
  G4double r0;  // Nuclear radius
  G4int fVerbose;
};

#endif
