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
#ifndef Tst68SteppingAction_H
#define Tst68SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst68SteppingAction : public G4UserSteppingAction {

public:

  Tst68SteppingAction();
  virtual ~Tst68SteppingAction();

public:

  virtual void UserSteppingAction( const G4Step* );
  // The main method to define.

  inline G4double getTotalEdepAllParticles() const;
  // Total deposited energy, in the same event, due to all particles,
  // inside the calorimeter.

  inline G4int getPrimaryParticleId() const;
  inline G4double getPrimaryParticleEnergy() const;
  // Id and Energy of the primary particle.

  void reset();
  // This method should be invoked at the end of the event,
  // to reset to zero the total deposit energy variables.

private:
 
  G4double totalEdepAllParticles;
  G4int    primaryParticleId;
  G4double primaryParticleEnergy;
  G4bool   isFirstStepOfTheEvent;

};


inline G4double Tst68SteppingAction::getTotalEdepAllParticles() const {
  return totalEdepAllParticles;
}


inline G4int Tst68SteppingAction::getPrimaryParticleId() const {
  return primaryParticleId;
}


inline G4double Tst68SteppingAction::getPrimaryParticleEnergy() const {
  return primaryParticleEnergy;
}


#endif
