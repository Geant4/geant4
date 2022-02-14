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
/// \file electromagnetic/TestEm1/include/EmStandardPhysics.hh
/// \brief Definition of the EmStandardPhysics class
//
// $Id: EmStandardPhysics.hh 66586 2012-12-21 10:48:39Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EmStandardPhysics_h
#define EmStandardPhysics_h 1

#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EmStandardPhysics : public G4VPhysicsConstructor
{
 public:
  PrimaryGeneratorAction();
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

  G4ParticleGun* GetParticleGun() { return fParticleGun; };

  void SetOptPhotonPolar();
  void SetOptPhotonPolar(G4double);
  void SetRandomDirection(G4bool val = true);
  G4bool GetPolarized() { return fPolarized; };
  G4double GetPolarization() { return fPolarization; }

 private:
  G4ParticleGun* fParticleGun;
  PrimaryGeneratorMessenger* fGunMessenger;
  G4bool fRandomDirection;
  G4bool fPolarized;
  G4double fPolarization;
};


