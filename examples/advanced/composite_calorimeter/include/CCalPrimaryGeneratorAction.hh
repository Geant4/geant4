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
///////////////////////////////////////////////////////////////////////////////
// File: CCalPrimaryGeneratorAction.hh
// Description: Sets up particle beam
//
//     By default 1 pi+  is shot from (0,0,0)
//     in (1,1,0.1) direction at 100 GeV
//     Use /gun/... commands to modify energy,origin,direction at run time.
//     or/and
//         /CCal/generator/random true/false to have random direction
//         /CCal/generator/scan   true/false to scan in eta/phi
//     Use 
//         /CCal/generator/minEnergy
//         /CCal/generator/maxEnergy
//         /CCal/generator/minPhi 
//         /CCal/generator/maxPhi 
//         /CCal/generator/minEta
//         /CCal/generator/maxEta 
//     to set the range in energy and direction of particles shot at random.
//     Use 
//         /CCal/generator/stepsPhi
//         /CCal/generator/stepsEta
//     to set number of steps in Phi and Eta for the scan
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalPrimaryGeneratorAction_h
#define CCalPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4Event.hh"
#include "G4VPrimaryGenerator.hh"

#include "CCalPrimaryGeneratorMessenger.hh"

enum generatorInputType {singleFixed, singleRandom, singleScan};

class CCalPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  CCalPrimaryGeneratorAction();
  ~CCalPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);
  
public:
  void SetVerboseLevel(G4int val);
  void SetRandom(G4String val);
  void SetScan(G4String val);
  void SetMinimumEnergy(G4double p);
  void SetMaximumEnergy(G4double p);
  void SetMinimumPhi(G4double p);
  void SetMaximumPhi(G4double p);
  void SetStepsPhi(G4int val);
  void SetMinimumEta(G4double p);
  void SetMaximumEta(G4double p);
  void SetStepsEta(G4int val);
  void SetGunPosition(const G4ThreeVector & pos) const;
  void SetRunNo(G4int val);

public:    
  G4ThreeVector GetParticlePosition() {return particleGun->GetParticlePosition();}
  G4double GetParticleEnergy() {return particleGun->GetParticleEnergy();}

private:
  CCalPrimaryGeneratorMessenger* gunMessenger;
  G4ParticleGun* particleGun;
  generatorInputType generatorInput;

  G4int verboseLevel;
  G4int n_particle;
  G4String particleName;
  G4double particleEnergy;
  G4ThreeVector particlePosition;
  G4ThreeVector particleDir;

  G4double energyMin,energyMax;
  G4double etaMin,etaMax;
  G4double phiMin,phiMax;
  G4int etaSteps,phiSteps;

  G4int isInitialized;
  G4double etaValue, phiValue;
  G4int scanSteps;

private:
  void initialize();
  void print(G4int val);

};

#endif
