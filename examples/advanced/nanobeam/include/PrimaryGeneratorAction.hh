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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <CLHEP/Matrix/Matrix.h>

#include "Randomize.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"

#include "DetectorConstruction.hh"

class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  PrimaryGeneratorAction(DetectorConstruction*);    
  ~PrimaryGeneratorAction();
  
  void GeneratePrimaries(G4Event*);

  G4ParticleGun* GetParticleGun() {return fParticleGun;};
  
  void SetEmission (G4int);

  CLHEP::HepMatrix GetMatrix(){return fBeamMatrix;};
  
  G4int fEmission;
 
private:
  
  G4double XYofAngle(G4double);	
  
  G4ParticleGun* fParticleGun;
  
  DetectorConstruction* fDetector;
  
  PrimaryGeneratorMessenger* fGunMessenger;     
  
  G4double fAngleMax;
  
  G4bool fShoot;
  
  // Matrix
  CLHEP::HepMatrix fBeamMatrix;	
};

#endif
