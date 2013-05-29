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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <CLHEP/Matrix/Matrix.h>

#include "globals.hh"
#include "Randomize.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(DetectorConstruction*);    
  ~PrimaryGeneratorAction();
  
public:
  
  void GeneratePrimaries(G4Event*);
  G4ParticleGun* GetParticleGun() {return particleGun;};
  void SetEmission (G4int);

  CLHEP::HepMatrix GetMatrix(){return beamMatrix;};
  
  G4int emission;
 
private:
  
  G4double XYofAngle(G4double);	
  
  G4ParticleGun* particleGun;
  DetectorConstruction* detector;
  PrimaryGeneratorMessenger* gunMessenger;     
  G4double angleMax;
  G4bool shoot;
  
   // Matrix
  CLHEP::HepMatrix beamMatrix;	
};

#endif
