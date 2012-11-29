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
// $Id$
//
/////////////////////////////////////////////////////////////////////////
//
// EventActionMessenger
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

#include "PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "HistoManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  histo = HistoManager::GetPointer();

  particleGun  = new G4ParticleGun(1);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  G4String pname = "";
  char* path = getenv("PRIMARYBEAM");
  if (path) pname = G4String(path); 
  if(pname != "") {
    particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(pname));
  } else {
    particleGun->SetParticleDefinition(G4Electron::Electron());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double z = -0.49*histo->GetWorldLength();
  G4double xySize = histo->GetBeamSizeXY();
  G4double x = 0.0;
  G4double y = 0.0;
  if(xySize > 0.0) {
    x = xySize*(2.0*G4UniformRand() - 1.0);
    y = xySize*(2.0*G4UniformRand() - 1.0);
  }
  particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
