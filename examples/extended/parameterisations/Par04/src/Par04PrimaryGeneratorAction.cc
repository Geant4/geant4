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
#include "Par04PrimaryGeneratorAction.hh"
#include <CLHEP/Units/SystemOfUnits.h>       // for GeV
#include <G4String.hh>                       // for G4String
#include <G4ThreeVector.hh>                  // for G4ThreeVector
#include <G4Types.hh>                        // for G4int
#include <G4VUserPrimaryGeneratorAction.hh>  // for G4VUserPrimaryGeneratorA...
#include <string>                            // for basic_string
#include "G4Event.hh"                        // for G4Event
#include "G4ParticleGun.hh"                  // for G4ParticleGun
#include "G4ParticleTable.hh"                // for G4ParticleTable
#include "G4SystemOfUnits.hh"                // for GeV
#include "Par04EventInformation.hh"          // for Par04EventInformation
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04PrimaryGeneratorAction::Par04PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun     = new G4ParticleGun(n_particle);
  // Default particle properties
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 1., 0.));
  fParticleGun->SetParticleEnergy(10. * GeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04PrimaryGeneratorAction::~Par04PrimaryGeneratorAction() { delete fParticleGun; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04PrimaryGeneratorAction::GeneratePrimaries(G4Event* aEvent)
{
  fParticleGun->GeneratePrimaryVertex(aEvent);
  aEvent->SetUserInformation(new Par04EventInformation());
}
