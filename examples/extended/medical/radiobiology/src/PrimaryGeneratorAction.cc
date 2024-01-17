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
//
/// \file radiobiology/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the RadioBio::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
  : G4VUserPrimaryGeneratorAction(), fDetector(det)
{
  fGPS = new G4GeneralParticleSource();

  G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fGPS->SetParticleDefinition(particle);

  // Hard code the primary direction along X
  //    and rotate the source
  fGPS->GetCurrentSource()->GetPosDist()->SetPosRot1(G4ThreeVector(0, 0, 1));
  fGPS->GetCurrentSource()->GetPosDist()->SetPosRot2(G4ThreeVector(0, 1, 0));
  fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(1, 0, 0));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Hard code position on the side and rotatio on the corrrect axis
  fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(
    G4ThreeVector(-fDetector->GetSizeX() / 2., 0., 0.));

  fGPS->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio