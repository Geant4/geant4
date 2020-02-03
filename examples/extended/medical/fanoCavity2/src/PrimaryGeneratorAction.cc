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
/// \file medical/fanoCavity2/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
:fParticleGun(0), fDetector(det)
{
   fParticleGun  = new G4ParticleGun(1);
  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(1*MeV);  
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::RunInitialisation(G4double effWall, G4double massR)
{
  //this function is called at beginning of run


  fCavityThickness = fDetector->GetCavityThickness();
  fEffWallThick = effWall;
  fMassWallRatio  = massR;
    
  fNwall = fNcavity = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //      
  G4double Zvertex;
  if (G4UniformRand() < fMassWallRatio) {
    Zvertex = 0.5*fCavityThickness + G4UniformRand()*fEffWallThick;
    if (G4UniformRand() < 0.5) Zvertex = -Zvertex;
    fNwall++;
  } else {
    Zvertex = (G4UniformRand() - 0.5)*fCavityThickness;
    fNcavity++;
  }
  
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., Zvertex));  
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
  analysisManager->FillH1(1,Zvertex);
  analysisManager->FillH1(2,fParticleGun->GetParticleEnergy()); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

