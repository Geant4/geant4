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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"

#include "G4VoxelLimits.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4AffineTransform.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle=
  particleTable->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(100*keV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  fpSourceFileUI = new G4UIcmdWithAString("/fpGun/SourceFile",this);
  fpSourceFileUI->SetGuidance("Select the secondary e- source data file.");

  fpVertexUI     = new G4UIcmdWithAnInteger("/fpGun/PrimariesPerEvent", this);
  fpVertexUI->SetGuidance("Number of primary particles per event.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fpVertexUI;
  delete fParticleGun;
  delete fpSourceFileUI;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  for(G4int i = 0; i < fVertex ; i++) {
    G4ThreeVector Position  = SamplePosition();
    G4double Energy         = SampleSpectrum();
    G4ThreeVector Direction = G4RandomDirection();

    fParticleGun->SetParticleMomentumDirection(Direction);
    fParticleGun->SetParticlePosition(Position);
    fParticleGun->SetParticleEnergy(Energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PrimaryGeneratorAction::SetNewValue(G4UIcommand * command,
                                G4String newValue)
{
  if(command == fpSourceFileUI){
    fSourceFile = newValue;
    LoadSpectrum();
  }

  if(command == fpVertexUI) {
    fVertex = fpVertexUI->GetNewIntValue(newValue);

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PrimaryGeneratorAction::LoadSpectrum() 
{
  fEnergyWeights.clear();
  fEnergies.clear();

  std::ifstream SpectrumFile;
  SpectrumFile.open(fSourceFile);
  G4double x, y;
  
  if(!SpectrumFile){
    std::cout << "Source File: " + fSourceFile + " not found!!!" << std::endl;
    exit(1);
  }

  else {
    while (SpectrumFile >> x >> y) {
      fEnergies.push_back(x);
      fEnergyWeights.push_back(y);
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double PrimaryGeneratorAction::SampleSpectrum() 
{
  G4double Energy = 1 * eV;
  G4int Bins = fEnergyWeights.size();

  while (true) {
    G4double X = ((fEnergies[Bins - 1] - fEnergies[0]) * G4UniformRand()) 
                 + fEnergies[0];
    G4double Y = SampleProbability(X);
    if (G4UniformRand() < Y) { 
      Energy = X * MeV;
      break;
    }
  }

  if (Energy > 1 * MeV) 
    Energy = 0.9999999 * MeV;

  return Energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double PrimaryGeneratorAction::SampleProbability(G4double X)  
{
  G4double Y = 0;

  for(std::size_t i = 0; i < fEnergies.size() - 1; i++) {
    if (fEnergies[i] < X && X < fEnergies[i+1]) {
      G4double M = (fEnergyWeights[i+1] - fEnergyWeights[i]) / 
                   (fEnergies[i+1] - fEnergies[i]);
      G4double B = fEnergyWeights[i] - (M*fEnergies[i]);
      Y = M*X + B;
    }
  }

  return Y;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4ThreeVector PrimaryGeneratorAction::SamplePosition() {
  G4VSolid* solid = G4LogicalVolumeStore::GetInstance()
                    ->GetVolume("PlasmidVolume")->GetSolid();
  G4VoxelLimits voxelLimits;
  G4AffineTransform dontUse;
  G4double testX, testY, testZ;
  G4double xmin, xmax, ymin, ymax, zmin, zmax;
  solid->CalculateExtent(kXAxis, voxelLimits, dontUse, xmin, xmax);
  solid->CalculateExtent(kYAxis, voxelLimits, dontUse, ymin, ymax);
  solid->CalculateExtent(kZAxis, voxelLimits, dontUse, zmin, zmax);

  while (true) { 
    testX = G4RandFlat::shoot(xmin, xmax);
    testY = G4RandFlat::shoot(ymin, ymax);
    testZ = G4RandFlat::shoot(zmin, zmax);
    if ( solid->Inside(G4ThreeVector(testX,testY,testZ)) == kInside ) {
      break;
    }
  }

  return G4ThreeVector(testX, testY, testZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....