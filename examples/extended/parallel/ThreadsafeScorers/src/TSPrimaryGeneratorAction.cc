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
/// \file parallel/ThreadsafeScorers/src/TSPrimaryGeneratorAction.cc
/// \brief Implementation of the TSPrimaryGeneratorAction class
//
//
//
//
/// Simple PrimaryGeneratorAction that produces a -Z surface flux of 1 MeV
///     neutrons into the world
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "TSPrimaryGeneratorAction.hh"
#include "TSDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4TiMemory.hh"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSPrimaryGeneratorAction::TSPrimaryGeneratorAction()
{
    fGun = new G4ParticleGun(1);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="neutron");

    fGun->SetParticleDefinition(particle);
    fGun->SetParticleEnergy(1.*MeV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSPrimaryGeneratorAction::~TSPrimaryGeneratorAction()
{
    delete fGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    TIMEMORY_AUTO_TIMER("");
    static TSDetectorConstruction* detector
            = TSDetectorConstruction::Instance();

    G4ThreeVector dir(0.,0., 1.);
    G4ThreeVector pos(detector->GetWorldDimensions().x()*(G4UniformRand()-0.5),
                      detector->GetWorldDimensions().y()*(G4UniformRand()-0.5),
                      -0.5*detector->GetWorldDimensions().z());

    dir /= dir.mag();
    fGun->SetParticleMomentumDirection(dir);
    fGun->SetParticlePosition(pos);
    fGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

