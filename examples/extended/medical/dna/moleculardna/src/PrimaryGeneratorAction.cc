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
/// \file PrimaryGeneratorAction.cc
/// \brief Primary Generator for Molecular DNA simulation

#include "PrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"

#include "PrimaryGeneratorSourceGRASCSV.hh"
#include "PrimaryGeneratorMessenger.hh"

// define a Mutex to avoid concurrent reading in multi-thread
#include "G4AutoLock.hh"
namespace { G4Mutex PrimaryGeneratorMutex = G4MUTEX_INITIALIZER; }

// instance of PrimaryGeneratorSource
PrimaryGeneratorSourceGRASCSV* PrimaryGeneratorAction::fPrimarySource = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4AutoLock lock(&PrimaryGeneratorMutex);
  fParticleGun = new G4GeneralParticleSource();

  fFirstEvent = true;
  fGunMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  G4AutoLock lock(&PrimaryGeneratorMutex);
  delete fParticleGun;
  if(fPrimarySource != nullptr)
  {
    delete fPrimarySource;
    fPrimarySource = nullptr;
  }
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(fMyInputFileName != "")
  {
    G4AutoLock lock(&PrimaryGeneratorMutex);
    // Read primaries from file
    // Before first event, instantiate file reader if fInputFileName is not empty
    if(fFirstEvent)
    {
      fPrimarySource = new PrimaryGeneratorSourceGRASCSV(fMyInputFileName);
      fFirstEvent = false;
    }

    auto *fpParticleGun = new G4ParticleGun();

    // Get a new primary particle.
    Primary* primary = fPrimarySource->GetPrimary();
    if(primary != nullptr)
    {
      G4String particleName = primary->GetName();
      G4ThreeVector pos = primary->GetPosition();
      G4ThreeVector momdir = primary->GetMomentumDirection();
      G4double energy = primary->GetEnergy();
      G4ParticleDefinition* particle = primary->GetParticleDefinition();
      //primary->Print();  // print of the data of the primary particle

      fpParticleGun->SetParticleDefinition( particle );
      fpParticleGun->SetParticlePosition( pos );
      fpParticleGun->SetParticleEnergy( energy*CLHEP::MeV );
      fpParticleGun->SetParticleMomentumDirection( momdir );
      fpParticleGun->GeneratePrimaryVertex( anEvent );
    }
    else
    {
      // If primary is NULL, the end of file has been reached or the file format is not consistent
      // A Geantino placed at kInfinity will be fired in this case
      G4cout << "WARNING: The phase space source is ended. Maybe you reach the end of file, or your file is broken. A Geantino will be generated." << G4endl;
      G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* particle = pTable->FindParticle( "geantino" );
      fpParticleGun->SetParticleDefinition( particle );
      G4ThreeVector pos = G4ThreeVector(kInfinity,kInfinity,kInfinity);
      fpParticleGun->SetParticlePosition( pos );
      fpParticleGun->SetParticleEnergy( 0 );
      fpParticleGun->SetParticleMomentumDirection( pos );
      fpParticleGun->GeneratePrimaryVertex( anEvent );
    }
  }
  else
  {
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
