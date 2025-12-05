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


#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4IAEAphspReader.hh"

#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4RunManager.hh"

#include <vector>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(const G4int threads)
  :fThreads(threads)
{
  fIAEAphspReaderName = "";
  fVerbose = 0;
  fMessenger = new PrimaryGeneratorMessenger(this);
  fParticleGun = new G4ParticleGun();
  fParticleGun->SetParticleDefinition(G4Gamma::Definition());

  fCounter = 0;
  fKinE = 50.0*MeV;
  fDE = 0.0;
  fX0 = 0.0;
  fY0 = 0.0;
  fZ0 = 0.0;
  fDX = 0.0;
  fDY = 0.0;
  fDZ = 0.0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if (fVerbose > 0) G4cout << "Destroying PrimaryGeneratorAction" << G4endl;
  delete fParticleGun;
  delete fMessenger;
  
  if (fIAEAphspReader) delete fIAEAphspReader;
  if (fVerbose > 0) G4cout << "PrimaryGeneratorAction destroyed" << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (fIAEAphspReader) {
    fIAEAphspReader->GeneratePrimaryVertex(anEvent);
  }
  
  else {

    fCounter++ ;

    // Simulation of beam kinetic energy
    G4double kinEnergy = fKinE;
    if(fDE > 0.0)
      kinEnergy = G4RandFlat::shoot(fKinE-fDE/2., fKinE+fDE/2.);

    fParticleGun->SetParticleEnergy(kinEnergy);

    // Simulation of beam position
    G4double x = fX0;
    G4double y = fY0;
    G4double z = fZ0;
    if (fDX > 0.0)
      x = G4RandFlat::shoot(fX0-fDX/2., fX0+fDX/2.);
    if (fDY > 0.0)
      y = G4RandFlat::shoot(fY0-fDY/2., fY0+fDY/2.);
    if (fDZ > 0.0)
      z = G4RandFlat::shoot(fZ0-fDZ/2., fZ0+fDZ/2.);

    fParticleGun->SetParticlePosition( G4ThreeVector(x,y,z) );

    // Simulation of beam direction
    G4double ux = 0.0;
    G4double uy = 0.0;
    G4double uz = 1.0;

    // Beam particles are randomly going upwards or downwards
    // This is done in order to let G4IAEAphspReader show how n_stat works
    if(G4UniformRand() < 0.5)
      uz = -uz;

    fParticleGun->SetParticleMomentumDirection( G4ThreeVector(ux,uy,uz) );

    if(fVerbose > 1) {
      G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
      G4String particleName = particle->GetParticleName();
      G4cout << G4endl
	     << "Event # " << fCounter
	     << "  ParticleGun vertex:  "
	     << "ParticleName = " << particleName
	     << "   PDGcode = " << particle->GetPDGEncoding()
	     << G4endl;
      G4cout << std::setprecision(6)
	     << "\t\t KinEnergy (MeV) = " << kinEnergy/MeV
	     << "   weight = " << fParticleGun->GetParticleWeight()
	     << G4endl
	     << "\t\t x (cm) = "  << x/cm
	     << "  y (cm) = "   << y/cm
	     << "  z (cm) = "   << z/cm
	     << "    ux = " << ux << "  uy = " << uy << "  uz = " << uz
	     << G4endl;
    }

    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetIAEAphspReader(const G4String filename)
{
  fIAEAphspReaderName = filename;
  fIAEAphspReader = new G4IAEAphspReader(filename, fThreads);
}
