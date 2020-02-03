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
/// \file field/field03/src/F03PrimaryGeneratorAction.cc
/// \brief Implementation of the F03PrimaryGeneratorAction class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F03PrimaryGeneratorAction.hh"

#include "F03DetectorConstruction.hh"
#include "F03PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ParticleDefinition*  F03PrimaryGeneratorAction::fgPrimaryParticle = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F03PrimaryGeneratorAction::F03PrimaryGeneratorAction(
                                            F03DetectorConstruction* det)
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(0),
    fDetector(det),
    fGunMessenger(0),
    fRndmFlag("off"),
    fXVertex(0.),
    fYVertex(0.),
    fZVertex(0.),
    fVertexDefined(false)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  // create a messenger for this class
  fGunMessenger = new F03PrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);

  fgPrimaryParticle = particle;

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  fParticleGun->SetParticleEnergy(1.*GeV);

  fZVertex=fDetector->GetAbsorberZpos()-0.5*(fDetector->GetAbsorberThickness());
  fParticleGun->SetParticlePosition(G4ThreeVector(fXVertex,fYVertex,fZVertex));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F03PrimaryGeneratorAction::~F03PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F03PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // this function is called at the begining of event
  //
  fgPrimaryParticle = fParticleGun->GetParticleDefinition();

  G4double x0,y0,z0;
  if (fVertexDefined)
  {
    x0 = fXVertex;
    y0 = fYVertex;
    z0 = fZVertex;
  }
  else
  {
    x0 = 0.;
    y0 = 0.;
    z0 = fDetector->GetAbsorberZpos()-0.5*(fDetector->GetAbsorberThickness());
  }

  G4double r0,phi0;
  if (fRndmFlag == "on")
  {
    r0 = (fDetector->GetAbsorberRadius())*std::sqrt(G4UniformRand());
    phi0 = twopi*G4UniformRand();
    x0 = r0*std::cos(phi0);
    y0 = r0*std::sin(phi0);
  }

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String F03PrimaryGeneratorAction::GetPrimaryName()
{
   return fgPrimaryParticle->GetParticleName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F03PrimaryGeneratorAction::SetXVertex(G4double x)
{
  fVertexDefined = true;
  fXVertex = x;
  G4cout << " X coordinate of the primary vertex = " << fXVertex/mm <<
            " mm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F03PrimaryGeneratorAction::SetYVertex(G4double y)
{
  fVertexDefined = true;
  fYVertex = y;
  G4cout << " Y coordinate of the primary vertex = " << fYVertex/mm <<
            " mm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F03PrimaryGeneratorAction::SetZVertex(G4double z)
{
  fVertexDefined = true;
  fZVertex = z;
  G4cout << " Z coordinate of the primary vertex = " << fZVertex/mm <<
            " mm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
