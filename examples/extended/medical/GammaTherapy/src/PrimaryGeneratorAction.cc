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
// $Id: PrimaryGeneratorAction.cc 103662 2017-04-20 14:58:33Z gcosmo $
//
/// \file medical/GammaTherapy/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//

//---------------------------------------------------------------------------
//
// ClassName:   PrimaryGeneratorAction
//
// Description: Generate primary beam
//
// Authors: V.Grichine, V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "Randomize.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* pDet):
  fDetector(pDet)
{
  InitializeMe();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorAction::InitializeMe()
{
  fVerbose = fDetector->GetVerbose();
  fMessenger = new PrimaryGeneratorMessenger(this);
  fParticleGun = new G4ParticleGun();
  fCounter = 0;
  fX0 = 0.0;
  fY0 = 0.0;
  fZ0 = 0.0;
  fSigmaX = 1.5*mm;
  fSigmaY = 1.5*mm;
  fSigmaZ = 0.0;
  fSigmaE = 0.0;
  fRMax2  = 2.5*2.5*mm*mm;
  fSigmaTheta = 0.0;
  //  fSigmaTheta = 0.17*degree;
  fMinCosTheta = 2.0;
  SetBeamEnergy(50.0*MeV);
  fPosition  = G4ThreeVector(fX0,fY0,fZ0);
  fDirection = G4ThreeVector(0.0,0.0,1.0);
  fGauss = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fCounter++ ;

  // Simulation of beam position
  G4double x = fX0;
  G4double y = fY0;
  G4double z = fDetector->GetGeneratorPosZ();
  do {
    if(0.0 < fSigmaX) { x = G4RandGauss::shoot(fX0,fSigmaX); }
    if(0.0 < fSigmaY) { y = G4RandGauss::shoot(fY0,fSigmaY); }
  } while (x*x + y*y > fRMax2);

  fPosition = G4ThreeVector(x,y,z);
  fParticleGun->SetParticlePosition(fPosition);

  // Simulation of beam direction
  G4double ux = fDirection.x();
  G4double uy = fDirection.y();
  G4double uz = fDirection.z();

  // Beam particles are uniformly distributed over phi, cosTheta
  if(1.0 > fMinCosTheta) {
    uz = fMinCosTheta + (1.0 - fMinCosTheta)*G4UniformRand() ;
    ux = std::sqrt((1.0 - uz)*(1.0 + uz)) ;
  } else if (fSigmaTheta > 0.0) {
    ux = G4RandGauss::shoot(0.0,fSigmaTheta);
    uz = std::sqrt((1.0 - ux)*(1.0 + ux));
  }

  G4double phi = twopi*G4UniformRand() ;
  uy = ux;
  ux *= std::cos(phi) ;
  uy *= std::sin(phi) ;
  fDirection = G4ThreeVector(ux,uy,uz) ;

  fParticleGun->SetParticleMomentumDirection(fDirection);

  // Simulation of beam kinetic energy
  G4double kinEnergy = fEnergy;

  if(fGauss == "flatE") {
    kinEnergy  = fEnergy - fSigmaE + 2.*fSigmaE*G4UniformRand();
  } else if(0.0 < fSigmaE) {
    kinEnergy  = fEnergy + G4RandGauss::shoot(0.0,fSigmaE);
  }
  fParticleGun->SetParticleEnergy(kinEnergy);

  if(fVerbose > 0) {
    G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
    G4String particleName = particle->GetParticleName() ;
    G4cout << "Event#  " << fCounter
           << "  Beam particle is generated by PrimaryGeneratorAction "
           << G4endl;
    G4cout << "ParticleName= " << particleName
           << "  PDGcode= " << particle->GetPDGEncoding()
           << std::setprecision(5)
           << "   KinEnergy(GeV)= "
           << fEnergy/GeV
           << "   x(mm)= "
           << x/mm
           << " y(mm)= "
           << y/mm
           << " z(mm)= "
           << z/mm
           << "   ux= "
           << ux
           << " uy= "
           << uy
           << " uz= "
           << uz
           << G4endl;
    }

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorAction::SetBeamEnergy(G4double val)
{
  fEnergy = val;
  if(fEnergy<fDetector->GetMaxEnergy()) fDetector->SetMaxEnergy(fEnergy);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



