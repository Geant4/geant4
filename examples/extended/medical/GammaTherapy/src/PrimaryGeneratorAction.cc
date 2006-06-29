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
#include "Histo.hh"

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
  theMessenger = new PrimaryGeneratorMessenger(this);
  particleGun = new G4ParticleGun();
  counter = 0;
  verbose = 0;
  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;
  sigmaX = 1.5*mm;
  sigmaY = 1.5*mm;
  sigmaZ = 0.0;
  sigmaE = 0.0;
  rMax2  = 2.5*2.5*mm*mm;
  sigmaTheta = 0.0;
//  sigmaTheta = 0.17*degree;
  minCosTheta = 2.0;
  SetBeamEnergy(50.0*MeV);
  position  = G4ThreeVector(x0,y0,z0);
  direction = G4ThreeVector(0.0,0.0,1.0);
  m_gauss = true;
  if(energy < (Histo::GetPointer())->GetMaxEnergy())
              (Histo::GetPointer())->SetMaxEnergy(energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete theMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  counter++ ;
  verbose = (Histo::GetPointer())->GetVerbose();

  // Simulation of beam position
  G4double x = x0;
  G4double y = y0;
  G4double z = fDetector->GetGeneratorPosZ();
  do {
    if(0.0 < sigmaX) x = G4RandGauss::shoot(x0,sigmaX);
    if(0.0 < sigmaY) y = G4RandGauss::shoot(y0,sigmaY);
  } while (x*x + y*y > rMax2);

  position  = G4ThreeVector(x,y,z);
  particleGun->SetParticlePosition(position);

  // Simulation of beam direction
  G4double ux = direction.x();
  G4double uy = direction.y();
  G4double uz = direction.z();

  // Beam particles are uniformly distributed over phi, cosTheta
  if(1.0 > minCosTheta) {
    uz = minCosTheta + (1.0 - minCosTheta)*G4UniformRand() ;
    ux = std::sqrt(1.0 - uz*uz) ;
  } else if (sigmaTheta > 0.0) {
    ux = G4RandGauss::shoot(0.0,sigmaTheta);
    uz = std::sqrt(1.0 - ux*ux);
  }

  G4double phi = twopi*G4UniformRand() ;
  uy = ux ;
  ux *= std::cos(phi) ;
  uy *= std::sin(phi) ;
  direction = G4ThreeVector(ux,uy,uz) ;

  direction = direction.unit();
  particleGun->SetParticleMomentumDirection(direction);
  G4ParticleDefinition* particle = particleGun->GetParticleDefinition();

  // Simulation of beam kinetic energy
  G4double kinEnergy = energy;

  if(m_gauss == "flatE") kinEnergy  = energy - sigmaE + 2.*sigmaE*G4UniformRand();
  else if(0.0 < sigmaE)  kinEnergy  = energy + G4RandGauss::shoot(0.0,sigmaE);

  particleGun->SetParticleEnergy(kinEnergy);

  G4String particleName = particle->GetParticleName() ;

  if(verbose > 0) {
    G4cout << "Event#  " << counter
           << "  Beam particle is generated by PrimaryGeneratorAction "
           << G4endl;
    G4cout << "ParticleName= " << particleName
           << "  PDGcode= " << particle->GetPDGEncoding()
           << std::setprecision(5)
	   << "   KinEnergy(GeV)= "
	   << energy/GeV
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

  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorAction::SetBeamSigmaE(G4double val)
{
  sigmaE = val;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorAction::SetBeamEnergy(G4double val)
{
  energy = val;
  if(energy < (Histo::GetPointer())->GetMaxEnergy())
              (Histo::GetPointer())->SetMaxEnergy(energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



