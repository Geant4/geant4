//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#define hTestPrimaryGeneratorAction_CPP 

//---------------------------------------------------------------------------
//
// ClassName:   hTestPrimaryGeneratorAction
//  
// Description: Generate primary beam 
//
// Authors:    0.6.04.01 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "hTestPrimaryGeneratorAction.hh"
#include "hTestPrimaryGeneratorMessenger.hh"
#include "Randomize.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "hTestHisto.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorAction::hTestPrimaryGeneratorAction(
			     hTestDetectorConstruction* det):
  theDet(det)
{
  InitializeMe();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::InitializeMe()
{
  verbose = theDet->GetVerbose();
  theMessenger = new hTestPrimaryGeneratorMessenger(this);
  particleGun = new G4ParticleGun();
  counter = 0;
  x0 = 0.0; 
  y0 = 0.0;
  z0 = 0.0;
  sigmaX = 0.0;
  sigmaY = 0.0;
  sigmaZ = 0.0;
  sigmaE = 0.0;
  minCosTheta = 1.0;
  SetBeamEnergy(10.0*MeV);
  position  = G4ThreeVector(x0,y0,z0);
  direction = G4ThreeVector(0.0,0.0,1.0);
  m_gauss = true;
  if(energy < (hTestHisto::GetPointer())->GetMaxEnergy())
              (hTestHisto::GetPointer())->SetMaxEnergy(energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorAction::~hTestPrimaryGeneratorAction()
{
  delete particleGun;
  delete theMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  counter++ ;
  verbose = theDet->GetVerbose();

  // Simulation of beam position
  G4double x = x0;
  G4double y = y0;
  G4double z = z0;
  if(0.0 < sigmaX) x += G4RandGauss::shoot(0.0,sigmaX);
  if(0.0 < sigmaY) y += G4RandGauss::shoot(0.0,sigmaY);
  if(0.0 < sigmaZ) z += G4RandGauss::shoot(0.0,sigmaZ);
  position  = G4ThreeVector(x,y,z);
  particleGun->SetParticlePosition(position);

  // Simulation of beam direction
  G4double ux = direction.x();
  G4double uy = direction.y();
  G4double uz = direction.z();

  // Beam particles are uniformly distributed over phi, cosTheta 
  if(1.0 > minCosTheta) {
    uz = minCosTheta + (1.0 - minCosTheta)*G4UniformRand() ;
    ux = sqrt(1.0 - uz*uz) ;
    uy = ux ;
    G4double phi = 360.0*deg*G4UniformRand() ;
    ux *= cos(phi) ;
    uy *= sin(phi) ;
    direction = G4ThreeVector(ux,uy,uz) ;
  }

  direction = direction.unit();
  particleGun->SetParticleMomentumDirection(direction);
  G4ParticleDefinition* particle = particleGun->GetParticleDefinition();
  G4double mass = particle->GetPDGMass();

  // Simulation of beam kinetic energy
  G4double kinEnergy = energy;

  if(m_gauss == "flatE") kinEnergy  = minE + (maxE-minE)*G4UniformRand();

  else if(m_gauss == "flatBeta") {
         G4double beta = minBeta + (maxBeta-minBeta)*G4UniformRand();
         kinEnergy = mass*(1./sqrt(1. - beta*beta) - 1.);
  }
  else if(0.0 < sigmaE) kinEnergy += G4RandGauss::shoot(0.0,sigmaE);
   

  if(0.0 > kinEnergy) kinEnergy = 0.0;
  particleGun->SetParticleEnergy(kinEnergy);

  G4String particleName = particle->GetParticleName() ;

  if(verbose > 0) {
    G4cout << "Event#  " << counter 
           << "  Beam particle is generated by hTestPrimaryGeneratorAction " 
           << G4endl;
    G4cout << "ParticleName= " << particleName 
           << "  PDGcode= " << particle->GetPDGEncoding()
           << G4std::setprecision(5) 
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
	   << endl;
    }

  particleGun->GeneratePrimaryVertex(anEvent);
  if(verbose > 1) G4cout << "hTestPrimaryGeneratorAction: BeamOn" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::SetBeamBeta(G4double val)
{
  G4ParticleDefinition* particle = particleGun->GetParticleDefinition();
  G4double mass = particle->GetPDGMass();
  if(val > 0. && val < 1.) energy = mass*(1./sqrt(1.-val*val) - 1.);
  G4cout << "hTestPrimaryGeneratorAction: KinEnergy(MeV)= " 
         << energy/MeV << G4endl;
  minE = energy;
  maxE = energy;
  minBeta = val;
  maxBeta = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::SetSigmaBeta(G4double val)
{
  G4ParticleDefinition* particle = particleGun->GetParticleDefinition();
  G4double mass = particle->GetPDGMass();
  if(val > 0. && val < 1.) {
    sigmaE = mass*(1./sqrt(1.-val*val) - 1.);
    G4double gamma = energy/mass + 1.;
    G4double beta0 = sqrt(1. - 1./(gamma*gamma));
    G4double beta  = beta0 + val;
    if (beta >= 1.) beta = 0.9999;
    maxBeta = beta;
    maxE = mass*(1./sqrt(1.-beta*beta) - 1.);
    beta  = beta0 - val;
    if (beta <= 0.) beta = 0.0001;
    minBeta = beta;
    minE = mass*(1./sqrt(1.-beta*beta) - 1.);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::SetBeamSigmaE(G4double val) 
{
  G4ParticleDefinition* particle = particleGun->GetParticleDefinition();
  G4double mass = particle->GetPDGMass();
  sigmaE = val; 
  minE = energy - sigmaE;
  G4double gamma = minE/mass + 1.;
  minBeta = sqrt(1. - 1./(gamma*gamma));
  maxE = energy + sigmaE;
  gamma = maxE/mass + 1.;
  maxBeta = sqrt(1. - 1./(gamma*gamma));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::SetBeamEnergy(G4double val) 
{
  G4ParticleDefinition* particle = particleGun->GetParticleDefinition();
  G4double mass = particle->GetPDGMass();
  energy = val;
  minE = energy - sigmaE;
  G4double gamma = minE/mass + 1.;
  minBeta = sqrt(1. - 1./(gamma*gamma));
  maxE = energy + sigmaE;
  gamma = maxE/mass + 1.;
  maxBeta = sqrt(1. - 1./(gamma*gamma));
  if(energy < (hTestHisto::GetPointer())->GetMaxEnergy())
              (hTestHisto::GetPointer())->SetMaxEnergy(energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....








