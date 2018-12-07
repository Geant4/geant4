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
/// \file B5PrimaryGeneratorAction.cc
/// \brief Implementation of the B5PrimaryGeneratorAction class

#include "B5PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5PrimaryGeneratorAction::B5PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),     
  fParticleGun(nullptr), fMessenger(nullptr), 
  fPositron(nullptr), fMuon(nullptr), fPion(nullptr), 
  fKaon(nullptr), fProton(nullptr),
  fMomentum(1000.*MeV),
  fSigmaMomentum(50.*MeV),
  fSigmaAngle(2.*deg),
  fRandomizePrimary(true)
{
  G4int nofParticles = 1;
  fParticleGun  = new G4ParticleGun(nofParticles);
  
  auto particleTable = G4ParticleTable::GetParticleTable();
  fPositron = particleTable->FindParticle("e+");
  fMuon = particleTable->FindParticle("mu+");
  fPion = particleTable->FindParticle("pi+");
  fKaon = particleTable->FindParticle("kaon+");
  fProton = particleTable->FindParticle("proton");
  
  // default particle kinematics
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-8.*m));
  fParticleGun->SetParticleDefinition(fPositron);
  
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5PrimaryGeneratorAction::~B5PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  G4ParticleDefinition* particle;  
  if (fRandomizePrimary) {
    G4int i = (int)(5.*G4UniformRand());
    switch(i) {
      case 0:
          particle = fPositron;
          break;
      case 1:
          particle = fMuon;
          break;
      case 2:
          particle = fPion;
          break;
      case 3:
          particle = fKaon;
          break;
      default:
          particle = fProton;
          break;
    }
    fParticleGun->SetParticleDefinition(particle);
  }
  else {
    particle = fParticleGun->GetParticleDefinition();
  }
  
  auto pp = fMomentum + (G4UniformRand()-0.5)*fSigmaMomentum;
  auto mass = particle->GetPDGMass();
  auto ekin = std::sqrt(pp*pp+mass*mass)-mass;
  fParticleGun->SetParticleEnergy(ekin);
  
  auto angle = (G4UniformRand()-0.5)*fSigmaAngle;
  fParticleGun->SetParticleMomentumDirection(
                  G4ThreeVector(std::sin(angle),0.,std::cos(angle)));
  
  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5PrimaryGeneratorAction::DefineCommands()
{
  // Define /B5/generator command directory using generic messenger class
  fMessenger 
    = new G4GenericMessenger(this, 
                             "/B5/generator/", 
                             "Primary generator control");
            
  // momentum command
  auto& momentumCmd
    = fMessenger->DeclarePropertyWithUnit("momentum", "GeV", fMomentum, 
        "Mean momentum of primaries.");
  momentumCmd.SetParameterName("p", true);
  momentumCmd.SetRange("p>=0.");                                
  momentumCmd.SetDefaultValue("1.");
  // ok
  //momentumCmd.SetParameterName("p", true);
  //momentumCmd.SetRange("p>=0.");                                
  
  // sigmaMomentum command
  auto& sigmaMomentumCmd
    = fMessenger->DeclarePropertyWithUnit("sigmaMomentum",
        "MeV", fSigmaMomentum, "Sigma momentum of primaries.");
  sigmaMomentumCmd.SetParameterName("sp", true);
  sigmaMomentumCmd.SetRange("sp>=0.");                                
  sigmaMomentumCmd.SetDefaultValue("50.");

  // sigmaAngle command
  auto& sigmaAngleCmd
    = fMessenger->DeclarePropertyWithUnit("sigmaAngle", "deg", fSigmaAngle, 
        "Sigma angle divergence of primaries.");
  sigmaAngleCmd.SetParameterName("t", true);
  sigmaAngleCmd.SetRange("t>=0.");                                
  sigmaAngleCmd.SetDefaultValue("2.");

  // randomizePrimary command
  auto& randomCmd
    = fMessenger->DeclareProperty("randomizePrimary", fRandomizePrimary);
  G4String guidance
    = "Boolean flag for randomizing primary particle types.\n";   
  guidance
    += "In case this flag is false, you can select the primary particle\n";
  guidance += "  with /gun/particle command.";                               
  randomCmd.SetGuidance(guidance);
  randomCmd.SetParameterName("flg", true);
  randomCmd.SetDefaultValue("true");
}

//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
