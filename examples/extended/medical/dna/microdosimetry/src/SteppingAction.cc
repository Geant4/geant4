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
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $ID$
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "Analysis.hh"

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "CommandLineParser.hh"

using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction() : G4UserSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
  G4double flagParticle=-1;
  G4double flagProcess=-1;
  G4double x,y,z,xp,yp,zp;

  // Process sub-types are listed in G4PhysicsListHelper.cc

  /*
  // The following alternative method avoids the usage of string comparison 

  if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
      == G4Electron::ElectronDefinition())
     flagParticle = 1; 

  if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
  == G4Proton::ProtonDefinition())
     flagParticle = 2; 

  if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
  == G4Alpha::AlphaDefinition())
     flagParticle = 4; 

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
  G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
  G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
  G4ParticleDefinition* heliumDef = instance->GetIon("helium");

  if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
  == instance->GetIon("hydrogen"))
     flagParticle = 3; 

  if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
  == instance->GetIon("alpha+"))
     flagParticle = 5; 

  if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
  == instance->GetIon("helium"))
     flagParticle = 6; 
 */
 
  const G4String& particleName = step->GetTrack()->GetDynamicParticle()->
      GetDefinition()->GetParticleName();
      
  const G4String& processName = step->GetPostStepPoint()->
      GetProcessDefinedStep()->GetProcessName();

  if (particleName == "gamma")         flagParticle = 0;
  else if (particleName == "e-")       flagParticle = 1;
  else if (particleName == "proton")   flagParticle = 2;
  else if (particleName == "hydrogen") flagParticle = 3;
  else if (particleName == "alpha")    flagParticle = 4;
  else if (particleName == "alpha+")   flagParticle = 5;
  else if (particleName == "helium")   flagParticle = 6;
  
  if (processName=="e-_G4DNAElectronSolvation")         flagProcess =10;
  else if (processName=="e-_G4DNAElastic")              flagProcess =11;
  else if (processName=="e-_G4DNAExcitation")           flagProcess =12;
  else if (processName=="e-_G4DNAIonisation")           flagProcess =13;
  else if (processName=="e-_G4DNAAttachment")           flagProcess =14;
  else if (processName=="e-_G4DNAVibExcitation")        flagProcess =15;
  else if (processName=="eCapture")                     flagProcess =16;
  else if (processName=="msc")                          flagProcess =17;
  else if (processName=="eIoni")                        flagProcess =130;

  else if (processName=="proton_G4DNAElastic")          flagProcess =21;
  else if (processName=="proton_G4DNAExcitation")       flagProcess =22;
  else if (processName=="proton_G4DNAIonisation")       flagProcess =23;
  else if (processName=="proton_G4DNAChargeDecrease")   flagProcess =24;
  else if (processName=="hIoni")                        flagProcess =230;

  else if (processName=="hydrogen_G4DNAElastic")        flagProcess =31;
  else if (processName=="hydrogen_G4DNAExcitation")     flagProcess =32;
  else if (processName=="hydrogen_G4DNAIonisation")     flagProcess =33;
  else if (processName=="hydrogen_G4DNAChargeIncrease") flagProcess =35;

  else if (processName=="alpha_G4DNAElastic")           flagProcess =41;
  else if (processName=="alpha_G4DNAExcitation")        flagProcess =42;
  else if (processName=="alpha_G4DNAIonisation")        flagProcess =43;
  else if (processName=="alpha_G4DNAChargeDecrease")    flagProcess =44;

  else if (processName=="alpha+_G4DNAElastic")          flagProcess =51;
  else if (processName=="alpha+_G4DNAExcitation")       flagProcess =52;
  else if (processName=="alpha+_G4DNAIonisation")       flagProcess =53;
  else if (processName=="alpha+_G4DNAChargeDecrease")   flagProcess =54;
  else if (processName=="alpha+_G4DNAChargeIncrease")   flagProcess =55;

  else if (processName=="helium_G4DNAElastic")          flagProcess =61;
  else if (processName=="helium_G4DNAExcitation")       flagProcess =62;
  else if (processName=="helium_G4DNAIonisation")       flagProcess =63;
  else if (processName=="helium_G4DNAChargeIncrease")   flagProcess =65;

  else if (processName=="GenericIon_G4DNAIonisation")   flagProcess =73;

  if (processName!="Transportation")
  {
    x=step->GetPreStepPoint()->GetPosition().x()/nanometer;
    y=step->GetPreStepPoint()->GetPosition().y()/nanometer;
    z=step->GetPreStepPoint()->GetPosition().z()/nanometer;
    xp=step->GetPostStepPoint()->GetPosition().x()/nanometer;
    yp=step->GetPostStepPoint()->GetPosition().y()/nanometer;
    zp=step->GetPostStepPoint()->GetPosition().z()/nanometer;

    CommandLineParser* parser = CommandLineParser::GetParser();
    Command* command(0);
    if((command = parser->GetCommandIfActive("-out"))==0) return;

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillNtupleDColumn(0, flagParticle);
    analysisManager->FillNtupleDColumn(1, flagProcess);
    analysisManager->FillNtupleDColumn(2, x);
    analysisManager->FillNtupleDColumn(3, y);
    analysisManager->FillNtupleDColumn(4, z);
    analysisManager->FillNtupleDColumn(5, step->GetTotalEnergyDeposit()/eV);
    analysisManager->FillNtupleDColumn(6,
                                          std::sqrt((x-xp)*(x-xp)+
                                           (y-yp)*(y-yp)+(z-zp)*(z-zp))/nm);
    analysisManager->FillNtupleDColumn(7, (step->GetPreStepPoint()->
                                           GetKineticEnergy() - 
                                           step->GetPostStepPoint()->
                                           GetKineticEnergy())/eV );
    analysisManager->AddNtupleRow();
  }
}    
