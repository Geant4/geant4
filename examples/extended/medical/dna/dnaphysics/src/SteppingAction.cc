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
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "Analysis.hh"

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4Alpha.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(): G4UserSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
  // Protection

  if (!step->GetPostStepPoint()) return;
  if (!step->GetPostStepPoint()->GetProcessDefinedStep()) return;

  //

  G4double flagParticle=-1.;
  G4double flagProcess=-1.;
  G4double x,y,z,xp,yp,zp;

  // Particle identification

  // The following method avoids the usage of string comparison 

  G4ParticleDefinition * partDef = 
   step->GetTrack()->GetDynamicParticle()->GetDefinition();

  if (partDef == G4Gamma::GammaDefinition())
     flagParticle = 0; 

  if (partDef == G4Electron::ElectronDefinition())
     flagParticle = 1; 

  if (partDef == G4Proton::ProtonDefinition())
     flagParticle = 2; 

  if (partDef == G4Alpha::AlphaDefinition())
     flagParticle = 4; 

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  
  // Usage example
  /*
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
  G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
  G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
  G4ParticleDefinition* heliumDef = instance->GetIon("helium");
  */

  if (partDef == instance->GetIon("hydrogen"))
     flagParticle = 3; 

  if (partDef == instance->GetIon("alpha+"))
     flagParticle = 5; 

  if (partDef == instance->GetIon("helium"))
     flagParticle = 6; 

  // Alternative method (based on string comparison)
  
  /*
  const G4String& particleName = step->GetTrack()->GetDynamicParticle()->
      GetDefinition()->GetParticleName();
      
  if (particleName == "gamma")         flagParticle = 0;
  else if (particleName == "e-")       flagParticle = 1;
  else if (particleName == "proton")   flagParticle = 2;
  else if (particleName == "hydrogen") flagParticle = 3;
  else if (particleName == "alpha")    flagParticle = 4;
  else if (particleName == "alpha+")   flagParticle = 5;
  else if (particleName == "helium")   flagParticle = 6;
  */

  // Process identification

  // Process sub-types are listed in G4PhysicsListHelper.cc
  // or in Geant4-DNA process class implementation files (*.cc)

  G4StepPoint * preStep = step->GetPreStepPoint();
  G4StepPoint * postStep = step->GetPostStepPoint();
  G4int procID = postStep->GetProcessDefinedStep()->GetProcessSubType();

  const G4String& processName = postStep->
      GetProcessDefinedStep()->GetProcessName();

  if (processName=="Capture") flagProcess =1;
  // (no subType and procID exists at the moment for this process)

  else if (flagParticle == 1)

  {
    if      (procID==58) flagProcess =10;
    else if (procID==51) flagProcess =11;
    else if (procID==52) flagProcess =12;
    else if (procID==53) flagProcess =13;
    else if (procID==55) flagProcess =14;
    else if (procID==54) flagProcess =15;
    else if (procID==10) flagProcess =110;
    else if (procID==1)  flagProcess =120;
    else if (procID==2)  flagProcess =130;
  }
        
  else if (flagParticle == 2)

  {
    if      (procID==51) flagProcess =21;
    else if (procID==52) flagProcess =22;
    else if (procID==53) flagProcess =23;
    else if (procID==56) flagProcess =24;
    else if (procID==10) flagProcess =210;
    else if (procID==1)  flagProcess =220;
    else if (procID==2)  flagProcess =230;
    else if (procID==8)  flagProcess =240;
  }
        
  else if (flagParticle == 3)

  {
    if      (procID==51) flagProcess =31;
    else if (procID==52) flagProcess =32;
    else if (procID==53) flagProcess =33;
    else if (procID==57) flagProcess =35;
  }
        
  else if (flagParticle == 4)

  {
    if      (procID==51) flagProcess =41;
    else if (procID==52) flagProcess =42;
    else if (procID==53) flagProcess =43;
    else if (procID==56) flagProcess =44;
    else if (procID==10) flagProcess =410;
    else if (procID==1)  flagProcess =420;
    else if (procID==2)  flagProcess =430;
    else if (procID==8)  flagProcess =440;
  }
        
  else if (flagParticle == 5)

  {
    if      (procID==51) flagProcess =51;
    else if (procID==52) flagProcess =52;
    else if (procID==53) flagProcess =53;
    else if (procID==56) flagProcess =54;
    else if (procID==57) flagProcess =55;
    else if (procID==10) flagProcess =510;
    else if (procID==1)  flagProcess =520;
    else if (procID==2)  flagProcess =530;
    else if (procID==8)  flagProcess =540;
  }
        
  else if (flagParticle == 6)

  {
    if      (procID==51) flagProcess =61;
    else if (procID==52) flagProcess =62;
    else if (procID==53) flagProcess =63;
    else if (procID==57) flagProcess =65;

  }

  else if (processName=="GenericIon_G4DNAIonisation") flagProcess =73;
  else if (processName=="msc")                        flagProcess =710;
  else if (processName=="CoulombScat")                flagProcess =720;
  else if (processName=="ionIoni")                    flagProcess =730;
  else if (processName=="nuclearStopping")            flagProcess =740;
  // (for all GenericIons)
      
  // Alternatively, using process names

  /*
  else if (processName=="e-_G4DNAElectronSolvation")    flagProcess =10;
  else if (processName=="e-_G4DNAElastic")              flagProcess =11;
  else if (processName=="e-_G4DNAExcitation")           flagProcess =12;
  else if (processName=="e-_G4DNAIonisation")           flagProcess =13;
  else if (processName=="e-_G4DNAAttachment")           flagProcess =14;
  else if (processName=="e-_G4DNAVibExcitation")        flagProcess =15;

  else if (processName=="proton_G4DNAElastic")          flagProcess =21;
  else if (processName=="proton_G4DNAExcitation")       flagProcess =22;
  else if (processName=="proton_G4DNAIonisation")       flagProcess =23;
  else if (processName=="proton_G4DNAChargeDecrease")   flagProcess =24;

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

  */

  if (processName!="Transportation")
  {
    x=preStep->GetPosition().x()/nanometer;
    y=preStep->GetPosition().y()/nanometer;
    z=preStep->GetPosition().z()/nanometer;

    xp=postStep->GetPosition().x()/nanometer;
    yp=postStep->GetPosition().y()/nanometer;
    zp=postStep->GetPosition().z()/nanometer;

    // get analysis manager

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    // fill ntuple
    analysisManager->FillNtupleDColumn(0, flagParticle);
    analysisManager->FillNtupleDColumn(1, flagProcess);
    analysisManager->FillNtupleDColumn(2, x);
    analysisManager->FillNtupleDColumn(3, y);
    analysisManager->FillNtupleDColumn(4, z);
    analysisManager->FillNtupleDColumn(5, 
                                     step->GetTotalEnergyDeposit()/eV);

    analysisManager->FillNtupleDColumn(6,
                                          std::sqrt((x-xp)*(x-xp)+
                                           (y-yp)*(y-yp)+(z-zp)*(z-zp)));

    analysisManager->FillNtupleDColumn(7,
                                          (preStep->
                                           GetKineticEnergy() -
                                           postStep->
                                           GetKineticEnergy())/eV );

    analysisManager->FillNtupleDColumn(8, preStep->
                                           GetKineticEnergy()/eV);

    analysisManager->FillNtupleDColumn(9, 
         preStep->GetMomentumDirection()
        *postStep->GetMomentumDirection() );

    analysisManager->FillNtupleIColumn(10, 
                                   G4EventManager::GetEventManager()->
                                   GetConstCurrentEvent()->GetEventID());

    analysisManager->FillNtupleIColumn(11, step->GetTrack()->GetTrackID());

    analysisManager->FillNtupleIColumn(12, step->GetTrack()->GetParentID());

    analysisManager->FillNtupleIColumn(13, 
                                   step->GetTrack()->GetCurrentStepNumber());

    analysisManager->AddNtupleRow();
  }
}    
