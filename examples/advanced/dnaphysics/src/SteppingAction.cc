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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
 G4double flagParticle=-1.;
 G4double flagProcess=-1.;
 G4double x,y,z,xp,yp,zp;
 
 // Process sub-types are listed in G4PhysicsListHelper.cc
 
/*
 // The following method avoids the usage of string comparison 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Electron::ElectronDefinition())       
    flagParticle = 1; 
 
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Proton::ProtonDefinition())       
    flagParticle = 2; 
 
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Alpha::AlphaDefinition())       
    flagParticle = 4; 

    G4DNAGenericIonsManager *instance;
    instance = G4DNAGenericIonsManager::Instance();
    G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
    G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
    G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
    G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
    G4ParticleDefinition* heliumDef = instance->GetIon("helium");
 
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() == instance->GetIon("hydrogen"))       
    flagParticle = 3; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() == instance->GetIon("alpha+"))       
    flagParticle = 5; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() == instance->GetIon("helium"))       
    flagParticle = 6; 
*/
 
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "e-")       flagParticle = 1;    
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "proton")   flagParticle = 2;
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "hydrogen") flagParticle = 3;
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha")    flagParticle = 4;
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha+")   flagParticle = 5;
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "helium")   flagParticle = 6;

 
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAElastic")		flagProcess =11;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAExcitation")		flagProcess =12;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAIonisation")		flagProcess =13;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAAttachment")		flagProcess =14;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAVibExcitation")	flagProcess =15;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="proton_G4DNAExcitation")	flagProcess =17;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="proton_G4DNAIonisation")	flagProcess =18;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="proton_G4DNAChargeDecrease")	flagProcess =19;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hydrogen_G4DNAExcitation")		flagProcess =20;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hydrogen_G4DNAIonisation")		flagProcess =21;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hydrogen_G4DNAChargeIncrease")	flagProcess =22;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha_G4DNAExcitation")		flagProcess =23;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha_G4DNAIonisation")		flagProcess =24;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha_G4DNAChargeDecrease")		flagProcess =25;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAExcitation")	flagProcess =26;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAIonisation")	flagProcess =27;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAChargeDecrease")	flagProcess =28;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAChargeIncrease")	flagProcess =29;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="helium_G4DNAExcitation")	flagProcess =30;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="helium_G4DNAIonisation")	flagProcess =31;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="helium_G4DNAChargeIncrease")	flagProcess =32;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!="Transportation")
 {  
   x=step->GetPreStepPoint()->GetPosition().x()/nanometer;
   y=step->GetPreStepPoint()->GetPosition().y()/nanometer;
   z=step->GetPreStepPoint()->GetPosition().z()/nanometer;
   xp=step->GetPostStepPoint()->GetPosition().x()/nanometer;
   yp=step->GetPostStepPoint()->GetPosition().y()/nanometer;
   zp=step->GetPostStepPoint()->GetPosition().z()/nanometer;
   
   // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // fill ntuple
  analysisManager->FillNtupleDColumn(0, flagParticle);
  analysisManager->FillNtupleDColumn(1, flagProcess);
  analysisManager->FillNtupleDColumn(2, x);
  analysisManager->FillNtupleDColumn(3, y);
  analysisManager->FillNtupleDColumn(4, z);
  analysisManager->FillNtupleDColumn(5, step->GetTotalEnergyDeposit()/eV);
  analysisManager->FillNtupleDColumn(6, std::sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp)+(z-zp)*(z-zp))/nm);
  analysisManager->FillNtupleDColumn(7, (step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy())/eV );
  analysisManager->AddNtupleRow();  
 }
}    
