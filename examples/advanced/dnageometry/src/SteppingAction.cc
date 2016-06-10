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
// and the DNA geometry given in the Geom_DNA example 
// shall cite the following Geant4-DNA collaboration publications:
// [1] NIM B 298 (2013) 47-54
// [2] Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(RunAction* run,DetectorConstruction* det,PrimaryGeneratorAction* pri, HistoManager* his)
:Run(run),Detector(det),Primary(pri),Histo(his)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
 G4double flagParticle=0.;
 G4double flagProcess=0.;
 G4double flagVolume=0.;
 G4double x,y,z,xp,yp,zp;
 G4double dE;
 
 dE=step->GetTotalEnergyDeposit()/eV;
 
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "e-")       flagParticle = 10;    
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "proton")   flagParticle = 20;
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "hydrogen") flagParticle = 30;
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha")    flagParticle = 40;
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha+")   flagParticle = 50;
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "helium")   flagParticle = 60;


 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAElastic")		flagProcess =11;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAExcitation")		flagProcess =12;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAIonisation")		flagProcess =13;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAAttachment")		flagProcess =14;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAVibExcitation")	flagProcess =15;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="eCapture")			flagProcess =16;
//  if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="msc")				flagProcess =17;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="proton_G4DNAExcitation")	flagProcess =21;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="proton_G4DNAIonisation")	flagProcess =22;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="proton_G4DNAChargeDecrease")	flagProcess =23;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hydrogen_G4DNAExcitation")	flagProcess =31;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hydrogen_G4DNAIonisation")	flagProcess =32;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hydrogen_G4DNAChargeIncrease")	flagProcess =33;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha_G4DNAExcitation")		flagProcess =41;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha_G4DNAIonisation")		flagProcess =42;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha_G4DNAChargeDecrease")	flagProcess =43;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAExcitation")	flagProcess =51;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAIonisation")	flagProcess =52;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAChargeDecrease")	flagProcess =53;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAChargeIncrease")	flagProcess =54;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="helium_G4DNAExcitation")	flagProcess =61;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="helium_G4DNAIonisation")	flagProcess =62;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="helium_G4DNAChargeIncrease")	flagProcess =63;

// if (step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hIoni")				flagProcess =24;
// if (step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName()=="eIoni")				flagProcess =18;
 
 if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="physi sugar 2") flagVolume=1;
  
 if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="physi sugar 4") flagVolume=2;
 
 if (flagVolume !=0 && dE!=0 )
 {
 
   x=step->GetPreStepPoint()->GetPosition().x()/nanometer;
   y=step->GetPreStepPoint()->GetPosition().y()/nanometer;
   z=step->GetPreStepPoint()->GetPosition().z()/nanometer;
   xp=step->GetPostStepPoint()->GetPosition().x()/nanometer;
   yp=step->GetPostStepPoint()->GetPosition().y()/nanometer;
   zp=step->GetPostStepPoint()->GetPosition().z()/nanometer;
   
   Histo->FillNtupleDColumn(0, flagParticle);
   Histo->FillNtupleDColumn(1, flagProcess);
   Histo->FillNtupleDColumn(2, flagVolume);
   Histo->FillNtupleDColumn(3, xp);
   Histo->FillNtupleDColumn(4, yp);
   Histo->FillNtupleDColumn(5, zp);
   Histo->FillNtupleDColumn(6, dE );
   Histo->FillNtupleDColumn(7, std::sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp)+(z-zp)*(z-zp)));
   
   Histo->AddNtupleRow();   
 }
 
 
}    
