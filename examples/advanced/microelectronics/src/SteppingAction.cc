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
// -------------------------------------------------------------------
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
 G4double flagParticle=0.;
 G4double flagProcess=0.;
 G4double x,y,z,xp,yp,zp;
 
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "e-")      	  flagParticle = 1;    
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "proton")  	  flagParticle = 2;
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "GenericIon")    flagParticle = 3;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="msc")				flagProcess =10;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4MicroElecElastic")		flagProcess =11;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4MicroElecInelastic")	flagProcess =12;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="eCapture")			flagProcess =13;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="p_G4MicroElecInelastic")	flagProcess =14;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="ion_G4MicroElecInelastic")	flagProcess =15;

 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hIoni")				flagProcess =16;
 if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="eIoni")				flagProcess =17;

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
