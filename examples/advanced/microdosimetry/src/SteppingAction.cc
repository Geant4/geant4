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
// $Id: SteppingAction.cc,v 1.3 2009-04-29 14:17:03 sincerti Exp $
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(RunAction* run,DetectorConstruction* det,PrimaryGeneratorAction* pri)
:Run(run),Detector(det),Primary(pri)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* s)
{ 
 G4int flagParticle=0;
 G4int flagProcess=0;
 
 if (s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "e-")       flagParticle = 1;    
 if (s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "proton")   flagParticle = 2;
 if (s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "hydrogen") flagParticle = 3;
 if (s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha")    flagParticle = 4;
 if (s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha+")   flagParticle = 5;
 if (s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "helium")   flagParticle = 6;

 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAElastic")		flagProcess =11;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAExcitation")		flagProcess =12;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="e-_G4DNAIonisation")		flagProcess =13;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="proton_G4DNAExcitation")	flagProcess =14;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="proton_G4DNAIonisation")	flagProcess =15;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="proton_G4DNAChargeDecrease")	flagProcess =16;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hydrogen_G4DNAIonisation")	flagProcess =17;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hydrogen_G4DNAChargeIncrease")	flagProcess =18;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha_G4DNAExcitation")		flagProcess =19;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha_G4DNAIonisation")		flagProcess =20;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha_G4DNAChargeDecrease")	flagProcess =21;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAExcitation")	flagProcess =22;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAIonisation")	flagProcess =23;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAChargeDecrease")	flagProcess =24;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="alpha+_G4DNAChargeIncrease")	flagProcess =25;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="helium_G4DNAExcitation")	flagProcess =26;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="helium_G4DNAIonisation")	flagProcess =27;
 if (s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="helium_G4DNAChargeIncrease")	flagProcess =28;

 if (
      s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!="initStep"
      &&
      s->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!="Transportation"
    )
 {
   FILE* myFile;
   myFile=fopen("track.txt","a");
   fprintf 
   ( myFile,"%i %i %e %e %e \n",
          flagParticle,
	  flagProcess,
	  (s->GetTrack()->GetPosition().x())/nanometer,
	  (s->GetTrack()->GetPosition().y())/nanometer,
	  (s->GetTrack()->GetPosition().z())/nanometer
    );
   fclose (myFile);
 }
}    
