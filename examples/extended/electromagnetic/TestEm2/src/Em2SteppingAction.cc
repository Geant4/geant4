// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2SteppingAction.cc,v 1.2 1999-12-15 14:49:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em2SteppingAction.hh"
#include "Em2DetectorConstruction.hh"
#include "Em2RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2SteppingAction::Em2SteppingAction(Em2DetectorConstruction* det,
                                     Em2RunAction* run)
:Em2Det(det),Em2Run(run)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2SteppingAction::~Em2SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
 G4Track* aTrack = aStep->GetTrack();

 // energy deposit
 //
 G4int SlideNb(0), RingNb(0);
 if (Em2Det->GetnRtot()>1)
    RingNb  = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1);
 if (Em2Det->GetnLtot()>1)
    SlideNb = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber();
             
 G4double dEStep = aStep->GetTotalEnergyDeposit();
 if (dEStep > 0.) Em2Run->fillPerStep(dEStep,SlideNb,RingNb);

 // particle flux
 //  
 if ((Em2Det->GetnLtot()>1)&&
     (aStep->GetPostStepPoint()->GetTouchable()->GetVolume()))
   {
     G4int next = aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber();
     if (next != SlideNb)
        Em2Run->particleFlux(aTrack->GetDefinition(), (SlideNb+next)/2);
   }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


