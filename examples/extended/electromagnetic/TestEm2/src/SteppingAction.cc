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
// $Id: SteppingAction.cc,v 1.3 2004/11/24 10:46:30 maire Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
                                     RunAction* run)
:Det(det),Run(run)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
 G4Track* aTrack = aStep->GetTrack();
 const G4VTouchable*  preStepTouchable= aStep->GetPreStepPoint()->GetTouchable();
 const G4VTouchable* postStepTouchable= aStep->GetPostStepPoint()->GetTouchable();

 // energy deposit
 //
 G4int SlideNb(0), RingNb(0);
 if (preStepTouchable->GetHistoryDepth()>0)
 {
   if (Det->GetnRtot()>1)
     RingNb  = preStepTouchable->GetReplicaNumber(1);   
///     RingNb  = preStepTouchable->GetCopyNumber(1);
   if (Det->GetnLtot()>1)
     SlideNb = preStepTouchable->GetReplicaNumber();   
///     SlideNb = preStepTouchable->GetCopyNumber();
 }
         
 G4double dEStep = aStep->GetTotalEnergyDeposit();
 if (dEStep > 0.)
   Run->fillPerStep(dEStep,SlideNb,RingNb);

 // particle flux
 //  
 if ((Det->GetnLtot()>1)&&
     (postStepTouchable->GetVolume()))
   {
     G4int next = postStepTouchable->GetReplicaNumber();   
///     G4int next = postStepTouchable->GetCopyNumber();
     if (next != SlideNb)
        Run->particleFlux(aTrack->GetDefinition(), (SlideNb+next)/2);
   }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


