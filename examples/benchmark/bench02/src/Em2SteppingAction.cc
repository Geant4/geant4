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
//
// $Id: Em2SteppingAction.cc,v 1.1 2002-05-21 11:53:20 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em2SteppingAction.hh"
#include "Em2DetectorConstruction.hh"
#include "Em2RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em2SteppingAction::Em2SteppingAction(Em2DetectorConstruction* det,
                                     Em2RunAction* run)
:Em2Det(det),Em2Run(run)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em2SteppingAction::~Em2SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em2SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
 G4Track* aTrack = aStep->GetTrack();
 const G4VTouchable*  preStepTouchable= aStep->GetPreStepPoint()->GetTouchable();
 const G4VTouchable* postStepTouchable= aStep->GetPostStepPoint()->GetTouchable();

 // energy deposit
 //
 G4int SlideNb(0), RingNb(0);
 if (preStepTouchable->GetHistoryDepth()>0)
 {
   if (Em2Det->GetnRtot()>1)
     RingNb  = preStepTouchable->GetReplicaNumber(1);
   if (Em2Det->GetnLtot()>1)
     SlideNb = preStepTouchable->GetReplicaNumber();
 }
         
 G4double dEStep = aStep->GetTotalEnergyDeposit();
 if (dEStep > 0.)
   Em2Run->fillPerStep(dEStep,SlideNb,RingNb);

 // particle flux
 //  
 if ((Em2Det->GetnLtot()>1)&&
     (postStepTouchable->GetVolume()))
   {
     G4int next = postStepTouchable->GetReplicaNumber();
     if (next != SlideNb)
        Em2Run->particleFlux(aTrack->GetDefinition(), (SlideNb+next)/2);
   }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


