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
// $Id: Tst26SteppingAction.cc,v 1.2 2003-02-01 18:14:59 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26SteppingAction.hh"
#include "Tst26DetectorConstruction.hh"
#include "Tst26RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26SteppingAction::Tst26SteppingAction(Tst26DetectorConstruction* det,
                                     Tst26RunAction* run)
:Tst26Det(det),Tst26Run(run)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26SteppingAction::~Tst26SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
 G4Track* aTrack = aStep->GetTrack();
 const G4VTouchable*  preStepTouchable= aStep->GetPreStepPoint()->GetTouchable();
 const G4VTouchable* postStepTouchable= aStep->GetPostStepPoint()->GetTouchable();

 // energy deposit
 //
 G4int SlideNb(0), RingNb(0);
 if (preStepTouchable->GetHistoryDepth()>0)
 {
   if (Tst26Det->GetnRtot()>1)
     RingNb  = preStepTouchable->GetReplicaNumber(1);
   if (Tst26Det->GetnLtot()>1)
     SlideNb = preStepTouchable->GetReplicaNumber();
 }
         
 G4double dEStep = aStep->GetTotalEnergyDeposit();
 if (dEStep > 0.)
   Tst26Run->fillPerStep(dEStep,SlideNb,RingNb);

 // particle flux
 //  
 if ((Tst26Det->GetnLtot()>1)&&
     (postStepTouchable->GetVolume()))
   {
     G4int next = postStepTouchable->GetReplicaNumber();
     if (next != SlideNb)
        Tst26Run->particleFlux(aTrack->GetDefinition(), (SlideNb+next)/2);
   }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


