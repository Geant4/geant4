// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17SteppingAction.cc,v 1.1 1999-11-30 18:01:56 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst17DetectorConstruction.hh"
#include "G4EnergyLossTables.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "Tst17SteppingAction.hh"
#include "Tst17PrimaryGeneratorAction.hh"
#include "Tst17EventAction.hh"
#include "Tst17RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "Tst17SteppingMessenger.hh"
#include "G4ios.hh"
#include <iomanip.h>
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst17SteppingAction::Tst17SteppingAction(Tst17DetectorConstruction* DET,
                                         Tst17EventAction* EA,
                                         Tst17RunAction* RA)
:detector (DET),eventaction (EA),runaction (RA),steppingMessenger(NULL),
 IDold(-1) ,evnoold(-1)
{
  steppingMessenger = new Tst17SteppingMessenger(this);

 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst17SteppingAction::~Tst17SteppingAction()
{
  delete steppingMessenger ;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 

  G4double Edep,Theta,Thetaback,Ttrans,Tback,Tsec,Egamma,xend,yend,zend,rend ;
  G4double Tkin ;
  G4int evno = eventaction->GetEventno() ; 

  IDnow = evno+10000*(aStep->GetTrack()->GetTrackID())+
          100000000*(aStep->GetTrack()->GetParentID()); 
  if(IDnow != IDold)
  {
   IDold=IDnow ;

   if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber"){

    if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()) == "e-") &&
       ((aStep->GetTrack()->GetTrackID() != 1) ||
       (aStep->GetTrack()->GetParentID() != 0)) ){

      eventaction->AddCharged() ;
      eventaction->AddE() ;
      Tsec = aStep->GetTrack()->GetKineticEnergy() ;  // !!!!!!!!!!!!
      Tsec += aStep->GetTotalEnergyDeposit() ;        // !!!!!!!!!!!!
      runaction->FillTsec(Tsec) ;
    }

    else if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()) == "e+") &&
       ((aStep->GetTrack()->GetTrackID() != 1) ||
       (aStep->GetTrack()->GetParentID() != 0)) ) 
    {
      eventaction->AddCharged() ;
      eventaction->AddP() ;
      Tsec = aStep->GetTrack()->GetKineticEnergy() ;  // !!!!!!!!!!!!
      Tsec += aStep->GetTotalEnergyDeposit() ;        // !!!!!!!!!!!!
      runaction->FillTsec(Tsec) ;
    }

    else if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()) == "gamma") &&
       ((aStep->GetTrack()->GetTrackID() != 1) ||
       (aStep->GetTrack()->GetParentID() != 0)) ){

      eventaction->AddNeutral() ;
    }
   }
  }

  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber"){


    if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()) == "e-") ||
       ((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()) == "e+")) {

          eventaction->CountStepsCharged() ;
    }

    if ((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()) == "gamma") {

          eventaction->CountStepsNeutral() ;
    }
  }

  if (((aStep->GetTrack()->GetTrackID() == 1) &&
       (aStep->GetTrack()->GetParentID()== 0)) ||
      (aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() ==
       Tst17PrimaryGeneratorAction::GetPrimaryName()) &&
      (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber") &&
      (aStep->GetTrack()->GetNextVolume()->GetName()=="World") &&
      (aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "Transportation") &&
      (aStep->GetTrack()->GetMomentumDirection().z()>0.)){
    
    Ttrans = aStep->GetTrack()->GetKineticEnergy() ;
    runaction->FillTt(Ttrans) ;
  }
       
  if (((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber") &&
       (aStep->GetTrack()->GetNextVolume()->GetName()=="World") &&
       (aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "Transportation") &&
       (aStep->GetTrack()->GetMomentumDirection().z()>0.) &&
       (aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "gamma"))){

    Egamma = aStep->GetTrack()->GetKineticEnergy() ;
    runaction->FillGammaSpectrum(Egamma) ;
  }
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

