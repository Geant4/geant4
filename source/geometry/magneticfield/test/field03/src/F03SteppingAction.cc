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
//
// $Id: F03SteppingAction.cc,v 1.4 2006-06-29 18:30:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F03DetectorConstruction.hh"
#include "G4EnergyLossTables.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "F03SteppingAction.hh"
#include "F03PrimaryGeneratorAction.hh"
#include "F03EventAction.hh"
#include "F03RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "F03SteppingMessenger.hh"
#include "G4ios.hh"
#include <iomanip>
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F03SteppingAction::F03SteppingAction(F03DetectorConstruction* DET,
                                     F03EventAction* EA,
                                     F03RunAction* RA)
:detector (DET),eventaction (EA),runaction (RA),steppingMessenger(NULL),
 IDold(-1) ,evnoold(-1)
{
  steppingMessenger = new F03SteppingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F03SteppingAction::~F03SteppingAction()
{
  delete steppingMessenger ;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F03SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 

  G4double Edep,Theta,Thetaback,Ttrans,Tback,Tsec,Egamma,xend,yend,zend,rend ;
  G4double Tkin ;
  G4int evno = eventaction->GetEventno() ; 

  IDnow = evno+10000*(aStep->GetTrack()->GetTrackID())+
          100000000*(aStep->GetTrack()->GetParentID()); 
  if(IDnow != IDold)
  {
   IDold=IDnow ;
   if(
      (((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e-") &&
       ((aStep->GetTrack()->GetTrackID() != 1) ||
       (aStep->GetTrack()->GetParentID() != 0)) )
       ||
      (((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e+") &&
       ((aStep->GetTrack()->GetTrackID() != 1) ||
       (aStep->GetTrack()->GetParentID() != 0)) )
     )
        runaction->Fillvertexz(aStep->GetTrack()->GetVertexPosition().x());

   if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber")
   {
    if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e-") &&
       ((aStep->GetTrack()->GetTrackID() != 1) ||
       (aStep->GetTrack()->GetParentID() != 0)) ) 
    {
        eventaction->AddCharged() ;
        eventaction->AddE() ;
        Tsec = aStep->GetTrack()->GetKineticEnergy() ;  // !!!!!!!!!!!!
        Tsec += aStep->GetTotalEnergyDeposit() ;        // !!!!!!!!!!!!
        runaction->FillTsec(Tsec) ;
    }
    else
    if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e+") &&
       ((aStep->GetTrack()->GetTrackID() != 1) ||
       (aStep->GetTrack()->GetParentID() != 0)) ) 
    {
        eventaction->AddCharged() ;
        eventaction->AddP() ;
        Tsec = aStep->GetTrack()->GetKineticEnergy() ;  // !!!!!!!!!!!!
        Tsec += aStep->GetTotalEnergyDeposit() ;        // !!!!!!!!!!!!
        runaction->FillTsec(Tsec) ;
    }
    else
    if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "gamma") &&
       ((aStep->GetTrack()->GetTrackID() != 1) ||
       (aStep->GetTrack()->GetParentID() != 0)) ) 
    {
        eventaction->AddNeutral() ;
    }
   }
  }

  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber")
  {
    if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e-") 
              ||
       ((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e+"))  
          eventaction->CountStepsCharged() ;

    if ((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "gamma") 
          eventaction->CountStepsNeutral() ;
  }

  if (
      (((aStep->GetTrack()->GetTrackID() == 1) &&
        (aStep->GetTrack()->GetParentID()== 0)) ||
        (aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
        GetParticleName() ==
        F03PrimaryGeneratorAction::GetPrimaryName())) 
        &&
        (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber")
        &&
        (aStep->GetTrack()->GetNextVolume()->GetName()=="World") &&
        (aStep->GetPostStepPoint()->GetProcessDefinedStep()
              ->GetProcessName() == "Transportation") &&
        (aStep->GetTrack()->GetMomentumDirection().z()>0.)
                                                        )
     {
       eventaction->SetTr();
       Theta = std::acos(aStep->GetTrack()->GetMomentumDirection().x()) ;
       runaction->FillTh(Theta) ;
       Ttrans = aStep->GetTrack()->GetKineticEnergy() ;
       runaction->FillTt(Ttrans) ;
       yend= aStep->GetTrack()->GetPosition().y() ;
       xend= aStep->GetTrack()->GetPosition().x() ;
       rend = std::sqrt(yend*yend+xend*xend) ;
       runaction->FillR(rend);
     }
       
  if (
      (((aStep->GetTrack()->GetTrackID() == 1) &&
        (aStep->GetTrack()->GetParentID()== 0)) ||
      (aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
       GetParticleName() ==
       F03PrimaryGeneratorAction::GetPrimaryName())) 
       &&
      (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber") &&
      (aStep->GetTrack()->GetNextVolume()->GetName()=="World") &&
      (aStep->GetPostStepPoint()->GetProcessDefinedStep()
               ->GetProcessName() == "Transportation") &&
      (aStep->GetTrack()->GetMomentumDirection().z()<0.)
                                                        )
     {
       eventaction->SetRef();
       Thetaback = std::acos(aStep->GetTrack()->GetMomentumDirection().x()) ;
       Thetaback -= 0.5*pi ;
       runaction->FillThBack(Thetaback) ;
       Tback  = aStep->GetTrack()->GetKineticEnergy() ;
       runaction->FillTb(Tback) ;
     }
 

  if (
      ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber") &&
      (aStep->GetTrack()->GetNextVolume()->GetName()=="World") &&
      (aStep->GetPostStepPoint()->GetProcessDefinedStep()
               ->GetProcessName() == "Transportation") &&
      (aStep->GetTrack()->GetMomentumDirection().x()>0.) &&
      (aStep->GetTrack()->GetDynamicParticle()->GetDefinition()
       ->GetParticleName() == "gamma"))
     )
     {
       Egamma = aStep->GetTrack()->GetKineticEnergy() ;
       runaction->FillGammaSpectrum(Egamma) ;
     }
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

