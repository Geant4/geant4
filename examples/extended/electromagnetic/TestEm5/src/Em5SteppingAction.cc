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
// $Id: Em5SteppingAction.cc,v 1.9 2002-06-06 17:23:22 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em5DetectorConstruction.hh"
#include "G4EnergyLossTables.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "Em5SteppingAction.hh"
#include "Em5PrimaryGeneratorAction.hh"
#include "Em5EventAction.hh"
#include "Em5RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "Em5SteppingMessenger.hh"
#include "G4ios.hh"
#include "g4std/iomanip"
#include "G4UImanager.hh"

#ifndef G4NOHIST
 #include "AIDA/IHistogram1D.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em5SteppingAction::Em5SteppingAction(Em5DetectorConstruction* DET,
                                     Em5EventAction* EA,
                                     Em5RunAction* RA)
:detector (DET),eventaction (EA),runaction (RA),steppingMessenger(NULL),
 IDold(-1) ,evnoold(-1)
{
  steppingMessenger = new Em5SteppingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em5SteppingAction::~Em5SteppingAction()
{
  delete steppingMessenger ;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 

#ifndef G4NOHIST
  const G4double cn = pi/(64800.*runaction->GetdTh()) ;
  const G4double cnback = pi/(64800.*runaction->GetdThback()) ;
  G4double wg ;
#endif
  G4double Theta,Thetaback,Ttrans,Tback,Tsec,Egamma,yend,zend,rend;
  G4int evno = eventaction->GetEventno(); 

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
#ifndef G4NOHIST
        if(runaction->GetHisto(8) != 0)
           runaction->GetHisto(8)->fill(
                               aStep->GetTrack()->GetVertexPosition().x());
#endif

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
#ifndef G4NOHIST
        if(runaction->GetHisto(7) != 0)
           runaction->GetHisto(7)->fill(Tsec/MeV) ;
#endif

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
#ifndef G4NOHIST
        if(runaction->GetHisto(7) != 0)
           runaction->GetHisto(7)->fill(Tsec) ;
#endif
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
        (aStep->GetTrack()->GetParentID()== 0))) 
        &&
        (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber")
        &&
        (aStep->GetTrack()->GetNextVolume()->GetName()=="World") &&
        (aStep->GetPostStepPoint()->GetProcessDefinedStep()
              ->GetProcessName() == "Transportation") &&
        (aStep->GetTrack()->GetMomentumDirection().x()>0.)
                                                        )
     {
       eventaction->SetTr();
       Theta = acos(aStep->GetTrack()->GetMomentumDirection().x()) ;
#ifndef G4NOHIST
        if((runaction->GetHisto(2) != 0) && (Theta > 0.))
        {
           wg = cn/sin(Theta) ;
           runaction->GetHisto(2)->fill(Theta/deg,wg) ;
        }
#endif

       Ttrans = aStep->GetTrack()->GetKineticEnergy() ;
#ifndef G4NOHIST
        if(runaction->GetHisto(4) != 0)
           runaction->GetHisto(4)->fill(Ttrans/MeV) ;
#endif
       yend= aStep->GetTrack()->GetPosition().y() ;
       zend= aStep->GetTrack()->GetPosition().z() ;
       rend = sqrt(yend*yend+zend*zend) ;
#ifndef G4NOHIST
        if(runaction->GetHisto(3) != 0)
           runaction->GetHisto(3)->fill(rend) ;
#endif
     }
       
  if (
      (((aStep->GetTrack()->GetTrackID() == 1) &&
        (aStep->GetTrack()->GetParentID()== 0))) 
       &&
      (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber") &&
      (aStep->GetTrack()->GetNextVolume()->GetName()=="World") &&
      (aStep->GetPostStepPoint()->GetProcessDefinedStep()
               ->GetProcessName() == "Transportation") &&
      (aStep->GetTrack()->GetMomentumDirection().x()<0.)
                                                        )
     {
       eventaction->SetRef();
       Thetaback = acos(aStep->GetTrack()->GetMomentumDirection().x()) ;
       Thetaback -= 0.5*pi ;
#ifndef G4NOHIST
        if((runaction->GetHisto(5) != 0) && (Thetaback > 0.))
        {
           wg = cnback/sin(Thetaback) ;
           runaction->GetHisto(5)->fill(Thetaback/deg,wg) ;
        }
#endif
       Tback  = aStep->GetTrack()->GetKineticEnergy() ;
#ifndef G4NOHIST
        if(runaction->GetHisto(6) != 0)
           runaction->GetHisto(6)->fill(Tback/MeV) ;
#endif
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
#ifndef G4NOHIST
        if(runaction->GetHisto(9) != 0)
           runaction->GetHisto(9)->fill(log10(Egamma/MeV)) ;
#endif
     }
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

