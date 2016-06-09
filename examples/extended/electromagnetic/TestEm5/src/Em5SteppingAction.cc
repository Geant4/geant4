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
// $Id: Em5SteppingAction.cc,v 1.10 2003/04/05 17:33:52 vnivanch Exp $
// GEANT4 tag $Name: geant4-05-01 $
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
  const G4Track* track = aStep->GetTrack();
#ifndef G4NOHIST
  const G4double cn = pi/(64800.*runaction->GetdTh()) ;
  const G4double cnback = pi/(64800.*runaction->GetdThback()) ;
  G4double wg ;
#endif
  G4double Theta,Thetaback,Ttrans,Tback,Tsec,Egamma,yend,zend,rend;
  G4int evno = eventaction->GetEventno(); 

  IDnow = evno+10000*(track->GetTrackID())+
          100000000*(track->GetParentID()); 

  if(IDnow != IDold)
  {
   IDold=IDnow ;
   if(
      (((track->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e-") &&
       ((track->GetTrackID() != 1) ||
       (track->GetParentID() != 0)) )
       ||
      (((track->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e+") &&
       ((track->GetTrackID() != 1) ||
       (track->GetParentID() != 0)) )
     )
#ifndef G4NOHIST
        if(runaction->GetHisto(8) != 0)
           runaction->GetHisto(8)->fill(
                               track->GetVertexPosition().x());
#endif


   if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber")
   {
    if(((track->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e-") &&
       ((track->GetTrackID() != 1) ||
       (track->GetParentID() != 0)) ) 
    {
        eventaction->AddCharged() ;
        eventaction->AddE() ;
        Tsec = track->GetKineticEnergy() ;  // !!!!!!!!!!!!
        Tsec += aStep->GetTotalEnergyDeposit() ;        // !!!!!!!!!!!!
#ifndef G4NOHIST
        if(runaction->GetHisto(7) != 0)
           runaction->GetHisto(7)->fill(Tsec/MeV) ;
#endif

    }
    else
    if(((track->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e+") &&
       ((track->GetTrackID() != 1) ||
       (track->GetParentID() != 0)) ) 
    {
        eventaction->AddCharged() ;
        eventaction->AddP() ;
        Tsec = track->GetKineticEnergy() ;  // !!!!!!!!!!!!
        Tsec += aStep->GetTotalEnergyDeposit() ;        // !!!!!!!!!!!!
#ifndef G4NOHIST
        if(runaction->GetHisto(7) != 0)
           runaction->GetHisto(7)->fill(Tsec) ;
#endif
    }
    else
    if(((track->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "gamma") &&
       ((track->GetTrackID() != 1) ||
       (track->GetParentID() != 0)) ) 
    {
        eventaction->AddNeutral() ;
    }
   }
  }


  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber")
  {
    if(((track->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e-") 
              ||
       ((track->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "e+"))  
          eventaction->CountStepsCharged() ;

    if ((track->GetDynamicParticle()->GetDefinition()->
        GetParticleName()) == "gamma") 
          eventaction->CountStepsNeutral() ;
  }

  G4String nextVolumeName = "OutOfWorld";
  if (track->GetNextVolume()) 
    nextVolumeName = track->GetNextVolume()->GetName();

  if (  nextVolumeName != "OutOfWorld"
        &&
        track->GetTrackID() == 1  
        &&
        aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber"
        &&
        nextVolumeName == "World" 
        &&
        track->GetMomentumDirection().x() > 0. )
     {
       eventaction->SetTr();
       Theta = acos(track->GetMomentumDirection().x()) ;
#ifndef G4NOHIST
        if((runaction->GetHisto(2) != 0) && (Theta > 0.))
        {
           wg = cn/sin(Theta) ;
           runaction->GetHisto(2)->fill(Theta/deg,wg) ;
        }
#endif

       Ttrans = track->GetKineticEnergy() ;
#ifndef G4NOHIST
        if(runaction->GetHisto(4) != 0)
           runaction->GetHisto(4)->fill(Ttrans/MeV) ;
#endif
       yend= track->GetPosition().y() ;
       zend= track->GetPosition().z() ;
       rend = sqrt(yend*yend+zend*zend) ;
#ifndef G4NOHIST
        if(runaction->GetHisto(3) != 0)
           runaction->GetHisto(3)->fill(rend) ;
#endif
     }

       
  if ( nextVolumeName != "OutOfWorld"
       &&
       track->GetTrackID() == 1
       && 
       aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber" 
       &&
       nextVolumeName =="World" 
       &&
       track->GetMomentumDirection().x() < 0. )

     {
       eventaction->SetRef();
       Thetaback = acos(track->GetMomentumDirection().x()) ;
       Thetaback -= 0.5*pi ;
#ifndef G4NOHIST
        if((runaction->GetHisto(5) != 0) && (Thetaback > 0.))
        {
           wg = cnback/sin(Thetaback) ;
           runaction->GetHisto(5)->fill(Thetaback/deg,wg) ;
        }
#endif
       Tback  = track->GetKineticEnergy() ;
#ifndef G4NOHIST
        if(runaction->GetHisto(6) != 0)
           runaction->GetHisto(6)->fill(Tback/MeV) ;
#endif
     }
 
  if ( nextVolumeName != "OutOfWorld"
       &&
       aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Absorber" 
       &&
       track->GetNextVolume()->GetName()=="World" 
       &&
       track->GetMomentumDirection().x()>0. 
       &&
       track->GetDynamicParticle()->GetDefinition()
       ->GetParticleName() == "gamma" )
     
     {
       Egamma = track->GetKineticEnergy() ;
#ifndef G4NOHIST
        if(runaction->GetHisto(9) != 0)
           runaction->GetHisto(9)->fill(log10(Egamma/MeV)) ;
#endif
     }

   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

