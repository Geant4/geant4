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
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17DetectorConstruction.hh"
#include "Test17SteppingAction.hh"
#include "Test17EventAction.hh"
#include "Test17RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17SteppingAction::Test17SteppingAction(Test17DetectorConstruction* DET,
                                           Test17EventAction* EA,
                                           Test17RunAction* RA)
:detector(DET),eventaction (EA),runaction(RA),
 IDnow(-2),IDold(-1),prim(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17SteppingAction::~Test17SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  G4double Edep,Theta,xend,Tkin;

  G4double Tsec = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4int evno = eventaction->GetEventNo() ; 

  IDnow = evno+10000*(aStep->GetTrack()->GetTrackID())+
          100000000*(aStep->GetTrack()->GetParentID()); 

  Tkin  = aStep->GetTrack()->GetKineticEnergy(); 
  Edep  = aStep->GetTotalEnergyDeposit();         
  Theta = acos(aStep->GetTrack()->GetMomentumDirection().x());
  xend  = aStep->GetPostStepPoint()->GetPosition().x()/mm;

  eventaction->AddE(Edep);

  // new particle
  if(IDnow != IDold) {
    IDold=IDnow ;
    if(IDnow == 10000) {
      runaction->FillEn(Tsec);
      runaction->FillDef(aStep->GetTrack()->
                         GetDynamicParticle()->GetDefinition());
    }

    // primary
    prim  = false;
    if(1 == aStep->GetTrack()->GetTrackID() ) {
      prim = true;

    // secondary
    } else {
      if(((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
            GetParticleName()) == "e-") 
                 ||
         ((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
            GetParticleName()) == "e+") ) {

            eventaction->AddCharged();

      } else if ((aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
		  GetParticleName()) == "gamma") {
            eventaction->AddNeutral();
      }
    }
  }

    // Primary 
  if( prim ) {
    eventaction->CountStepsCharged(aStep->GetStepLength());

    // Stopping or end of track (only once per event)
    if(0.0 == Tkin || 0.0 > xend || xend >= 100.0*mm) {

      if(0.5*(eventaction->TrackLength()) > xend) {

        G4cout << "EvtNo= " << evno
               << "  length= " << eventaction->TrackLength() 
               << "  xend= " << xend << G4endl;
        eventaction->CountEvent(false);

      } else {

        if(eventaction->EventVerbose() > 0) {
          G4cout << "EvtNo= " << evno
                 << "  length= " << eventaction->TrackLength() 
                 << "  xend= " << xend << G4endl;
	}
        eventaction->CountEvent(true);
        runaction->EndOfTrackCharged(xend) ;
      }
      prim = false;
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





