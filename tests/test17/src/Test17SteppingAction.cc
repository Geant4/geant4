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
  G4double Edep,xend,Tkin;

  G4double Tsec = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4int evno = eventaction->GetEventNo(); 
  G4int runo = runaction->RunID(); 
  G4int trno = aStep->GetTrack()->GetTrackID(); 

  IDnow = trno + evno*10000 + 100000000*runo; 

  Tkin  = aStep->GetTrack()->GetKineticEnergy(); 
  Edep  = aStep->GetTotalEnergyDeposit();         
  xend  = aStep->GetPostStepPoint()->GetPosition().x()/mm;

  eventaction->AddE(Edep);

  // new particle
  if(IDnow != IDold) {
    IDold=IDnow ;
    if(trno == 1) {
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





