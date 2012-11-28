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
#include <sstream>

#include "Tst68EventAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Timer.hh"
#include "G4Event.hh"
#include "G4ios.hh"
#include "G4EventManager.hh"
#include "Tst68SteppingAction.hh"
#include "G4RunManager.hh"

Tst68EventAction::Tst68EventAction() : 
  theSteppingAction( 0 ), numberOfEvents( 0 ), sumTotalDepositedEnergy( 0.0 ) {
  instanciateSteppingAction();
  eventTimer = new G4Timer;
}


Tst68EventAction::~Tst68EventAction() {
  delete eventTimer;
  if ( numberOfEvents > 0 ) {
    G4cout << G4endl
           << " ----------------------------------- " << G4endl
           << " Average deposited energy = " 
           << sumTotalDepositedEnergy / static_cast< G4double >( numberOfEvents )
           << "  MeV " << G4endl
           << " ----------------------------------- " << G4endl
           << G4endl;
  }
}


void Tst68EventAction::BeginOfEventAction( const G4Event* ) {
  //G4cout << "\n---> Begin of event: " << evt->GetEventID() << G4endl;
  eventTimer->Start();
}


void Tst68EventAction::EndOfEventAction( const G4Event* evt ) {

  int incidentParticleId = theSteppingAction->getPrimaryParticleId();
  G4double incidentParticleEnergy = theSteppingAction->getPrimaryParticleEnergy();
  if ( evt->GetEventID() == 0 ) {
    G4cout << "\t ---------------------------------------" << G4endl
           << "\t Beam Particle PDG Id = " << incidentParticleId << G4endl
           << "\t Beam Particle Kinetic Energy = " 
           << incidentParticleEnergy / GeV << " GeV" << G4endl
           << "\t ---------------------------------------" << G4endl;
  }

  G4double totalEdepAllParticles = theSteppingAction->getTotalEdepAllParticles();
  sumTotalDepositedEnergy += totalEdepAllParticles;

  theSteppingAction->reset();

  //G4cout << " ---  Tst68EventAction::EndOfEventAction  ---  event= " 
  G4cout << " ---  EndOfEventAction  ---  event= " 
	 << evt->GetEventID() << "   t=";
  eventTimer->Stop();
  G4double timeEventInSeconds = 0.0;
  if ( eventTimer->IsValid() ) {
    timeEventInSeconds = eventTimer->GetUserElapsed();
    G4cout << timeEventInSeconds << "s";
  } else {
    G4cout << "UNDEFINED";
  }

  G4double r = G4UniformRand();

  union double_ull_t { double d; unsigned long long u;} x; 
  x.d = r;
  G4cout << std::hex << "   random=" << x.u << std::dec << G4endl;

  numberOfEvents++;
}


void Tst68EventAction::instanciateSteppingAction() {      
  G4UserSteppingAction* theUserAction = const_cast< G4UserSteppingAction* >
    ( G4RunManager::GetRunManager()->GetUserSteppingAction() );
  if (theUserAction == 0) {
    theSteppingAction = new Tst68SteppingAction;  
    G4RunManager::GetRunManager()->SetUserAction( theSteppingAction );
  }   
}

