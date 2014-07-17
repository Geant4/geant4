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
// $Id: Tst69EventAction.cc 66512 2012-12-19 10:26:52Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////
//
// Tst69EventAction
//
// Created: 17.07.2014 D. Mancusi
//
// Modified:
// 17.07.2014 Adapted from test46 (D. Mancusi)
//
////////////////////////////////////////////////////////////////////////
// 

#include "Tst69EventAction.hh"
#include "G4Event.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst69EventAction::Tst69EventAction():
  actionStart(clock()),
  currentEventStart(0),
  lastEventStart(0),
  lastEventPrint(0),
  lastEventIDPrint(0),
  eventPrintModulo(50),
  eventID(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst69EventAction::~Tst69EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst69EventAction::BeginOfEventAction(const G4Event* anEvent)
{
  currentEventStart = clock();
  eventID = anEvent->GetEventID();
  const G4int secondsSinceLastPrint = (currentEventStart-lastEventPrint)/CLOCKS_PER_SEC;
  bool forcePrint = false;
  if(secondsSinceLastPrint>nSecondsReport) {
    G4cout << "Tst69EventAction: event " << eventID << " starting after more than "
      << nSecondsReport << " s since last printout ("
      << secondsSinceLastPrint
      << " s), printout forced"
      << G4endl;
    forcePrint = true;
  }
  if(forcePrint || (eventID-lastEventIDPrint)%eventPrintModulo==0) {
    G4cout << "Tst69EventAction: starting event " << eventID << ", "
      << (currentEventStart-actionStart)/CLOCKS_PER_SEC
      << " s after action creation, "
      << secondsSinceLastPrint
      << " s after last printout. "
      << G4endl;
    lastEventIDPrint = eventID;
    lastEventPrint = currentEventStart;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst69EventAction::EndOfEventAction(const G4Event*)
{
  const G4int secondsSinceEventStart = (clock()-currentEventStart)/CLOCKS_PER_SEC;
  if(secondsSinceEventStart>nSecondsPerEventMax) {
    G4cout << "Tst69EventAction: event " << eventID
      << " completed after more than " << nSecondsPerEventMax
      << " s (" << secondsSinceEventStart << " s)"
      << G4endl;
  }
  lastEventStart = currentEventStart;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
