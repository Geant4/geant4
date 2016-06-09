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
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "RunAction.hh"
#include "EventActionMessenger.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
{
  runAct = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
  eventMessenger = new EventActionMessenger(this);
  printModulo = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    CLHEP::HepRandom::showEngineStatus();
}
 
 // initialisation per event
 EnergyAbs = EnergyGap = 0.;
 TrackLAbs = TrackLGap = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  //accumulates statistic
  //
  runAct->fillPerEvent(EnergyAbs, EnergyGap, TrackLAbs, TrackLGap);
  
  //print per event (modulo n)
  //
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    G4cout << "---> End of event: " << evtNb << G4endl;        

    G4cout
       << "   Absorber: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyAbs,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyGap,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLGap,"Length")
       << G4endl;
          
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
