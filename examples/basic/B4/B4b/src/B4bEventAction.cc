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
// $Id$
// 
/// \file B4bEventAction.cc
/// \brief Implementation of the B4bEventAction class

#include "B4bEventAction.hh"
#include "B4bRunData.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bEventAction::B4bEventAction()
 : G4UserEventAction(),
   fMessenger(0),
   fPrintModulo(1)
{
  // Define /B4/event commands using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/B4/event/", "Event control");

  // Define /B4/event/setPrintModulo command
  G4GenericMessenger::Command& setPrintModulo 
    = fMessenger->DeclareProperty("setPrintModulo", 
                                  fPrintModulo, 
                                 "Print events modulo n");
  setPrintModulo.SetRange("value>0");                                
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bEventAction::~B4bEventAction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bEventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength) const
{
  // print event statistics
  G4cout
     << "   Absorber: total energy: " 
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: " 
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bEventAction::BeginOfEventAction(const G4Event* evt)
{  

  G4int eventID = evt->GetEventID();
  if ( eventID % fPrintModulo == 0 )  { 
    G4cout << "\n---> Begin of event: " << eventID << G4endl;
    //CLHEP::HepRandom::showEngineStatus();
  }

  B4bRunData::GetInstance()->Reset();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bEventAction::EndOfEventAction(const G4Event* event)
{
  //accumulates statistic
  //
  B4bRunData::GetInstance()->FillPerEvent();
  
  //print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  if ( eventID % fPrintModulo == 0) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    PrintEventStatistics(
      B4bRunData::GetInstance()->GetEdep(kAbs),
      B4bRunData::GetInstance()->GetTrackLength(kAbs),
      B4bRunData::GetInstance()->GetEdep(kGap),
      B4bRunData::GetInstance()->GetTrackLength(kGap));
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
