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
/// \file B3/B3a/include/EventAction.hh
/// \brief Definition of the B3a::EventAction class

#ifndef B3aEventAction_h
#define B3aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

namespace B3a
{

class RunAction;

/// Event action class
///
/// In EndOfEventAction() there is collected information event per event
/// from Hits Collections, and accumulated statistic for
/// RunAction::EndOfRunAction().

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    ~EventAction() override = default;

    void  BeginOfEventAction(const G4Event*) override;
    void    EndOfEventAction(const G4Event*) override;

  private:
    RunAction* fRunAction = nullptr;
    G4int fCollID_cryst = -1;
    G4int fCollID_patient = -1;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


