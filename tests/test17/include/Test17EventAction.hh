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
// Class Description:
// Event action definition.
// Class Description - end
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17EventAction_h
#define Test17EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Test17RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17EventAction : public G4UserEventAction
{
public: // Without description

    Test17EventAction(Test17RunAction* Test17RA);
   ~Test17EventAction();

    void  BeginOfEventAction(const G4Event*);
    void  EndOfEventAction(const G4Event*);

    G4int GetEventNo() const {return evtNo;};
    void  setEventVerbose(G4int level) {verbose = level;};
    G4int EventVerbose() {return verbose;};
    
    void CountStepsCharged(G4double step) {totLength += step;};
    G4double TrackLength() const {return totLength;} ;
    void AddCharged() {nCharged++;};
    void AddNeutral() {nNeutral++;};
    void AddE(G4double e) {edep += e;};
    void CountEvent(G4bool val) {good = val;};
    
  private:
    Test17RunAction* runaction;
    G4int verbose;
    G4int evtNo;
    G4int nCharged;
    G4int nNeutral;
    G4double edep;
    G4double totLength;
    G4bool good;
};

#endif

    
