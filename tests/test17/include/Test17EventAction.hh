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

    
