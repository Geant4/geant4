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
// $Id: TiaraTimedEventAction.cc,v 1.3 2005/12/15 14:24:33 ahoward Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TiaraTimedEventAction.hh"

#include <sys/times.h>
#include <time.h>



#include "G4Event.hh"
#include "G4EventManager.hh"
#include "TiaraCellScorerStore.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TiaraTimedEventAction::TiaraTimedEventAction(G4int time):
  fScorerStore(0),
  fEvStartTime(0),
  fCurrentRunTime(0),
  fMaxRunTime(time),
  fTimeFromPreviousRuns(0),
  fFilename("none")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TiaraTimedEventAction::~TiaraTimedEventAction() 
{}

void TiaraTimedEventAction::BeginOfEventAction(const G4Event*)
{
  struct tms time;
  times(&time);
  fEvStartTime = time.tms_utime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TiaraTimedEventAction::
SetScorerStore(TiaraCellScorerStore *scorerStore){
  fScorerStore = scorerStore;
}

void TiaraTimedEventAction::EndOfEventAction(const G4Event* )
{
  
  if (fScorerStore) {
    fScorerStore->EndOfEventAction();
  }

  struct tms time;
  times(&time);


  G4int dt(time.tms_utime - fEvStartTime);
  fCurrentRunTime += dt;

  if (fCurrentRunTime >= fMaxRunTime) {
    G4cout << "TiaraTimedEventAction::EndOfEventAction: aborting after"
	   << " running for: " << fCurrentRunTime << " (mu)seconds" 
	   << G4endl;
    G4RunManager::GetRunManager()->AbortRun();
    fTimeFromPreviousRuns+=fCurrentRunTime;
    fCurrentRunTime = 0;
    if (fFilename!="none") {
      CLHEP::HepRandom::saveEngineStatus(fFilename);
    }
  }  


}  

G4int TiaraTimedEventAction::GetTotalProcessedTime() const {
  return fTimeFromPreviousRuns + fCurrentRunTime;
}

void TiaraTimedEventAction::SetRnadomNumFilename(const G4String &fname) {
  fFilename = fname;
}
