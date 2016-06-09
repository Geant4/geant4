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
// $Id: TiaraTimedEventAction.cc,v 1.4 2006/06/29 15:45:40 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
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
