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
// $Id: B01EventAction.cc,v 1.1 2002-10-22 14:09:05 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B01EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "G4CellScorer.hh"
#include "G4RunManager.hh"

#include <sys/times.h>
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B01EventAction::B01EventAction(const G4CellScorer *scorer, G4int time)
  :
  fCScorer(scorer),
  fEvStartTime(0),
  fProcessTime(0),
  //  fOut("FluxWeightedEnergy.txt",G4std::ios::out|G4std::ios::app),
  fMaxRunTime(time)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B01EventAction::~B01EventAction()
{
  fOut.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B01EventAction::BeginOfEventAction(const G4Event* evt)
{
  struct G4std::tms time = {0};
  G4std::times(&time);
  fEvStartTime = time.tms_utime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B01EventAction::EndOfEventAction(const G4Event* evt)
{

  struct G4std::tms time = {0};
  G4std::times(&time);
  fProcessTime += time.tms_utime - fEvStartTime;
  if (fCScorer) {
    G4CellScoreValues v=fCScorer->GetCellScoreValues();
    G4std::G4cout << "fSumSLWE= " << v.fSumSLWE 
	   << " fSumSLW= " << v.fSumSLW
	   << " pro.time= " << fProcessTime
	   << G4endl;
  }
  if (fProcessTime >= fMaxRunTime) {
    G4RunManager::GetRunManager()->AbortRun();
    G4std::G4cout << "B01EventAction::EndOfEventAction: aborted after"
	   << " running for: " << fProcessTime << " (mu)seconds" 
	   << G4endl;
  }  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

