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
// $Id: Tst33TimedEventAction.cc,v 1.2 2002-10-31 08:32:45 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst33TimedEventAction.hh"

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

Tst33TimedEventAction::Tst33TimedEventAction(G4int time)
  :
  fCScorer(0),
  fEvStartTime(0),
  fProcessTime(0),
  //  fOut("FluxWeightedEnergy.txt",G4std::ios::out|G4std::ios::app),
  fMaxRunTime(time),
  fOld_lwe(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst33TimedEventAction::~Tst33TimedEventAction()
{
  fOut.close();
}

void Tst33TimedEventAction::Clear() {
  fCScorer = 0;
}
void Tst33TimedEventAction::SetCell_19_Scorer(const G4CellScorer *scorer){
  fCScorer = scorer;
  fSig.Init();
  fOld_lwe = 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst33TimedEventAction::BeginOfEventAction(const G4Event* evt)
{
  struct G4std::tms time = {0};
  G4std::times(&time);
  fEvStartTime = time.tms_utime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst33TimedEventAction::EndOfEventAction(const G4Event* evt)
{

  struct G4std::tms time = {0};
  G4std::times(&time);
  fProcessTime += time.tms_utime - fEvStartTime;
  if (fCScorer) {
    G4CellScoreValues v=fCScorer->GetCellScoreValues();
    G4double lwe(v.fSumSLWE);
    G4double lwediff(lwe - fOld_lwe); // energy weighted step length 
                                      // in the last event
    fOld_lwe = lwe;

    fSig.Xin(lwediff);
    //    G4std::G4cout << "fSumSLWE= " << v.fSumSLWE 
    //	   << " fSumSLW= " << v.fSumSLW
    //	   << " pro.time= " << fProcessTime
    //	   << G4endl;
  }
  if (fProcessTime >= fMaxRunTime) {
    G4RunManager::GetRunManager()->AbortRun();
    G4std::G4cout << "Tst33TimedEventAction::EndOfEventAction: aborted after"
	   << " running for: " << fProcessTime << " (mu)seconds" 
	   << G4endl;
    CalculateFOM();
    fSig.Init();
    fProcessTime = 0;
    fOld_lwe = 0;
  }  
  
}

void Tst33TimedEventAction::CalculateFOM(){
  G4int entries(fSig.GetEntries());
  G4double fom(0);
  G4double mean(fSig.GetMean());
  G4double R(-1);
  G4double sigma(fSig.GetSigma());

  if (sigma>0. && mean > 0.) {
    R = sigma/(mean*G4std::sqrt(entries));
    fom = 1./(R*R);
  }
  G4std::G4cout << "Tst33TimedEventAction::CalculateFOM(): FOM: " 
		<< G4endl;
  G4std::G4cout << "    mean= " << mean 
		<< ", error= " << R 
		<< ", fom*time= " << fom
		<< ", time= " << fProcessTime
		<< ", n= " << fSig.GetEntries() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

