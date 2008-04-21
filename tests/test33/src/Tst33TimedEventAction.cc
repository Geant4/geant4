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
// $Id: Tst33TimedEventAction.cc,v 1.14 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

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

Tst33TimedEventAction::Tst33TimedEventAction(G4int time)
  :
  fCScorer(0),
  fProcessTime(0.),
  fMaxRunTime(time),
  fOld_lwe(0)
{
}


Tst33TimedEventAction::~Tst33TimedEventAction()
{
  fOut.close();
}

void Tst33TimedEventAction::Clear() {
  fCScorer = 0;
}
void Tst33TimedEventAction::SpecialCellScorer(const G4CellScorer *scorer){
  fCScorer = scorer;
  fSig.Init();
  fOld_lwe = 0;
}


void Tst33TimedEventAction::BeginOfEventAction(const G4Event*)
{
  fTimer.Start();
}


void Tst33TimedEventAction::EndOfEventAction(const G4Event*)
{
  fTimer.Stop();

  fProcessTime += fTimer.GetUserElapsed()*100;
  if (fCScorer) {
    G4CellScoreValues v=fCScorer->GetCellScoreValues();
    G4double lwe(v.fSumSLWE);
    G4double lwediff(lwe - fOld_lwe); // energy weighted step length 
                                      // in the last event
    fOld_lwe = lwe;

    fSig.Xin(lwediff);

  }
  if (fProcessTime >= fMaxRunTime) {
    G4RunManager::GetRunManager()->AbortRun();
    G4cout << "Tst33TimedEventAction::EndOfEventAction: aborted after"
	   << " running for: " << fProcessTime << " (mu)seconds" 
	   << G4endl;
    CalculateFOM();
    fSig.Init();
    fProcessTime = 0.;
    fOld_lwe = 0;
  }  
}

void Tst33TimedEventAction::CalculateFOM(){
  G4int entries=fSig.GetEntries();
  G4double fom=0;
  G4double mean=fSig.GetMean();
  G4double R=-1;
  G4double sigma=fSig.GetSigma();

  if (sigma>0. && mean > 0.) {
    R = sigma/(mean*std::sqrt(static_cast<G4double>(entries)));
    fom = 1./(R*R);
  }
  G4cout << "Tst33TimedEventAction::CalculateFOM(): FOM: " 
		<< G4endl;
  G4cout << "    mean= " << mean 
		<< ", error= " << R 
		<< ", fom*time= " << fom
		<< ", time= " << fProcessTime
		<< ", n= " << fSig.GetEntries() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

