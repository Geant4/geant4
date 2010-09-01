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
// $Id: Test2Run.cc,v 1.2 2010-09-01 08:03:10 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Test2Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

#define NPRIM 7

Test2Run::Test2Run() {

  G4String detName = "ScoringWorld";
  G4String primNameSum[NPRIM] = {"eDep",
				 "trackLengthGamma",
				 "nStepGamma",
				 "trackLengthElec",
				 "nStepElec",
				 "trackLengthPosi",
				 "nStepPosi"};
  G4SDManager * SDMan = G4SDManager::GetSDMpointer();
  SDMan->SetVerboseLevel(1);
  G4String fullName;
  for(G4int i = 0; i < NPRIM; i++) {
    fullName = detName+"/"+primNameSum[i];
    colIDSum[i] = SDMan->GetCollectionID(fullName);
  }

  //
  sdID =  SDMan->GetCollectionID(fullName="PhantomCollection");
  sdTotalEdep = 0.;
  sdTotalTrackLengthGamma = 0.;
  sdTotalTrackLengthElec = 0.;
  sdTotalTrackLengthPosi = 0.;
  SDMan->SetVerboseLevel(0);
}

Test2Run::~Test2Run() {;}

#include "Test2PhantomHit.hh"

void Test2Run::RecordEvent(const G4Event* evt) {

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;
  for(G4int i = 0; i < NPRIM; i++) {
    G4THitsMap<G4double>* evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(colIDSum[i]));
    mapSum[i] += *evtMap;
  }

  Test2PhantomHitsCollection * phantomHC = (Test2PhantomHitsCollection*)(HCE->GetHC(sdID));
  if(phantomHC) {
    G4int nent = phantomHC->entries();
    for(G4int i = 0; i < nent; i++) {
      Test2PhantomHit * sdHit = (Test2PhantomHit*)(phantomHC->GetHit(i));
      sdTotalEdep += sdHit->GetEdep();
      if(sdHit->GetParticleName() == "gamma") {
	sdTotalTrackLengthGamma += sdHit->GetTrackLength();
      } else if(sdHit->GetParticleName() == "e-") {
	sdTotalTrackLengthElec += sdHit->GetTrackLength();
      } else if(sdHit->GetParticleName() == "e+") {
	  sdTotalTrackLengthPosi += sdHit->GetTrackLength();
      }
      /*
      G4cout << sdHit->GetEdep() << " at "
	     << sdHit->GetX() << ", " << sdHit->GetY() << ", " << sdHit->GetZ()
	     << G4endl;
      */
    }
  }
}

G4double Test2Run::GetTotal(const G4THitsMap<G4double> &map) const
{
  G4double total = 0.;
  std::map<G4int,G4double*>::iterator itr = map.GetMap()->begin();
  for(; itr != map.GetMap()->end(); itr++) {
    total += *(itr->second);
  }
  return total;
}


  
