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
// $Id: Test2RunSD.cc,v 1.1 2010-11-03 08:48:57 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Test2RunSD.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

#define NPRIM 7

Test2RunSD::Test2RunSD(const G4String& detName, const G4String& hcname,
		       std::vector<G4String>& hcnameVec) {

  G4SDManager * SDMan = G4SDManager::GetSDMpointer();
  SDMan->SetVerboseLevel(1);
  //
  G4String fullname = detName+"/"+hcname;
  fSdID =  SDMan->GetCollectionID(fullname);
  if ( fSdID >=0 ) HitSum = new Test2SDHitSum(detName,hcnameVec);
  G4cout << "ASO " << fullname <<" "<<fSdID<< G4endl;
}

Test2RunSD::~Test2RunSD() {
  if ( HitSum ) delete HitSum;
}

#include "Test2PhantomHit.hh"

void Test2RunSD::RecordEvent(const G4Event* evt) {

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;

  if ( fSdID >= 0 ){
    Test2PhantomHitsCollection * phantomHC = 
      (Test2PhantomHitsCollection*)(HCE->GetHC(fSdID));
    if(phantomHC) {
      G4int nent = phantomHC->entries();
      for(G4int i = 0; i < nent; i++) {
	Test2PhantomHit * sdHit = (Test2PhantomHit*)(phantomHC->GetHit(i));
	HitSum->Analyze(sdHit);
      }
    }
  }
}

void Test2RunSD::DumpQuantitiesToFile(){
  if ( fSdID >= 0 ){
    HitSum->DumpQuantitiesToFile();
  }
}

G4double Test2RunSD::GetTotal(G4int i) const
{
  if ( fSdID >= 0 ){
    G4double total = 0.;
    total = HitSum->GetTotal(i);
    return total;
  }
  return 0.0;
}

