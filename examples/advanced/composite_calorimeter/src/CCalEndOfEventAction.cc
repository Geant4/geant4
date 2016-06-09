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
///////////////////////////////////////////////////////////////////////////////
// File: CCalEndOfEventAction.cc
// Description: CCalEndOfEventAction provides User actions at end of event
///////////////////////////////////////////////////////////////////////////////
#include "CCalEndOfEventAction.hh"
#include "CCaloSD.hh"
#include "CCalPrimaryGeneratorAction.hh"
#include "CCalG4HitCollection.hh"
#include "CCalG4Hit.hh"
#include "CCaloOrganization.hh"
#include "CCalSDList.hh"
#include "CCalSteppingAction.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4UserSteppingAction.hh"

#include <iostream>
#include <vector>
#include <map>

#ifdef G4ANALYSIS_USE
#include "CCalAnalysis.hh"
#endif

//#define debug
//#define ddebug


CCalEndOfEventAction::CCalEndOfEventAction (CCalPrimaryGeneratorAction* pg): 
  isInitialized(false),SDnames(0),numberOfSD(0) {

  primaryGenerator = pg;
#ifdef debug
  G4cout << "Instantiate CCalEndOfEventAction" << G4endl;
#endif

  G4cout << "Now Instantiate stepping action" << G4endl;
  instanciateSteppingAction();
  
  G4cout << "Get Calorimter organisation" << G4endl;
  theOrg = new CCaloOrganization;
  G4cout << "end of instantiation of EndofEventAction" << G4endl;
}


CCalEndOfEventAction::~CCalEndOfEventAction() {

  if (theOrg)
    delete theOrg;
  if (SDnames)
    delete[] SDnames;
  G4cout << "CCalEndOfEventAction deleted" << G4endl;
}


void CCalEndOfEventAction::initialize() {

  isInitialized = true;
  numberOfSD = CCalSDList::getInstance()->getNumberOfCaloSD();
#ifdef debug
  G4cout << "CCalEndOfEventAction look for " << numberOfSD 
       << " calorimeter-like SD" << G4endl;
#endif
  if (numberOfSD > 0) {
    G4int n = numberOfSD;
    if (n < 2) n = 2;
    SDnames = new nameType[n];
  }
  for (int i=0; i<numberOfSD; i++) {
    SDnames[i] = G4String(CCalSDList::getInstance()->getCaloSDName(i));
#ifdef debug
    G4cout << "CCalEndOfEventAction: found SD " << i << " name "
	 << SDnames[i] << G4endl;
#endif
  }       
}


void CCalEndOfEventAction::StartOfEventAction(const G4Event* evt) { 
  G4cout << "\n---> Begin of event: " << evt->GetEventID() << G4endl;
}


void CCalEndOfEventAction::EndOfEventAction(const G4Event* evt){

#ifdef debug
  G4cout << G4endl << "=== Begin of EndOfEventAction === " << G4endl;
#endif

  if (!isInitialized) initialize();

  theSteppingAction->endOfEvent();
  
  //
  // Look for the Hit Collection 
  //  
  G4HCofThisEvent* allHC = evt->GetHCofThisEvent();
  if (allHC == 0) {
#ifdef debug
    G4cout << "CCalEndOfEventAction: No Hit Collection in this event" 
	 << G4endl;
#endif
    return;
  }
  	
  //
  // hits info
  //
  
  //Now make summary
  float hcalE[28], ecalE[49], fullE=0., edec=0, edhc=0;
  int i = 0;
  for (i = 0; i < 28; i++) {hcalE[i]=0.;}
  for (i = 0; i < 49; i++) {ecalE[i]=0.;}

  float* edep = new float[numberOfSD];
  int nhit=0;
  for (i = 0; i < numberOfSD; i++){

    //
    // Look for the Hit Collection
    //
    edep[i] = 0;
    int caloHCid = G4SDManager::GetSDMpointer()->GetCollectionID(SDnames[i]);

    if (caloHCid >= 0) {
      CCalG4HitCollection* theHC = 
	(CCalG4HitCollection*) allHC->GetHC(caloHCid);
    
      if (theHC != 0) {

	G4int nentries = theHC->entries();
#ifdef debug
	G4cout << " There are " << nentries << " hits in " << SDnames[i] 
	       << " :" << G4endl;
#endif

	if (nentries > 0) {
  
	  int j;
	  for (j=0; j<nentries; j++){
#ifdef ddebug
	    G4cout << "Hit " << j;
#endif
	    CCalG4Hit* aHit =  (*theHC)[j];
	    float En = aHit->getEnergyDeposit();
	    int unitID = aHit->getUnitID();
	    int id=-1;
	    if (unitID > 0 && unitID < 29) {
	      id = unitID - 1; // HCal
	      hcalE[id] += En/GeV;
	    } else {
	      int i0 = unitID/4096;
	      int i1 = (unitID/64)%64;
	      int i2 = unitID%64;
	      if (i0 == 1 && i1 < 8 && i2 < 8) {
		id = i1*7 + i2; // ECal
		ecalE[id] += En/GeV;
	      }
	    }
#ifdef ddebug
	    G4cout << " with Energy = " << En/MeV << " MeV in Unit " << unitID 
		   << " " << id << G4endl;
#endif
	    fullE   += En/GeV;
	    edep[i] += En/GeV;
	    nhit++;
	  }
#ifdef ddebug
	  G4cout << " ===> Total Energy Deposit in this Calorimeter = " 
		 << edep[i]*1000.0 << "  MeV " << G4endl; 
#endif
	}
      }
    }
    if (SDnames[i] == "HadronCalorimeter") {
      edhc = edep[i];
    } else if (SDnames[i] == "CrystalMatrix") {
      edec = edep[i];
    }
  }

  delete[] edep;

#ifdef G4ANALYSIS_USE
  G4ThreeVector pos = primaryGenerator->GetParticlePosition();
  float ener = primaryGenerator->GetParticleEnergy()/GeV;
  float x    = pos.x()/mm;
  float y    = pos.y()/mm;
  float z    = pos.z()/mm;

  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->InsertEnergy(fullE);
  analysis->InsertEnergyHcal(hcalE);
  analysis->InsertEnergyEcal(ecalE);
  analysis->setNtuple(hcalE, ecalE, ener, x, y, z, fullE, edec, edhc);
  analysis->EndOfEvent(nhit);
  for (i = 0; i < numberOfSD; i++){
    int caloHCid = G4SDManager::GetSDMpointer()->GetCollectionID(SDnames[i]);
    if (caloHCid >= 0) {
      CCalG4HitCollection* theHC = 
	(CCalG4HitCollection*) allHC->GetHC(caloHCid);
      if (theHC != 0) {
	G4int nentries = theHC->entries();
	if (nentries > 0) {
	  for (G4int k=0; k<nentries; k++) {
	    CCalG4Hit* aHit =  (*theHC)[k];
	    analysis->InsertTimeProfile(k,aHit->getTimeSlice(),
					aHit->getEnergyDeposit()/GeV);
	  }
	}
      }
    }
  }
#endif
}


void CCalEndOfEventAction::instanciateSteppingAction(){
	
  G4UserSteppingAction* theUA = const_cast<G4UserSteppingAction*>(G4RunManager::GetRunManager()->GetUserSteppingAction());
        
  if (theUA == 0) {
#ifdef debug
    G4cout << " CCalEndOfEventAction::instanciateSteppingAction creates"
	 << " CCalSteppingAction" << G4endl;
#endif
    theSteppingAction = new CCalSteppingAction;  
    G4RunManager::GetRunManager()->SetUserAction(theSteppingAction);
  }   
	
}
