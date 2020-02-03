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
#include "CCalEventAction.hh"
#include "CCaloSD.hh"
#include "CCalPrimaryGeneratorAction.hh"
#include "CCalSteppingAction.hh"
#include "CCalG4HitCollection.hh"
#include "CCalG4Hit.hh"
#include "CCaloOrganization.hh"
#include "CCalSDList.hh"
#include "CCalAnalysis.hh"


#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4UserSteppingAction.hh"
#include "G4UImanager.hh"

#include <iostream>
#include <vector>
#include <map>

#include "CCalAnalysis.hh"

//#define debug
//#define ddebug


CCalEventAction::CCalEventAction(CCalPrimaryGeneratorAction* pg,
                                 CCalSteppingAction* sa): 
  isInitialized(false),fPrimaryGenerator(pg),
  fSteppingAction(sa),SDnames(nullptr),numberOfSD(0) 
{  
#ifdef debug
  G4cout << "Instantiate CCalEndOfEventAction" << G4endl;
#endif 
  
  G4cout << "Get Calorimter organisation" << G4endl;
  theOrg = new CCaloOrganization;
}

CCalEventAction::~CCalEventAction() {
  delete theOrg;
  delete[] SDnames;
}

void CCalEventAction::initialize() {

  isInitialized = true;
  numberOfSD = CCalSDList::getInstance()->getNumberOfCaloSD();
#ifdef debug
  G4cout << "CCalEndOfEventAction look for " << numberOfSD 
       << " calorimeter-like SD" << G4endl;
#endif
  if (numberOfSD > 0) {
    G4int n = numberOfSD;
    n = std::min(n, 2);
    SDnames = new nameType[n];
  }
  for (G4int i=0; i<numberOfSD; ++i) {
    SDnames[i] = G4String(CCalSDList::getInstance()->getCaloSDName(i));
#ifdef debug
    G4cout << "CCalEndOfEventAction: found SD " << i << " name "
         << SDnames[i] << G4endl;
#endif
  }       
}

void CCalEventAction::BeginOfEventAction(const G4Event* evt) { 

  if (!isInitialized) { initialize(); }
  G4cout << " --- Begin of event: " << evt->GetEventID() << G4endl;
  /*
  if(15 == evt->GetEventID()) {
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/tracking/verbose 2");
  }
  */
}

void CCalEventAction::EndOfEventAction(const G4Event* evt){

  fSteppingAction->endOfEvent();  
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
  G4float hcalE[28], ecalE[49], fullE=0., edec=0, edhc=0;
  G4int i = 0;
  for (i = 0; i < 28; i++) {hcalE[i]=0.;}
  for (i = 0; i < 49; i++) {ecalE[i]=0.;}

  G4float* edep = new G4float[numberOfSD];
  G4int nhit=0;
  for (i = 0; i < numberOfSD; ++i){

    //
    // Look for the Hit Collection
    //
    edep[i] = 0;
    G4int caloHCid = G4SDManager::GetSDMpointer()->GetCollectionID(SDnames[i]);

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
  
          G4int j;
          for (j=0; j<nentries; j++){
#ifdef ddebug
            G4cout << "Hit " << j;
#endif
            CCalG4Hit* aHit =  (*theHC)[j];
            G4float En = aHit->getEnergyDeposit();
            G4int unitID = aHit->getUnitID();
            G4int id=-1;
            if (unitID > 0 && unitID < 29) {
              id = unitID - 1; // HCal
              hcalE[id] += En/GeV;
            } else {
              G4int i0 = unitID/4096;
              G4int i1 = (unitID/64)%64;
              G4int i2 = unitID%64;
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

  G4ThreeVector pos = fPrimaryGenerator->GetParticlePosition();
  G4float ener = fPrimaryGenerator->GetParticleEnergy()/GeV;
  G4float x    = pos.x()/mm;
  G4float y    = pos.y()/mm;
  G4float z    = pos.z()/mm;

  //Save results
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  //1)
  static G4int IDenergy = -1;
  if (IDenergy < 0)
    IDenergy = man->GetH1Id("h4000");
  man->FillH1(IDenergy,fullE); 
  //2)
  G4double totalFilledEnergyHcal = 0.0;
  static G4int IDhcalE = -1;
  if (IDhcalE < 0)
    IDhcalE = man->GetH1Id("h100");
  for (G4int j=0; j<28; j++) {
    man->FillH1(IDhcalE+j,hcalE[j]);
#ifdef debug
    G4cout << "Fill Hcal histo " << j << " with " << hcalE[j] << G4endl;
#endif    
    totalFilledEnergyHcal += hcalE[j];  
  }
#ifdef debug
    G4cout << 
      "CCalAnalysis::InsertEnergyHcal: Total filled Energy Hcal histo " 
           << totalFilledEnergyHcal << G4endl;
#endif

    //3)
    static G4int IDecalE = -1;
    if (IDecalE < 0)
      IDecalE = man->GetH1Id("h200");
    G4double totalFilledEnergyEcal = 0.0;
  for (G4int j=0; j<49; j++) {
    man->FillH1(IDecalE+j,ecalE[j]);
#ifdef debug
    G4cout << "Fill Ecal histo " << j << " with " << ecalE[j] << G4endl;
#endif    
    totalFilledEnergyEcal += ecalE[j];  
  }
#ifdef debug
  G4cout << 
    "CCalAnalysis::InsertEnergyEal: Total filled Energy Ecal histo " 
         << totalFilledEnergyEcal << G4endl;
#endif
  // 4)
  G4int counter=0;
  for (G4int j=0; j<28; j++) 
    {
      man->FillNtupleFColumn(counter,hcalE[j]);
      counter++;
    }
  for (G4int j=0; j<49; j++) 
    {
      man->FillNtupleFColumn(counter,ecalE[j]);
      counter++;
    }
  man->FillNtupleFColumn(counter,ener);
  man->FillNtupleFColumn(counter+1,x);
  man->FillNtupleFColumn(counter+2,y);
  man->FillNtupleFColumn(counter+3,z);
  man->FillNtupleFColumn(counter+4,fullE);
  man->FillNtupleFColumn(counter+5,edec);
  man->FillNtupleFColumn(counter+6,edhc);
  man->AddNtupleRow();  
#ifdef debug
  G4cout << "CCalAnalysis:: Fill Ntuple " << G4endl;
#endif
    
  // 5)
  static G4int IDtimeProfile = -1;
  if (IDtimeProfile < 0)
    IDtimeProfile = man->GetH1Id("h901");
  for (i = 0; i < numberOfSD; i++){
    G4int caloHCid = G4SDManager::GetSDMpointer()->GetCollectionID(SDnames[i]);
    if (caloHCid >= 0) {
      CCalG4HitCollection* theHC = 
        (CCalG4HitCollection*) allHC->GetHC(caloHCid);
      if (theHC != 0) {
        G4int nentries = theHC->entries();
        if (nentries > 0) {
          for (G4int k=0; k<nentries; k++) {
            CCalG4Hit* aHit =  (*theHC)[k];
            man->FillH1(IDtimeProfile,aHit->getTimeSlice(),aHit->getEnergyDeposit()/GeV);

#ifdef debug
            G4cout << "CCalAnalysis:: Fill Time Profile with Hit " << k
                   << " Edeposit " << edep << " Gev" << G4endl;
#endif
          }
        }
      }
    }
  }
}


