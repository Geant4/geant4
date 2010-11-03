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
// $Id: Test2RunPS.cc,v 1.1 2010-11-03 08:48:57 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Test2RunPS.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

#include <fstream>

Test2RunPS::Test2RunPS(G4String& detName,std::vector<G4String>& hcnameVec) {

  G4SDManager * SDMan = G4SDManager::GetSDMpointer();
  SDMan->SetVerboseLevel(1);
  DTName = detName;
  G4String fullName;
  for(size_t i = 0; i < hcnameVec.size() ; i++) {
    fullName = detName+"/"+hcnameVec[i];
    G4int id =  SDMan->GetCollectionID(fullName);
    if ( id >= 0 ) {
      HCName.push_back(fullName);
      CollID.push_back(id);
      MapSum.push_back(new G4THitsMap<G4double>(detName,hcnameVec[i]));
    }
  }
}

Test2RunPS::~Test2RunPS() {
  //--- Clear HitsMap for RUN                                                   
  G4int Nmap = MapSum.size();
  for ( G4int i = 0; i < Nmap; i++){
    if(MapSum[i] ) MapSum[i]->clear();
  }
  HCName.clear();
  CollID.clear();
  MapSum.clear();
}

#include "Test2PhantomHit.hh"

void Test2RunPS::RecordEvent(const G4Event* evt) {

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;

  G4int Nmap = MapSum.size();  
  for(G4int i = 0; i < Nmap; i++) {
    G4THitsMap<G4double>* evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(CollID[i]));
    *MapSum[i] += *evtMap;
  }
}

void Test2RunPS::DumpQuantitiesToFile(){
  G4int Nmap = MapSum.size();  
  for(G4int i = 0; i < Nmap; i++) {
    DumpQuantitiyToFile(i);
  }
}

void     Test2RunPS::DumpQuantitiyToFile(G4int i){
  if ( MapSum.size() > size_t(i+1) && MapSum[i] ) {
    G4String filename = HCName[i];
    DumpQuantitiyToFile(MapSum[i],filename);
  }
}

void     Test2RunPS::DumpQuantitiyToFile(const G4THitsMap<G4double> *map, 
					 G4String& fileName){
  std::ofstream ofile(fileName);
  if ( ofile ){
    std::map<G4int,G4double*>::iterator itr = map->GetMap()->begin();
    for(; itr != map->GetMap()->end(); itr++) {
      ofile << (itr->first) << "\t" << *(itr->second) << G4endl;
    }
  }
  ofile.close();
}

G4double Test2RunPS::GetTotal(G4int i) const{
  if ( MapSum.size() > size_t(i+1) && MapSum[i] ) {
    return GetTotal(MapSum[i]);
  }else{
    return 0;
  }
}

G4double Test2RunPS::GetTotal(const G4THitsMap<G4double> *map) const
{
  G4double total = 0.;
  std::map<G4int,G4double*>::iterator itr = map->GetMap()->begin();
  for(; itr != map->GetMap()->end(); itr++) {
    total += *(itr->second);
  }
  return total;
}




  
