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
#include "Test2SDHitSum.hh"
#include <fstream>
#include <map>

Test2SDHitSum::Test2SDHitSum(const G4String& detName,
			     std::vector<G4String>& hcnameVec) {
  DTName = detName;
  G4String fullName;
  for(size_t i = 0; i < hcnameVec.size() ; i++) {
      fullName = detName+"/"+hcnameVec[i];
      G4cout << "AAAA" << fullName << G4endl;
      HCName.push_back(fullName);
      MapSum.push_back(new G4THitsMap<G4double>(detName,hcnameVec[i]));
  }
}

Test2SDHitSum::~Test2SDHitSum() {
  //--- Clear HitsMap for RUN                                        
  G4int Nmap = MapSum.size();
  for ( G4int i = 0; i < Nmap; i++){
    if(MapSum[i] ) MapSum[i]->clear();
  }
  HCName.clear();
  MapSum.clear();
}


void Test2SDHitSum::Analyze(Test2PhantomHit* hc){
  G4double one = 1.0;
  G4int index = hc->GetID();
  //
  // eDep
  //
  G4double eDep = hc->GetEdep();
  MapSum[0]->add(index,eDep);

  //
  // trackLength
  //
  G4double tlen = hc->GetTrackLength();
  //
  // trackLengthGamma
  if ( hc->GetParticleName() == "gamma" ){
    MapSum[1]->add(index,tlen);
  }
  //
  // trackLengthElec
  if ( hc->GetParticleName() == "e-" ){
    MapSum[2]->add(index,tlen);
  }
  //
  // trackLengthPosi
  if ( hc->GetParticleName() == "e+" ){
    MapSum[3]->add(index,tlen);
  }
  //
  // nStepGamma
  if ( hc->GetParticleName() == "gamma" ){
    MapSum[4]->add(index,one);
  }
  //
  // nStepElec
  if ( hc->GetParticleName() == "e-" ){
    MapSum[5]->add(index,one);
  }
  //
  // nStepPosi
  if ( hc->GetParticleName() == "e+" ){
    MapSum[6]->add(index,one);
  }

}

void Test2SDHitSum::DumpQuantitiesToFile(){
  G4int Nmap = MapSum.size();  
  G4cout << "Nmap " << Nmap<<G4endl;
  for(G4int i = 0; i < Nmap; i++) {
    DumpQuantitiyToFile(i);
  }
}

void Test2SDHitSum::DumpQuantitiyToFile(G4int i){
  if ( MapSum.size() > size_t(i+1) && MapSum[i] ) {
    G4String filename = HCName[i];
    DumpQuantitiyToFile(MapSum[i],filename);
  }
}

void Test2SDHitSum::DumpQuantitiyToFile(const G4THitsMap<G4double> *map, 
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

G4double Test2SDHitSum::GetTotal(G4int i) const{
  if ( MapSum.size() > size_t(i+1) && MapSum[i] ) {
    return GetTotal(MapSum[i]);
  }else{
    return 0;
  }
}

G4double Test2SDHitSum::GetTotal(const G4THitsMap<G4double>* map) const{
  G4double total = 0.;
  std::map<G4int,G4double*>::iterator itr = map->GetMap()->begin();
  for(; itr != map->GetMap()->end(); itr++) {
    total += *(itr->second);
  }
  return total;
}


