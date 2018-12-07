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
/// \file biasing/B02/src/B02RunAction.cc
/// \brief Implementation of the B02RunAction class
//
//
// 
#include "B02RunAction.hh"
#include "B02Run.hh"

//-- In order to obtain detector information.
#include "G4RunManager.hh"
#include "B02DetectorConstruction.hh"
#include "B02ImportanceDetectorConstruction.hh"
#include "G4THitsMap.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// B02RunAction
//  
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor
B02RunAction::B02RunAction(): 
  G4UserRunAction(),
  //  fFieldName(15),
  fFieldValue(14)
{
  // - Prepare data member for B02Run.
  //   vector represents a list of MultiFunctionalDetector names.
  fSDName.push_back(G4String("ConcreteSD"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Destructor.
B02RunAction::~B02RunAction()
{
  fSDName.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* B02RunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  //  Detail description can be found in B02Run.hh/cc.
  return new B02Run(fSDName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02RunAction::EndOfRunAction(const G4Run* aRun)
{
    G4cout << " ###### EndOfRunAction  " <<G4endl;
  //- B02Run object.
  B02Run* b02Run = (B02Run*)aRun;
  //--- Dump all socred quantities involved in B02Run.
  // re02Run->DumpAllScorer();
  //---
  G4RunManager* mgr = G4RunManager::GetRunManager();
  //
  
  for ( G4int i = 0; i < (G4int)fSDName.size(); i++ ){
    const G4VUserDetectorConstruction* vdet = 
                                       mgr->GetUserDetectorConstruction();
    B02ImportanceDetectorConstruction* bdet = 
                                      (B02ImportanceDetectorConstruction*)vdet;
    //
    
    //---------------------------------------------
    // Dump accumulated quantities for this RUN.
    //  (Display only central region of x-y plane)
    //      0       ConcreteSD/Collisions
    //      1       ConcreteSD/CollWeight
    //      2       ConcreteSD/Population
    //      3       ConcreteSD/TrackEnter
    //      4       ConcreteSD/SL
    //      5       ConcreteSD/SLW
    //      6       ConcreteSD/SLWE
    //      7       ConcreteSD/SLW_V
    //      8       ConcreteSD/SLWE_V
    //---------------------------------------------
    G4THitsMap<G4double>* Collisions = 
                          b02Run->GetHitsMap(fSDName[i]+"/Collisions");
    G4THitsMap<G4double>* CollWeight = 
                          b02Run->GetHitsMap(fSDName[i]+"/CollWeight");
    G4THitsMap<G4double>* Population = 
                          b02Run->GetHitsMap(fSDName[i]+"/Population");
    G4THitsMap<G4double>* TrackEnter = 
                          b02Run->GetHitsMap(fSDName[i]+"/TrackEnter");
    G4THitsMap<G4double>* SL = b02Run->GetHitsMap(fSDName[i]+"/SL");
    G4THitsMap<G4double>* SLW = b02Run->GetHitsMap(fSDName[i]+"/SLW");
    G4THitsMap<G4double>* SLWE = b02Run->GetHitsMap(fSDName[i]+"/SLWE");
    G4THitsMap<G4double>* SLW_V = b02Run->GetHitsMap(fSDName[i]+"/SLW_V");
    G4THitsMap<G4double>* SLWE_V = b02Run->GetHitsMap(fSDName[i]+"/SLWE_V");

    if (IsMaster())
      {
        G4cout << 
          "\n--------------------End of Global Run-----------------------" << 
        G4endl;
        G4cout << 
          " Number of event processed : "<< aRun->GetNumberOfEvent() << G4endl;
      }
    else
      {
        G4cout << 
          "\n--------------------End of Local Run------------------------" << 
        G4endl;
        G4cout << 
          " Number of event processed : "<< aRun->GetNumberOfEvent() << G4endl;
      }      
    
    G4cout << "=============================================================" 
           <<G4endl;
    G4cout << "=============================================================" 
           <<G4endl;

    std::ostream *myout = &G4cout;
    PrintHeader(myout);

    for ( G4int iz = 0; iz < 20; iz++){ 
      G4double* SumCollisions = (*Collisions)[iz];
      G4double* SumCollWeight = (*CollWeight)[iz];
      G4double* Populations   = (*Population)[iz];
      G4double* TrackEnters   = (*TrackEnter)[iz];
      G4double* SLs   = (*SL)[iz];
      G4double* SLWs   = (*SLW)[iz];
      G4double* SLWEs   = (*SLWE)[iz];
      G4double* SLW_Vs   = (*SLW_V)[iz];
      G4double* SLWE_Vs   = (*SLWE_V)[iz];
      if ( !SumCollisions ) SumCollisions = new G4double(0.0);
      if ( !SumCollWeight ) SumCollWeight = new G4double(0.0);
      if ( !Populations   ) Populations   = new G4double(0.0);
      if ( !TrackEnters   ) TrackEnters   = new G4double(0.0);
      if ( !SLs   ) SLs   = new G4double(0.0);
      if ( !SLWs   ) SLWs   = new G4double(0.0);
      if ( !SLWEs   ) SLWEs   = new G4double(0.0);
      if ( !SLW_Vs   ) SLW_Vs   = new G4double(0.0);
      if ( !SLWE_Vs   ) SLWE_Vs   = new G4double(0.0);
      G4double NumWeightedEnergy =0.0;
      G4double FluxWeightedEnergy=0.0; 
      G4double AverageTrackWeight=0.0;
      if ( *SLW_Vs !=0. ) NumWeightedEnergy  = (*SLWE_Vs)/(*SLW_Vs);
      if ( *SLWs   !=0. ) FluxWeightedEnergy  = (*SLWEs)/(*SLWs);
      if ( *SLs    !=0. ) AverageTrackWeight  = (*SLWs)/(*SLs);
      G4String cname = bdet->GetCellName(iz);
      G4cout 
        << std::setw(fFieldValue) << cname << " |"
        << std::setw(fFieldValue) << (*TrackEnters) << " |"
        << std::setw(fFieldValue) << (*Populations) << " |"
        << std::setw(fFieldValue) << (*SumCollisions) << " |"
        << std::setw(fFieldValue) << (*SumCollWeight) << " |"
        << std::setw(fFieldValue) << NumWeightedEnergy << " |"
        << std::setw(fFieldValue) << FluxWeightedEnergy << " |"
        << std::setw(fFieldValue) << AverageTrackWeight << " |"
        << std::setw(fFieldValue) << (*SLs) << " |"
        << std::setw(fFieldValue) << (*SLWs) << " |"
        << std::setw(fFieldValue) << (*SLW_Vs) << " |"
        << std::setw(fFieldValue) << (*SLWEs) << " |"
        << std::setw(fFieldValue) << (*SLWE_Vs) << " |"
        << G4endl;
    }
    G4cout << "============================================="<<G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02RunAction::PrintHeader(std::ostream *out)
{
  std::vector<G4String> vecScoreName;
  vecScoreName.push_back("Tr.Entering");
  vecScoreName.push_back("Population");
  vecScoreName.push_back("Collisions");
  vecScoreName.push_back("Coll*WGT");
  vecScoreName.push_back("NumWGTedE");
  vecScoreName.push_back("FluxWGTedE");
  vecScoreName.push_back("Av.Tr.WGT");
  vecScoreName.push_back("SL");
  vecScoreName.push_back("SLW");
  vecScoreName.push_back("SLW_v");
  vecScoreName.push_back("SLWE");
  vecScoreName.push_back("SLWE_v");

  // head line
//   std::string vname;
//   vname = FillString("Volume", ' ', fFieldName+1);
  //*out << vname << '|';
  *out << std::setw(fFieldValue) << "Volume" << " |";
  for (std::vector<G4String>::iterator it = vecScoreName.begin();
      it != vecScoreName.end(); it++) {
//      vname = FillString((*it),
//                         ' ', 
//                         fFieldValue+1, 
//                         false);
//      *out << vname << '|';
      *out << std::setw(fFieldValue) << (*it) << " |";
  }
  *out << G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::string B02RunAction::FillString(const std::string &name, 
                                       char c, G4int n, G4bool back)
{
  std::string fname("");
  G4int k = n - name.size();
  if (k > 0) {
    if (back) {
      fname = name;
      fname += std::string(k,c);
    }
    else {
      fname = std::string(k,c);
      fname += name;
    }
  }
  else {
    fname = name;
  }
  return fname;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
