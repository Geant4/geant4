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
// $Id: Tst33RunAction.cc,v 1.3 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
#include "Tst33RunAction.hh"
#include "Tst33Run.hh"

//-- In order to obtain detector information.
#include "G4RunManager.hh"
//
//#include "Tst33DetectorConstruction.hh"
//test
#include "Tst33SlobedConcreteShield.hh"
#include "G4THitsMap.hh"

#include "G4UnitsTable.hh"

#include "G4GeometryCell.hh"

//=======================================================================
// Tst33RunAction
//  
//
//
//=======================================================================
// Constructor
Tst33RunAction::Tst33RunAction():
  FieldName(13),
  FieldValue(13)
//   FieldName(25),
//   FieldValue(14)
{
  // - Prepare data member for Tst33Run.
  //   vector represents a list of MultiFunctionalDetector names.

  theSDName.push_back(G4String("PhantomSD"));
//double - i.e. same logical volume as PhantomSD   theSDName.push_back(G4String("PhantomSDColl"));
  theSDName.push_back(G4String("PhantomSDMinus"));
  theSDName.push_back(G4String("PhantomSDPlus"));
}

// Destructor.
Tst33RunAction::~Tst33RunAction()
{
  theSDName.clear();
}

//
//== 
G4Run* Tst33RunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  //  Detail description can be found in Tst33Run.hh/cc.
  return new Tst33Run(theSDName);
}

//
//==
void Tst33RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

//
//== 
void Tst33RunAction::EndOfRunAction(const G4Run* aRun)
{
    G4cout << " ###### EndOfRunAction  " <<G4endl;
  //- Tst33Run object.
  Tst33Run* re02Run = (Tst33Run*)aRun;
  //--- Dump all socred quantities involved in Tst33Run.
  //  re02Run->DumpAllScorer();
  //---
  G4RunManager* mgr = G4RunManager::GetRunManager();
  //
  
  for ( G4int i = 0; i < (G4int)theSDName.size(); i++ ){
    const G4VUserDetectorConstruction* vdet = mgr->GetUserDetectorConstruction();
    //test    Tst33DetectorConstruction* bdet = (Tst33DetectorConstruction*)vdet;
    Tst33SlobedConcreteShield* bdet = (Tst33SlobedConcreteShield*)vdet;
    //
    
    //---------------------------------------------
    // Dump accumulated quantities for this RUN.
    //  (Display only central region of x-y plane)
    //      0       PhantomSD/Collisions
    //      1       PhantomSD/CollWeight
    //      2       PhantomSD/Population
    //      3       PhantomSD/TrackEnter
    //      4       PhantomSD/SL
    //      5       PhantomSD/SLW
    //      6       PhantomSD/SLWE
    //      7       PhantomSD/SLW_V
    //      8       PhantomSD/SLWE_V
    //---------------------------------------------
    G4THitsMap<G4double>* Collisions = re02Run->GetHitsMap(theSDName[i]+"/Collisions");
    G4THitsMap<G4double>* CollWeight = re02Run->GetHitsMap(theSDName[i]+"/CollWeight");
    //x    G4THitsMap<G4double>* Population = re02Run->GetHitsMap(theSDName[i]+"/Population");
    G4THitsMap<G4double>* TrackEnter = re02Run->GetHitsMap(theSDName[i]+"/TrackEnter");
    G4THitsMap<G4double>* SL = re02Run->GetHitsMap(theSDName[i]+"/SL");
    G4THitsMap<G4double>* SLW = re02Run->GetHitsMap(theSDName[i]+"/SLW");
    G4THitsMap<G4double>* SLWE = re02Run->GetHitsMap(theSDName[i]+"/SLWE");
    G4THitsMap<G4double>* SLW_V = re02Run->GetHitsMap(theSDName[i]+"/SLW_V");
    G4THitsMap<G4double>* SLWE_V = re02Run->GetHitsMap(theSDName[i]+"/SLWE_V");


    G4cout << "=============================================================" <<G4endl;
    G4cout << " Number of event processed : "<< aRun->GetNumberOfEvent() << " for hit collection: " << theSDName[i] << G4endl;
    G4cout << "=============================================================" <<G4endl;

    std::ostream *myout = &G4cout;
    PrintHeader(myout);

    for ( G4int iz = 0; iz < 20; iz++){ 
      G4double* SumCollisions = (*Collisions)[iz];
      G4double* SumCollWeight = (*CollWeight)[iz];
      //x      G4double* Populations   = (*Population)[iz];
      G4double* TrackEnters   = (*TrackEnter)[iz];
      G4double* SLs   = (*SL)[iz];
      G4double* SLWs   = (*SLW)[iz];
      G4double* SLWEs   = (*SLWE)[iz];
      G4double* SLW_Vs   = (*SLW_V)[iz];
      G4double* SLWE_Vs   = (*SLWE_V)[iz];
      if ( !SumCollisions ) SumCollisions = new G4double(0.0);
      if ( !SumCollWeight ) SumCollWeight = new G4double(0.0);
      //x      if ( !Populations   ) Populations   = new G4double(0.0);
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
//       G4GeometryCell gCell(bdet->GetGeometryCell(i, ""));
//       G4String cname = gCell.GetPhysicalVolume().GetName();
//test      G4String cname = "bob";
      //.GetPhysicalVolume()->GetName();
      G4String cname = bdet->GetCellName(iz);
      G4cout 
	<< std::setw(FieldValue) << cname << " |"
	<< std::setw(FieldValue) << (*TrackEnters) << " |"
	//x	<< std::setw(FieldValue) << (*Populations) << " |"
	<< std::setw(FieldValue) << (*SumCollisions) << " |"
	<< std::setw(FieldValue) << (*SumCollWeight) << " |"
	<< std::setw(FieldValue) << NumWeightedEnergy << " |"
	<< std::setw(FieldValue) << FluxWeightedEnergy << " |"
	<< std::setw(FieldValue) << AverageTrackWeight << " |"
	<< std::setw(FieldValue) << (*SLs) << " |"
	<< std::setw(FieldValue) << (*SLWs) << " |"
	<< std::setw(FieldValue) << (*SLW_Vs) << " |"
	<< std::setw(FieldValue) << (*SLWEs) << " |"
	<< std::setw(FieldValue) << (*SLWE_Vs) << " |"
	<< G4endl;
    }
    G4cout << "============================================="<<G4endl;
  }
}
//
// --

void Tst33RunAction::PrintHeader(std::ostream *out)
{
  std::vector<G4String> vecScoreName;
  vecScoreName.push_back("Tr.Entering");
  //x  vecScoreName.push_back("Population");
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
  //std::string vname = FillString("Volume", ' ', FieldName+1);
  //*out << vname << '|';
  std::string vname;
  *out << std::setw(FieldValue) << "Volume" << " |";
  for (std::vector<G4String>::iterator it = vecScoreName.begin();
       it != vecScoreName.end(); it++) {
      //vname = FillString((*it),
//		       ' ', 
//		       FieldValue+1, 
//		       false);
//    *out << vname << '|';
      *out << std::setw(FieldValue) << (*it) << " |";
  }
  *out << G4endl;  
}

std::string Tst33RunAction::FillString(const std::string &name, 
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
