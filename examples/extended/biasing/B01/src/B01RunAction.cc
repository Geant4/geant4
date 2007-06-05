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
// $Id: B01RunAction.cc,v 1.1 2007-06-05 18:20:09 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
#include "B01RunAction.hh"
#include "B01Run.hh"

//-- In order to obtain detector information.
#include "G4RunManager.hh"
#include "B01DetectorConstruction.hh"
#include "G4THitsMap.hh"

#include "G4UnitsTable.hh"
//=======================================================================
// B01RunAction
//  
//
//
//=======================================================================
// Constructor
B01RunAction::B01RunAction():
  FieldName(25),
  FieldValue(14)
{
  // - Prepare data member for B01Run.
  //   vector represents a list of MultiFunctionalDetector names.
  theSDName.push_back(G4String("PhantomSD"));
}

// Destructor.
B01RunAction::~B01RunAction()
{
  theSDName.clear();
}

//
//== 
G4Run* B01RunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  //  Detail description can be found in B01Run.hh/cc.
  return new B01Run(theSDName);
}

//
//==
void B01RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

//
//== 
void B01RunAction::EndOfRunAction(const G4Run* aRun)
{
    G4cout << " ###### EndOfRunAction  " <<G4endl;
  //- B01Run object.
  B01Run* re02Run = (B01Run*)aRun;
  //--- Dump all socred quantities involved in B01Run.
  // re02Run->DumpAllScorer();
  //---
  G4RunManager* mgr = G4RunManager::GetRunManager();
  //
  
  for ( G4int i = 0; i < (G4int)theSDName.size(); i++ ){
    const G4VUserDetectorConstruction* vdet = mgr->GetUserDetectorConstruction();
    B01DetectorConstruction* bdet = (B01DetectorConstruction*)vdet;
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
    G4cout << " Number of event processed : "<< aRun->GetNumberOfEvent() << G4endl;
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

void B01RunAction::PrintHeader(std::ostream *out)
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

std::string B01RunAction::FillString(const std::string &name, 
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
