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
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALHadModuleSD.cc,v 1.3 2002-12-12 19:16:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FCALHadModuleSD.hh"

#include "FCALCalorHit.hh"

#include "FCALTestbeamSetup.hh"
#include "FCALHadModule.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"
#include "iostream.h"
#include "fstream.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALHadModuleSD::FCALHadModuleSD(G4String name) : G4VSensitiveDetector(name)
{
   HadModule = new FCALHadModule(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALHadModuleSD::~FCALHadModuleSD()
{
  delete HadModule;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::Initialize(G4HCofThisEvent*HCE)
{
  if(InitF2 == 0) {
    HadModule->InitializeGeometry();
    InitF2++;
  }
  for (G4int j=0; j<2330; j++) { EvisF2Tile[j]=0.;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FCALHadModuleSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep==0.) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());  
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();


    if(strcmp(physVol->GetName(),"F2LArGapPhysical")==0){
    G4int F2LArGapId = physVol->GetCopyNo();
    G4int F2TileId = HadModule->GetF2TileID(F2LArGapId);
    EvisF2Tile[F2TileId] = EvisF2Tile[F2TileId] + edep;
    };

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  G4int NF2Tile = 0;
  G4int AddTileP[300];
  G4double EvisTileP[300];

  for (G4int i=0; i<2330; i++){
    if(EvisF2Tile[i] > 0.) {
      NF2Tile++;
      AddTileP[NF2Tile] = i;
      EvisTileP[NF2Tile] = EvisF2Tile[i];
    };};

  G4cout << "Number of F2 tiles with Positive energy : " << NF2Tile <<  endl;

  // Write data in File
  //-------------------
  G4String FileName = "HadModule_802_1mm.dat";
  G4int iostemp;
  if(InitF2 == 1) {
    iostemp = ios::out;
    InitF2++;
  } else {
    iostemp = ios::out|ios::app; // ios::app;  
  };
  
  ofstream HadDatafile(FileName, iostemp);
  // EmDatafile.precision(5);

  HadDatafile << NF2Tile << endl;
  for (i=1; i <= NF2Tile; i++) {
    HadDatafile << AddTileP[i] << " " << EvisTileP[i]/MeV << endl;
  }
  HadDatafile.close();



    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

