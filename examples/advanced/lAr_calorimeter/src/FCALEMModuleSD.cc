// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALEMModuleSD.cc,v 1.2 2002-10-01 13:53:17 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FCALEMModuleSD.hh"

#include "ExN03CalorHit.hh"

#include "FCALTestbeamSetup.hh"
#include "FCALEMModule.hh"
#include "ExN03SteppingAction.hh"

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

FCALEMModuleSD::FCALEMModuleSD(G4String name) : G4VSensitiveDetector(name)
{
  EmModule = new FCALEMModule();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALEMModuleSD::~FCALEMModuleSD()
{
  delete EmModule;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALEMModuleSD::Initialize(G4HCofThisEvent*HCE)
{
  if (Init_state==0) 
    {
      EmModule->InitializeGeometry();
      Init_state++;
    };

  for (G4int j=0 ; j<1131 ; j++) {EvisF1Tile[j] = 0.;}    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FCALEMModuleSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep==0.) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());  
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();


  if(strcmp(physVol->GetName(),"F1LArGapPhysical")==0){
    G4int F1LArGapId = physVol->GetCopyNo();
    G4int F1TileId = EmModule->GetF1TileID(F1LArGapId);
    EvisF1Tile[F1TileId] = EvisF1Tile[F1TileId] + edep;
  };

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALEMModuleSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  G4int NF1Tile = 0;
  G4int AddTileP[300];
  G4double EvisTileP[300];

  for (G4int i=1; i <= 1130; i++){
    if(EvisF1Tile[i] > 0.) {
      NF1Tile++;
      AddTileP[NF1Tile] = i;
      EvisTileP[NF1Tile] = EvisF1Tile[i];
    };};

  G4cout << "Number of F1 Tiles with Positive energy : " << NF1Tile <<  endl;

  // Write in File
  //--------------
  G4String FileName = "EmModule_802_1mm.dat";
  G4int iostemp;
  if(Init_state == 1) {
    iostemp = ios::out;
    Init_state++;
  } else {
    iostemp = ios::out|ios::app; // ios::app;  
  };
  
  ofstream EmDatafile(FileName, iostemp);
  // EmDatafile.precision(5);

  EmDatafile << NF1Tile << endl;
  for (i=1; i <= NF1Tile; i++) {
    EmDatafile << AddTileP[i] << " " << EvisTileP[i]/MeV << endl;
  }
  EmDatafile.close();

}
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALEMModuleSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALEMModuleSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALEMModuleSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

