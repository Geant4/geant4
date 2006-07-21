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
// $Id: FCALEMModuleSD.cc,v 1.12 2006-07-21 11:45:53 ribon Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FCALEMModuleSD.hh"

#include "FCALCalorHit.hh"

#include "FCALTestbeamSetup.hh"
#include "FCALEMModule.hh"
#include "FCALSteppingAction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"
#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALEMModuleSD::FCALEMModuleSD(G4String name) : G4VSensitiveDetector(name),
						Init_state(0)
{
  EmModule = new FCALEMModule();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALEMModuleSD::~FCALEMModuleSD()
{
  delete EmModule;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALEMModuleSD::Initialize(G4HCofThisEvent*)
{
  if (Init_state==0) 
    {
      EmModule->InitializeGeometry();
      Init_state++;
    };

  for (G4int j=0 ; j<1131 ; j++) {EvisF1Tile[j] = 0.;}    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FCALEMModuleSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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

void FCALEMModuleSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int NF1Tile = 0;
  G4int AddTileP[300];
  G4double EvisTileP[300];
  G4int i=0;
  for (i=1; i <= 1130; i++){
    if(EvisF1Tile[i] > 0.) {
      NF1Tile++;
      AddTileP[NF1Tile] = i;
      EvisTileP[NF1Tile] = EvisF1Tile[i];
    }
  }

  G4cout << "Number of F1 Tiles with Positive energy : " << NF1Tile <<  G4endl;

  // Write in File
  //--------------
  const char * FileName = "EmModule_802_1mm.dat";
  std::ios::openmode iostemp;
  if(Init_state == 1) {
    iostemp = std::ios::out;
    Init_state++;
  } else {
    iostemp = std::ios::out|std::ios::app; // std::ios::app;  
  }
  
  std::ofstream EmDatafile(FileName, iostemp);
  // EmDatafile.precision(5);

  EmDatafile << NF1Tile << std::endl;
  for (i=1; i <= NF1Tile; i++) {
    EmDatafile << AddTileP[i] << " " << EvisTileP[i]/MeV << std::endl;
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

