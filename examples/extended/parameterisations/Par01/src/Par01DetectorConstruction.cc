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
/// \file Par01/src/Par01DetectorConstruction.cc
/// \brief Implementation of the Par01DetectorConstruction class
//
//
//
#include "Par01DetectorConstruction.hh"
#include "Par01CalorimeterSD.hh"
#include "Par01EMShowerModel.hh"
#include "Par01PiModel.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ProductionCuts.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01DetectorConstruction::Par01DetectorConstruction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01DetectorConstruction::~Par01DetectorConstruction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Par01DetectorConstruction::Construct()
{
  G4cout << "\nPar01DetectorConstruction....\n" << G4endl;
  
  //--------- Material definition ---------
  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();

  // Build materials
  G4Material* air    = nistManager->FindOrBuildMaterial("G4_AIR");
  G4Material* csi    = nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  G4Material* helium = nistManager->FindOrBuildMaterial("G4_He");
  G4Material* iron   = nistManager->FindOrBuildMaterial("G4_Fe");

  
  //--------- G4VSolid, G4LogicalVolume, G4VPhysicalVolume  ---------
  
  //--------------
  // World:
  //--------------
  G4Box *WorldBox= new G4Box("WorldBox",400*cm, 400*cm, 400*cm);
  G4LogicalVolume *WorldLog=new G4LogicalVolume(WorldBox,air,
                                                "WorldLogical", 0, 0, 0);
  G4PVPlacement *WorldPhys=new G4PVPlacement(0,G4ThreeVector(),
                                             "WorldPhysical",
                                             WorldLog,
                                             0,false,0);
  // Size of detectors:
  G4double detectSize = 125*cm;
  
  //-----------------------------
  // "Drift Chamber":
  // Not used in parameterisation.
  //-----------------------------
  // -- Logical volume:
  G4Box *driftChamberBox
    = new G4Box("DriftChamberSolid", detectSize, detectSize, 40*cm);
  G4LogicalVolume *driftChamberLog
    = new G4LogicalVolume(driftChamberBox,helium,
                          "DriftChamberLogical", 0, 0, 0);
  // -- Placement:
  // G4PVPlacement *driftChamberPhys  =
      new G4PVPlacement(0,G4ThreeVector(0., 0., 50*cm),
                        "DriftChamberPhysical",
                        driftChamberLog,
                        WorldPhys,false,0);
  
  //--------------------------
  // "Calorimeter": used in
  // parameterisation below
  //--------------------------
  // -- Logical volume:
  G4Box *calorimeterBox
    = new G4Box("CalorimeterSolid", detectSize, detectSize, 20*cm);
  G4LogicalVolume *calorimeterLog = new G4LogicalVolume(calorimeterBox,air,
                                                        "CalorimeterLogical", 0, 0, 0);
  // -- Placement:
  G4PVPlacement *calorimeterPhys  = new G4PVPlacement(0,G4ThreeVector(0., 0., 120*cm),
                                                      "CalorimeterPhysical",
                                                      calorimeterLog,
                                                      WorldPhys,false,0);
  
  //--------------------------------------
  // The calorimeter is filled with
  // crystals:
  //--------------------------------------
  // -- Logical volume:
  G4double CrystalX = 2.5*cm;
  G4double CrystalY = CrystalX;
  G4double CrystalZ = 20*cm;
  G4Box *CrystalSolid = new G4Box("CrystalSolid", CrystalX, CrystalY, CrystalZ);
  fCrystalLog       = new G4LogicalVolume(CrystalSolid,csi,
                                            "CrystalLogical", 0, 0, 0);
  
  G4String tName1("Crystal");        // Allow all target physicals to share
  // same name (delayed copy)
  
  // -- and placements inside the calorimeter:
  G4int copyNo=0;
  G4double xTlate, yTlate;
  fnX = 48;
  fnY = 48;
  for (G4int j = 0; j < fnY; j++)
    {
      yTlate = -detectSize + 3*CrystalY + j*2*CrystalY;
      for (G4int i = 0; i < fnX; i++)
        {
          xTlate = -detectSize + 3*CrystalX + i*2*CrystalX;
          new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
                            tName1,
                            fCrystalLog,
                            calorimeterPhys,false,copyNo++);
        }
    }


  //--------------------------
  // "Hadron Calorimeter": used
  // in parameterisation with
  // a parallel geometry
  //--------------------------
  // -- Logical volume:
  G4Box *hadCaloBox
    = new G4Box("HadCaloSolid", detectSize, detectSize, 50*cm);
  G4LogicalVolume *hadCaloLog = new G4LogicalVolume(hadCaloBox,air,
                                                     "HadCaloLogical", 0, 0, 0);
  // -- Placement:
  G4PVPlacement *hadCaloPhys  = new G4PVPlacement(0,G4ThreeVector(0., 0., 200*cm),
                                                   "HadCaloPhysical",
                                                   hadCaloLog,
                                                   WorldPhys,false,0);
  
  //--------------------------------------
  // The calorimeter is filled with
  // towers:
  //--------------------------------------
  // -- Logical volume:
  G4double TowerX = 5*cm;
  G4double TowerY = TowerX;
  G4double TowerZ = 45*cm;
  G4Box *TowerSolid = new G4Box("TowerSolid", TowerX, TowerY, TowerZ);
  fTowerLog       = new G4LogicalVolume(TowerSolid,iron,
                                          "TowerLogical", 0, 0, 0);
  
  G4String tName2("Tower");
  
  // -- and placements inside the calorimeter:
  copyNo=0;
  fnXhad = 23;
  fnYhad = 23;
  for (G4int jj = 0; jj < fnYhad; jj++)
    {
      yTlate = -detectSize + 3*TowerY + jj*2*TowerY;
      for (G4int i = 0; i < fnXhad; i++)
         {
          xTlate = -detectSize + 3*TowerX + i*2*TowerX;
           new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
                             tName2,
                             fTowerLog,
                             hadCaloPhys,false,copyNo++);
         }
    }
  
  
 
  
  // -- Makes the calorimeterLog volume becoming a G4Region: 
   G4Region* caloRegion = new G4Region("EM_calo_region");
   caloRegion->AddRootLogicalVolume(calorimeterLog);
   std::vector<double> cuts; 
   cuts.push_back(1.0*mm);cuts.push_back(1.0*mm);cuts.push_back(1.0*mm);cuts.push_back(1.0*mm);
   caloRegion->SetProductionCuts(new G4ProductionCuts());
   caloRegion->GetProductionCuts()->SetProductionCuts(cuts);

   //  Makes had. calo a region to:
   G4Region* hadRegion = new G4Region("HAD_calo_region");
   hadRegion->AddRootLogicalVolume(hadCaloLog);
   cuts.clear();
   cuts.push_back(1.0*cm);cuts.push_back(1.0*cm);cuts.push_back(1.0*cm);cuts.push_back(1.0*cm);
   hadRegion->SetProductionCuts(new G4ProductionCuts());
   hadRegion->GetProductionCuts()->SetProductionCuts(cuts);

  //--------- Visualization attributes -------------------------------
  WorldLog->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VisAttributes * driftchamberTubeVisAtt
    = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  driftchamberTubeVisAtt->SetForceWireframe(true);
  driftChamberLog->SetVisAttributes(driftchamberTubeVisAtt);
  
  G4VisAttributes * calorimeterBoxVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  calorimeterBoxVisAtt->SetForceWireframe(true);
  calorimeterLog->SetVisAttributes(calorimeterBoxVisAtt);
  
  G4VisAttributes * crystalVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  crystalVisAtt->SetForceWireframe(true);
  fCrystalLog->SetVisAttributes(crystalVisAtt);

  G4VisAttributes * hadCaloBoxVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  hadCaloBoxVisAtt->SetForceWireframe(true);
  hadCaloLog->SetVisAttributes(hadCaloBoxVisAtt);
  
  G4VisAttributes * towerVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
  towerVisAtt->SetForceWireframe(true);
  fTowerLog->SetVisAttributes(towerVisAtt);
  
  //------------------------------------------------------------------
  
  
  //-----------------------
  // Returns the pointer to
  // the physical world:
  //-----------------------
  return WorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01DetectorConstruction::ConstructSDandField()
{
  //--------- Sensitive detector -------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String calorimeterSDname = "Par01/Calorimeter";
  Par01CalorimeterSD* CalorimeterSD = new Par01CalorimeterSD( calorimeterSDname,
                                                              fnX*fnY,
                                                              "CalCollection" );
  SDman->AddNewDetector( CalorimeterSD );
  fCrystalLog->SetSensitiveDetector(CalorimeterSD); 
  
  G4String hadCalorimeterSDname = "Par01/HadronCalorimeter";
  Par01CalorimeterSD* HadCalorimeterSD = new Par01CalorimeterSD( hadCalorimeterSDname,
                                                                 fnXhad*fnYhad,
                                                                 "HadCollection" );
  SDman->AddNewDetector( HadCalorimeterSD );
  fTowerLog->SetSensitiveDetector(HadCalorimeterSD); 
  
  // --------------- fast simulation ----------------------------
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  
  G4Region* caloRegion = regionStore->GetRegion("EM_calo_region");
  // builds a model and sets it to the envelope of the calorimeter:
  new Par01EMShowerModel("emShowerModel",caloRegion);
}
