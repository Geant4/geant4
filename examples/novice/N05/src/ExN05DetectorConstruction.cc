// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05DetectorConstruction.cc,v 1.2 1999-12-15 14:49:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "ExN05DetectorConstruction.hh"
#include "ExN05CalorimeterSD.hh"
#include "ExN05DetectorMessenger.hh"
#include "ExN05EMShowerModel.hh"
#include "ExN05PionShowerModel.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4FastSimulationManager.hh"


ExN05DetectorConstruction::ExN05DetectorConstruction()
{
  theUserLimitsForCrystal = NULL; 
  fUseUserLimits = false;
  theMaxTimeCutsInCrystal   = 1000.0 * ns;
  theMinEkineCutsInCrystal  = 1.0 * MeV;
  theMinRangeCutsInCrystal  = 1.0 * mm;
  theMessenger =  new ExN05DetectorMessenger(this);
}

ExN05DetectorConstruction::~ExN05DetectorConstruction()
{
  if (theUserLimitsForCrystal != NULL) delete theUserLimitsForCrystal;
  delete theMessenger;
}

G4VPhysicalVolume* ExN05DetectorConstruction::Construct()
{
  G4cout << "\nExN05DetectorConstruction....\n" << G4endl;
  
  //--------- Material definition ---------
  
  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;
  
  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N",  iz=7.,  a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen",   symbol="O",  iz=8.,  a);
  a = 126.9*g/mole;
  G4Element* elI = new G4Element(name="Iodine",   symbol="I",  iz=53., a);
  a = 132.9*g/mole;
  G4Element* elCs= new G4Element(name="Cesium",   symbol="Cs", iz=55., a);

  
  density = 4.51*g/cm3;
  G4Material* CsI = new G4Material(name="CsI", density, nel = 2);
  CsI->AddElement(elI, .5);
  CsI->AddElement(elCs,.5);
  a = 4.0*g/mole;
  density = 0.1786e-03*g/cm3;
  G4Material* He  = new G4Material(name="He", z=2., a, density);
  
  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);
  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);
  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Fer", z=26., a, density);
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);


  
  
  
  //--------- G4VSolid, G4LogicalVolume, G4VPhysicalVolume  ---------
  
  //--------------
  // World:
  //--------------
  G4Box *WorldBox= new G4Box("WorldBox",400*cm, 400*cm, 400*cm);
  G4LogicalVolume *WorldLog=new G4LogicalVolume(WorldBox,Air,
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
  G4LogicalVolume *driftChamberLog = new G4LogicalVolume(driftChamberBox,He,
							 "DriftChamberLogical", 0, 0, 0);
  // -- Placement:
  G4PVPlacement *driftChamberPhys  = new G4PVPlacement(0,G4ThreeVector(0., 0., 50*cm),
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
  G4LogicalVolume *calorimeterLog = new G4LogicalVolume(calorimeterBox,Air,
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
  theCrystalLog       = new G4LogicalVolume(CrystalSolid,CsI,
					    "CrystalLogical", 0, 0, 0);
  // create UserLimits
  if (theUserLimitsForCrystal != NULL) delete theUserLimitsForCrystal;
  theUserLimitsForCrystal = new G4UserLimits( DBL_MAX,  //step max
					      DBL_MAX,  // track max
					      theMaxTimeCutsInCrystal,
					      theMinEkineCutsInCrystal,
					      theMinRangeCutsInCrystal );
  
  // attach UserLimits   
  if (fUseUserLimits) {
    theCrystalLog->SetUserLimits(theUserLimitsForCrystal);
  }
  
  
  G4String tName1("Crystal");	// Allow all target physicals to share
  // same name (delayed copy)
  
  // -- and placements inside the calorimeter:
  G4PVPlacement *crystalPhys;
  G4int copyNo=0;
  G4double xTlate, yTlate;
  G4int nX = 48;
  G4int nY = 48;
  for (G4int j = 0; j < nY; j++)
    {
      yTlate = -detectSize + 3*CrystalY + j*2*CrystalY;
      for (G4int i = 0; i < nX; i++)
	{
	  xTlate = -detectSize + 3*CrystalX + i*2*CrystalX;
	  crystalPhys
	    = new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
				tName1,
				theCrystalLog,
				calorimeterPhys,false,copyNo++);
	}
    }

  
  //--------------------------
  // "Hadron Calorimeter": used
  // in parameterisation with
  // a ghost volume
  //--------------------------
  // -- Logical volume:
  G4Box *hadCaloBox
    = new G4Box("HadCaloSolid", detectSize, detectSize, 50*cm);
  G4LogicalVolume *hadCaloLog = new G4LogicalVolume(hadCaloBox,Air,
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
  theTowerLog       = new G4LogicalVolume(TowerSolid,Fe,
					  "TowerLogical", 0, 0, 0);
  
  G4String tName2("Tower");
  
  // -- and placements inside the calorimeter:
  G4PVPlacement *towerPhys;
  copyNo=0;
  G4int nXhad = 23;
  G4int nYhad = 23;
  for (G4int jj = 0; jj < nYhad; jj++)
    {
      yTlate = -detectSize + 3*TowerY + jj*2*TowerY;
      for (G4int i = 0; i < nXhad; i++)
 	{
	  xTlate = -detectSize + 3*TowerX + i*2*TowerX;
 	  towerPhys
 	    = new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
 				tName2,
 				theTowerLog,
 				hadCaloPhys,false,copyNo++);
 	}
    }
  
  
  //------------------------------------------------------------
  //
  //

  //--------------- Ghost Volumes ------------------------------------
  // Ghost volumes will be placed in an automatic copy of the world volume.
  // The placement of the ghost logical volume are specified by the
  // G4FastSimulationManager object attached to it.
  // Those placement are given compared to the world.
  // -- Solid:
   G4Box *ghostBox
     = new G4Box("GhostBox", detectSize+5*cm, detectSize+5*cm, 80*cm);
   // -- Logical volume:
   G4LogicalVolume* ghostLogical
     = new G4LogicalVolume(ghostBox,Air,
  			  "GhostLogical", 
  			  0, 0, 0);
  
   // G4FastSimulationManager doesn't exist yet: we set it
   // (not needed if we set a G4VFastSimulationModel which
   // takes care of creating one if needed)
   new G4FastSimulationManager(ghostLogical);
  
   ghostLogical->GetFastSimulationManager()->
     AddGhostPlacement(0,G4ThreeVector(0., 0., 175*cm));
  
  //--------- Sensitive detector -------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String calorimeterSDname = "ExN05/Calorimeter";
  ExN05CalorimeterSD* CalorimeterSD = new ExN05CalorimeterSD( calorimeterSDname, nX*nY, "CalCollection" );
  SDman->AddNewDetector( CalorimeterSD );
  theCrystalLog->SetSensitiveDetector(CalorimeterSD); 

   G4String hadCalorimeterSDname = "ExN05/HadronCalorimeter";
   ExN05CalorimeterSD* HadCalorimeterSD = new ExN05CalorimeterSD( hadCalorimeterSDname, nXhad*nYhad, "HadCollection" );
   SDman->AddNewDetector( HadCalorimeterSD );
   theTowerLog->SetSensitiveDetector(HadCalorimeterSD); 
  
  //------------------ Parameterisation ------------------------------
  // builds a model and sets it to the envelope of the calorimeter:
  ExN05EMShowerModel* emShowerModel =
    new ExN05EMShowerModel("emShowerModel",calorimeterLog);
  
  // and a model attached to the ghost:
  ExN05PionShowerModel* ghostPionShowerModel =
    new ExN05PionShowerModel("ghostPionShowerModel",ghostLogical);
  
  //--------- Visualization attributes -------------------------------
  WorldLog->SetVisAttributes(G4VisAttributes::Invisible);

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
  theCrystalLog->SetVisAttributes(crystalVisAtt);

  G4VisAttributes * hadCaloBoxVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  hadCaloBoxVisAtt->SetForceWireframe(true);
  hadCaloLog->SetVisAttributes(hadCaloBoxVisAtt);
  
  G4VisAttributes * towerVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
  towerVisAtt->SetForceWireframe(true);
  theTowerLog->SetVisAttributes(towerVisAtt);
  
  G4VisAttributes * ghostVisAtt
    = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  ghostVisAtt->SetForceWireframe(true);
  ghostLogical->SetVisAttributes(ghostVisAtt);
  
  //------------------------------------------------------------------
  
  
  //-----------------------
  // Returns the pointer to
  // the physical world:
  //-----------------------
  return WorldPhys;
}

void  ExN05DetectorConstruction::SetMaxTimeInCrystal(G4double value)
{ 
  theMaxTimeCutsInCrystal = value;  
  if (theUserLimitsForCrystal != NULL) {
    theUserLimitsForCrystal->SetUserMaxTime(value);
  }
}
void  ExN05DetectorConstruction::SetMinEkineInCrystal(G4double value)
{ 
  theMinEkineCutsInCrystal = value;  
  if (theUserLimitsForCrystal != NULL) {
    theUserLimitsForCrystal->SetUserMinEkine(value);
  }
}
void  ExN05DetectorConstruction::SetMinRangeInCrystal(G4double value)
{ 
  theMinRangeCutsInCrystal = value; 
  if (theUserLimitsForCrystal != NULL) {
    theUserLimitsForCrystal->SetUserMinRange(value);
  }
}

void  ExN05DetectorConstruction::UseUserLimits(G4bool isUse)
{
  fUseUserLimits = isUse;
  if ( fUseUserLimits && (theUserLimitsForCrystal!= NULL)) {
    theCrystalLog->SetUserLimits(theUserLimitsForCrystal);
  }    
} 
