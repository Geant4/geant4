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
//   Name of file:       FCALHadModule.cc
//   Author:             Mathieu Fontaine           Rachid Mazini
//                       fontainerlps.umontreal.ca  Rachid.Mazinircern.ch
//   Language:           C++
//   Tested on :         g++
//   Prerequisites:      None
//   Purpose:            Source file defining the geometry of HadModule 0 of the
//                       FCAL.
//   Developped:         10-March-2000   M.F.
//
//-----------------------------------------------------------------------------

#include <fstream>

#include "FCALHadModule.hh"

#include "FCALMaterialConsultant.hh"
#include "FCALHadModuleSD.hh"
#include "G4SDManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"


FCALHadModule::FCALHadModule() :
  FcalHadModuleSD(0)
{
  F2LArGapID = new G4int[2600];
  F2LArIX     = new G4int[2600];
  F2LArJY     = new G4int[2600];
  F2LArITile  = new G4int[2600];
  F2LArGapPosX = new G4double[2600];
  F2LArGapPosY = new G4double[2600];
}

FCALHadModule::~FCALHadModule() {
  delete [] F2LArGapID;
  delete [] F2LArGapPosX;
  delete [] F2LArGapPosY;
  delete [] F2LArIX;
  delete [] F2LArJY;
  delete [] F2LArITile;
}


void FCALHadModule::InitializeGeometry() {

#include "FCALHadModuleParameters.input"

  std::ifstream File
   ("geom_data/FCal2Electrodes.dat");
   
   if(!File)  G4cerr << "Failed to open file FCal2Electrode data file" << G4endl;
   File.seekg(0);

   NF2LarGap = 0;
   while(!(File.eof())) {
     NF2LarGap++;
     File >> F2LArGapID[NF2LarGap] >> F2LArGapPosX[NF2LarGap] >> F2LArGapPosY[NF2LarGap]
	  >> F2LArIX[NF2LarGap] >>  F2LArJY[NF2LarGap] >> F2LArITile[NF2LarGap];
   };
   
   G4cout << "*********" << " Number of Rods in FCAL2 : " << NF2LarGap-1 << G4endl;
}


G4LogicalVolume * FCALHadModule::Construct()
{
  //-----------------------------
  // construction of materials
  //-----------------------------
  
  FCALMaterialConsultant * FCALMaterials = 
    FCALMaterialConsultant::GetInstance();
 
  G4VisAttributes * ColorOfTungsten = new G4VisAttributes(G4Colour(.5,.5,.5));
  G4VisAttributes * ColorOfCopper =new G4VisAttributes(G4Colour(0.58,0.15,0.05));
  G4VisAttributes * ColorOfLarg = new  G4VisAttributes(G4Colour(0.,0.,1.));


  //----------------------------
  //   Read Parameters
  //----------------------------
  InitializeGeometry();


  //-----------------------------------------
  // the logical to be returned (mother)
  //-----------------------------------------

  G4Tubs * SolidHadModule =
    new G4Tubs("HadModuleSolid", HadModuleRMin, HadModuleRMax, HadModuleLenght,
	       HadModuleStartPhi,HadModuleDPhi);
  G4LogicalVolume * LogicalHadModule = 
    new G4LogicalVolume(SolidHadModule, FCALMaterials->Material("Copper"),
			"HadModuleLogical");
 
  LogicalHadModule->SetSmartless(FCAL2HadSmart);
  
  LogicalHadModule->SetVisAttributes(ColorOfCopper);
  //  LogicalHadModule->SetVisAttributes(G4VisAttributes::GetInvisible());


  //-----------------------------------------
  //  Tungsten Absorber
  //-----------------------------------------
  G4Tubs * SolidWAbsorber = 
    new G4Tubs("WAbsorberSolid", WAbsorberRMin, WAbsorberRMax, WAbsorberLenght,
	       WAbsorberStartPhi, WAbsorberDPhi);      
  G4LogicalVolume * LogicalWAbsorber = 
    new G4LogicalVolume(SolidWAbsorber, FCALMaterials->Material("FCAL2WFeNi"),
			"SolidWLogical");
//  G4VPhysicalVolume * PhysicalWAbsorber =
    new G4PVPlacement(0, G4ThreeVector(), LogicalWAbsorber, "WAbsorberPhysical",
		      LogicalHadModule, 0, 0);

  LogicalWAbsorber->SetVisAttributes(ColorOfTungsten);
  // LogicalWAbsorber->SetVisAttributes(G4VisAttributes::GetInvisible());


  // -----------------
  //  Copper Plates
  //------------------
  G4Tubs * SolidCuPlate = 
    new G4Tubs("CuPlateSolid",HadModuleRMin, HadModuleRMax, CuPlateLenght, 
	       HadModuleStartPhi, HadModuleDPhi);
  G4LogicalVolume * LogicalCuPlate =
    new G4LogicalVolume(SolidCuPlate, FCALMaterials->Material("Copper"), "CuPlateLogical");

//  G4VPhysicalVolume * PhysicalCuPlateA =
    new G4PVPlacement(0, G4ThreeVector(0.,0.,CuPlateAPosZ), LogicalCuPlate, 
		      "CuPlateAPhysical", LogicalHadModule, 0, 0);
//  G4VPhysicalVolume * PhysicalCuPlateB =
    new G4PVPlacement(0, G4ThreeVector(0.,0.,CuPlateBPosZ), LogicalCuPlate, 
		      "CuPlateBPhysical", LogicalHadModule, 0, 0);

  LogicalCuPlate->SetVisAttributes(ColorOfCopper);
  //  LogicalCuPlate->SetVisAttributes(G4VisAttributes::GetInvisible());

  //------------------------------------------
  // Had Module (F2)  Main and A/B Cable Troff 
  //------------------------------------------
  G4Tubs * SolidF2TroffMain = 
    new G4Tubs("F2TroffMainSolid", F2TroffRmin, F2TroffRmax, F2TroffMainLenght, 
	       F2TroffStartPhi, F2TroffDphi);
  G4LogicalVolume * LogicalF2TroffMain =
    new G4LogicalVolume(SolidF2TroffMain, FCALMaterials->Material("FCAL2CuArKap"),
			"F2TroffMainLogical");
  
  G4Tubs * SolidF2TroffAB = 
    new G4Tubs("F2TroffABSolid", F2TroffRmin, F2TroffRmax, F2TroffABLenght, 
	       F2TroffStartPhi, F2TroffDphi);
  G4LogicalVolume * LogicalF2TroffAB =
    new G4LogicalVolume(SolidF2TroffAB, FCALMaterials->Material("FCAL2CuArKap"),
			"F2TroffABLogical");
  
  G4ThreeVector F2TroffMainTrans(0.,0.,0.);
  G4ThreeVector F2TroffABTrans(0.,0.,0.);
  G4RotationMatrix F2TroffRot;
  G4int i=0;
    for(i=0 ; i < NCableTroff ; i++)
      {
//      G4VPhysicalVolume * PhysicalF2TroffMain =
	new G4PVPlacement(G4Transform3D(F2TroffRot,F2TroffMainTrans), LogicalF2TroffMain,
			  "F2TroffMainPhysical", LogicalWAbsorber,0,i+1);
      
//      G4VPhysicalVolume * PhysicalF2TroffAB = 
	new G4PVPlacement(G4Transform3D(F2TroffRot,F2TroffABTrans), LogicalF2TroffAB, 
			  "F2TroffAPhysical", LogicalCuPlate, 0, i+1);
      
      F2TroffRot.rotateZ(F2TroffRotZ);
    }

  LogicalF2TroffMain->SetVisAttributes(ColorOfCopper);
  //  LogicalF2TroffMain->SetVisAttributes(G4VisAttributes::GetInvisible());
  LogicalF2TroffAB->SetVisAttributes(ColorOfCopper);
  // LogicalF2TroffAB->SetVisAttributes(G4VisAttributes::GetInvisible());


   //----------------------
   //  LArg Gaps  + F2 Rod
   //----------------------
   G4Tubs * SolidF2LArGap = 
     new G4Tubs("F2LArGapSolid", F2LArGapRmin, F2LArGapRmax, F2LArGapLenght, 
                 F2LArGapStartPhi, F2LArGapDphi);
   G4LogicalVolume * LogicalF2LArGap = 
      new G4LogicalVolume(SolidF2LArGap, FCALMaterials->Material("LiquidArgon"),
                          "F2LArGapLogical");

     LogicalF2LArGap->SetVisAttributes(ColorOfLarg);
   // LogicalF2LArGap->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4Tubs * SolidF2Rod =
      new G4Tubs("F2RodSolid", F2RodRmin, F2RodRmax, F2RodLenght, F2RodStartPhi, F2RodDphi);
   G4LogicalVolume * LogicalF2Rod = 
      new G4LogicalVolume(SolidF2Rod, FCALMaterials->Material("Tungsten"),"F2RodLogical");
//    G4VPhysicalVolume * PhysicalF2Rod = 
      new G4PVPlacement(0,G4ThreeVector(),LogicalF2Rod,"F2RodPhysical",LogicalF2LArGap,0, 0);

    LogicalF2Rod->SetVisAttributes(ColorOfTungsten);
    // LogicalF2Rod->SetVisAttributes(G4VisAttributes::GetInvisible());

    //---------------------------------
    // Electrod (Rod + LArg) placement
    //---------------------------------
    for(i=1; i < NF2LarGap; i++){ 
//      G4VPhysicalVolume * PhysicalF2LArGap =
	new G4PVPlacement(0,G4ThreeVector(F2LArGapPosX[i]*cm,F2LArGapPosY[i]*cm,0.*cm),
			  LogicalF2LArGap,"F2LArGapPhysical",
			  LogicalHadModule, 0, i); 
    };

    LogicalF2LArGap->SetVisAttributes(ColorOfLarg);
    // LogicalF2LArGap->SetVisAttributes(G4VisAttributes::GetInvisible());


    // Sensitive Volumes
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
    if(!FcalHadModuleSD)
      {
	FcalHadModuleSD = new FCALHadModuleSD("FCALTB/HadModuleSD");
	SDman->AddNewDetector(FcalHadModuleSD);
      }
    LogicalF2LArGap->SetSensitiveDetector(FcalHadModuleSD);


   return LogicalHadModule;

}

G4int FCALHadModule::GetF2TileID(G4int TileID) 
{
  return F2LArITile[TileID];
}



