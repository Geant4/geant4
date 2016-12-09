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
//   Author:             Rachid Mazini
//                       Rachid.Mazini@cern.ch
//   Language:           C++
//   Tested on :         g++ (egcs.2.1.1, RH6.1)
//   Prerequisites:      None
//   Purpose:            Source file defining the geometry of EMModule 0 of the
//                       FCAL.
//   Developped:         10-March-2000   R.M.
//   
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <cstdlib>

#include "FCALEMModule.hh"

#include "FCALMaterialConsultant.hh"

#include "G4SDManager.hh"
#include "FCALEMModuleSD.hh"

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

FCALEMModule::FCALEMModule() : 
  FcalEmModuleSD(0)
 {
  F1LArGapID   = new G4int[2400];
  F1LArIX      = new G4int[2400];
  F1LArJY      = new G4int[2400];
  F1LArITile   = new G4int[2400];
  F1LArGapPosX = new G4double[2400];
  F1LArGapPosY = new G4double[2400];
}


FCALEMModule::~FCALEMModule(){
  delete [] F1LArGapID;
  delete [] F1LArGapPosX;
  delete [] F1LArGapPosY;
  delete [] F1LArIX;
  delete [] F1LArJY;
  delete [] F1LArITile;
}


void FCALEMModule::InitializeGeometry() {
#include "FCALEMModuleParameters.input"  
  std::ifstream File
    ("geom_data/FCal1Electrodes.dat");
  
  if(!File)  G4cerr << "Failed to open file FCal1Electrodes data file  " << G4endl;
  File.seekg(0);
  
  NF1LarGap = 0;
  while(!(File.eof())) {
    NF1LarGap++;
    File >> F1LArGapID[NF1LarGap] >> F1LArGapPosX[NF1LarGap] >> F1LArGapPosY[NF1LarGap]
	 >> F1LArIX[NF1LarGap] >>  F1LArJY[NF1LarGap] >> F1LArITile[NF1LarGap];
  };   
  G4cout << "********" << " Number of Rods in FCAL1 : " << NF1LarGap-1 << G4endl;;
}



G4LogicalVolume * FCALEMModule::Construct()
{
  //-----------------------------
  // construction of materials
  //----------------------------- 
  FCALMaterialConsultant *FCALMaterials = 
    FCALMaterialConsultant::GetInstance();

  G4VisAttributes * ColorOfEMModule = new G4VisAttributes(G4Colour(1.,0.,0.5));
//  G4VisAttributes * ColorOfLArg = new G4VisAttributes(G4Colour(0.,0.,1.));

  //----------------------------
  //      Read Parameters
  //----------------------------
  InitializeGeometry();

//-----------------------------------------
// Logical to be returned (FCAL EM module)
//-----------------------------------------
  G4Tubs * SolidEmModule =
    new G4Tubs("EmModuleSold", EmModuleRMin, EmModuleRMax, EmModuleLenght,
	       EmModuleStartPhi,EmModuleDPhi);
  G4LogicalVolume * LogicalEmModule = 
    new G4LogicalVolume(SolidEmModule, FCALMaterials->Material("Copper"),
			"EmModuleLogical");
 
  LogicalEmModule->SetSmartless(FCALEmSmart);

   LogicalEmModule->SetVisAttributes(ColorOfEMModule);
  //  LogicalEmModule->SetVisAttributes(G4VisAttributes::GetInvisible());


//---------------------
//  FCAL Cable Troff  
//---------------------
  G4Tubs * SolidF1CableTroff =
    new G4Tubs("F1CableTroffSolid", F1CableTroffRMin, F1CableTroffRMax,
	       F1CableTroffLenght, F1CableTroffStartPhi, F1CableTroffDPhi);
  G4LogicalVolume * LogicalF1CableTroff =
    new G4LogicalVolume(SolidF1CableTroff, FCALMaterials->Material("FCAL1CuArKap"),
			"F1CableTroffLogical");

  G4ThreeVector F1CableTroffTrans(0.,0.,0.);
  G4RotationMatrix F1CableTroffRot;

  G4int i=0;
  for(i=0 ; i < NCableTroff ; i++)
    {
//      G4VPhysicalVolume * PhysicalF1CableTroff =
	new G4PVPlacement(G4Transform3D(F1CableTroffRot,F1CableTroffTrans),
			  LogicalF1CableTroff,"F1CableTroffPhysical",
			  LogicalEmModule,0,i+1);

      F1CableTroffRot.rotateZ(F1CableTroffRotZ);
    }

  LogicalF1CableTroff->SetVisAttributes(ColorOfEMModule);
  // LogicalF1CableTroff->SetVisAttributes(G4VisAttributes::GetInvisible());


   //----------------------
   //    LArg gaps
   //----------------------

  G4Tubs * SolidF1LArGap = 
    new G4Tubs("F1LArGapSolid",F1LArGapRmin, F1LArGapRmax, F1LArGapLenght, 
	       F1LArGapStartPhi,F1LArGapDPhi);
	
   G4LogicalVolume * LogicalF1LArGap = 
     new G4LogicalVolume(SolidF1LArGap, FCALMaterials->Material("LiquidArgon"),
			 "LArg Gap");
  
   for(i=1; i < NF1LarGap; i++){
//     G4VPhysicalVolume * PhysicalF1LArGap =
       new G4PVPlacement(0,G4ThreeVector(F1LArGapPosX[i]*cm,F1LArGapPosY[i]*cm,0.*cm),
			  LogicalF1LArGap,"F1LArGapPhysical", LogicalEmModule, 0, i); 
   };

   // LogicalF1LArGap->SetVisAttributes(ColorOfLArg);
   LogicalF1LArGap->SetVisAttributes(G4VisAttributes::GetInvisible());
   

    // Sensitive Volumes
   G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
   if(!FcalEmModuleSD)
     {
       FcalEmModuleSD = new FCALEMModuleSD("FCALTB/EmModuleSD");
       SDman->AddNewDetector(FcalEmModuleSD);
     }
   LogicalF1LArGap->SetSensitiveDetector(FcalEmModuleSD);

   

   return LogicalEmModule;

}


G4int FCALEMModule::GetF1TileID(G4int GapID) 
{ return F1LArITile[GapID]; }

G4double FCALEMModule::GetF1LArGapPosX(G4int GapID)
{ return F1LArGapPosX[GapID]; }


