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
//   Author:             Mathieu Fontaine           Rachid Mazini
//                       fontaine@lps.umontreal.ca  Rachid.Mazini@cern.ch
//   Language:           C++
//   Tested on :         g++
//   Prerequisites:      None
//   Purpose:            Source file defining the differents volumes
//                       in the cryostat
//   Developped:         10-March-2000   M.F.
//
//-----------------------------------------------------------------------------

#include "FCALCryostatVolumes.hh"

#include "FCALMaterialConsultant.hh"

#include "FCALEMModule.hh"
#include "FCALHadModule.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

FCALCryostatVolumes::FCALCryostatVolumes()
{
#include "FCALCryostatVolumesParameters.input"
}

FCALCryostatVolumes::~FCALCryostatVolumes() {;}

G4LogicalVolume * FCALCryostatVolumes::Construct()
{

  //-----------------------------
  // construction of materials
  //-----------------------------
  
  FCALMaterialConsultant * FCALMaterials = 
    FCALMaterialConsultant::GetInstance();


//-----------------------------------------
//  G4VisAttributes * ColorOfIron = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  G4VisAttributes * ColorOfLead = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  G4VisAttributes * ColorOfAir  = new G4VisAttributes(G4Colour(1.,1.,1.));
//  G4VisAttributes * ColorOfLarg = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));


  //-----------------------------
  // Cryostat
  //-----------------------------
  G4Tubs * SolidCryostat = 
    new G4Tubs("CryostatSolid", CryostatRMin, CryostatRMax, CryostatLenght,
	       StartingPhi, DPhi);
  G4LogicalVolume * LogicalCryostat = 
    new G4LogicalVolume(SolidCryostat,FCALMaterials->Material("Iron"),
			"CryostatLogical");
 
  // LogicalCryostat->SetVisAttributes(ColorOfIron);
  LogicalCryostat->SetVisAttributes(G4VisAttributes::GetInvisible());
 

  //------------------------------
  // Insulation
  //------------------------------
  G4Tubs * SolidInsulation = 
    new G4Tubs("InsulationSolid", InsulationRMin, InsulationRMax, 
	       InsulationLenght, StartingPhi, DPhi);
  G4LogicalVolume * LogicalInsulation = 
    new G4LogicalVolume(SolidInsulation, FCALMaterials->Material("Air"),
			"InsulationLogical");
//  G4VPhysicalVolume * PhysicalInsulation = 
    new G4PVPlacement(0, G4ThreeVector(), LogicalInsulation, "InsulationPhysical",
		      LogicalCryostat, 0, 0);
  
  LogicalInsulation->SetVisAttributes(ColorOfAir);
  // LogicalInsulation->SetVisAttributes(G4VisAttributes::GetInvisible());


  //-------------------------------------
  //  Air to replace Iron inside Cryostat
  //-------------------------------------
  /*  
  G4Tubs * SolidAirCryostat = 
    new G4Tubs("AirCryostatSolid", CryostatRMin, LArgRMax, CryostatLength, 
	       StartingPhi, DPhi);
  G4LogicalVolume * LogicalAirCryostat =
    new G4LogicalVolume(SolidAirCryostat, FCALMaterials->Material("Air"),
			"AirCryostatLogical");
  G4VPhysicalVolume * PhysicalAirCryostat =
    new G4PVPlacement(0, 0, LogicalAirCryostat, "AirCryostatPhysical",
		      LogicalCryostat, 0, 0);

   LogicalAirCryostat->SetVisAttributes(ColorOfAir);
  // LogicalAirCryostat->SetVisAttributes(G4VisAttributes::GetInvisible());	 
  */


  //--------------------
  // Liquid Argon
  //--------------------
    G4Tubs * SolidLArg = 
      new G4Tubs("LArgSolid", LArgRMin, LArgRMax, LArgLenght,StartingPhi,DPhi);
    G4LogicalVolume * LogicalLArg = 
      new G4LogicalVolume(SolidLArg, FCALMaterials->Material("LiquidArgon"),
			"LArgLogical");
    G4VPhysicalVolume * PhysicalLArg = 
      new G4PVPlacement(0,G4ThreeVector(LArgPosX, LArgPosY, LArgPosZ), 
			LogicalLArg, "LArgPhysical", LogicalCryostat, 0,0);

    // LogicalLArg->SetVisAttributes(ColorOfLarg);
    LogicalLArg->SetVisAttributes(G4VisAttributes::GetInvisible());

  //-------------------
  // Front Excluder
  //-------------------
  G4Box * SolidFrontExcluder = 
    new G4Box("FrontExcluderSolid", FrontExcluderSizeX, FrontExcluderSizeY,
	      FrontExcluderSizeZ);  
  G4LogicalVolume * LogicalFrontExcluder =
    new G4LogicalVolume(SolidFrontExcluder, FCALMaterials->Material("Air")
			, "FrontExcluderLogical");

//  G4VPhysicalVolume * PhysicalFrontExcluder =
    new G4PVPlacement(0,G4ThreeVector(FrontExcluderPosX, FrontExcluderPosY,
		      FrontExcluderPosZ), "FrontExcluderPhysical",
		      LogicalFrontExcluder, PhysicalLArg, 0,0);

  LogicalFrontExcluder->SetVisAttributes(ColorOfLead);
  // LogicalFrontExcluder->SetVisAttributes(G4VisAttributes::GetInvisible());


  //--------------------
  // Back Excluder
  //--------------------
  G4Trd * SolidBackExcluder =
    new G4Trd("BackExcluderSolid", BackExcluderSize1X, BackExcluderSize2X,
	      BackExcluderSize1Y, BackExcluderSize2Y, BackExcluderSizeZ);
  G4LogicalVolume * LogicalBackExcluder = 
    new G4LogicalVolume(SolidBackExcluder, FCALMaterials->Material("Air"),
			"BackExcluderLogical");

  G4RotationMatrix * BackExcluderRotationMatrix = new G4RotationMatrix();
  BackExcluderRotationMatrix->rotateX(BackExcluderRotX);

//  G4VPhysicalVolume * PhysicalBackExcluder =
    new G4PVPlacement(BackExcluderRotationMatrix,
		      G4ThreeVector(BackExcluderPosX, BackExcluderPosY,
		      BackExcluderPosZ), "BackExcluder", LogicalBackExcluder, 
		      PhysicalLArg, 0,0);

  LogicalBackExcluder->SetVisAttributes(ColorOfLead);
  // LogicalBackExcluder->SetVisAttributes(G4VisAttributes::GetInvisible());


  //------------------------
  // fcal envelope
  //------------------------
  G4Tubs * SolidFCALEnvelope = 
    new G4Tubs("FCALEnvelopeSolid", FCALEnvelopeRMin, FCALEnvelopeRMax, 
	       FCALEnvelopeLenght, FCALEnvelopeStartPhi, FCALEnvelopeDPhi);
  
  G4LogicalVolume * LogicalFCALEnvelope = 
    new G4LogicalVolume(SolidFCALEnvelope, FCALMaterials->Material("LiquidArgon"),
			"FCALEnvelopeLogical");

  G4RotationMatrix * FCALRotationMatrix = new G4RotationMatrix();
  FCALRotationMatrix->rotateX(FCALEnvelopeRotX);
  //  FCALRotationMatrix->rotateY(FCALEnvelopeRotY);

//  G4VPhysicalVolume *  PhysicalFCALEnvelopp = 
    new G4PVPlacement(FCALRotationMatrix, 
		      G4ThreeVector(FCALEnvelopePosX,FCALEnvelopePosY,FCALEnvelopePosZ)
		      , LogicalFCALEnvelope, "FCALEnvelopePhysical", LogicalLArg, 0,0);

  //LogicalFCALEnvelope->SetVisAttributes(ColorOfIron);
  LogicalFCALEnvelope->SetVisAttributes(G4VisAttributes::GetInvisible());

  //-----------------------------
  // FCAL electromagnetic Module
  //-----------------------------
  EmModule = new FCALEMModule();
  G4LogicalVolume * LogicalFCALEmModule  = EmModule->Construct();

  G4RotationMatrix * EmModuleRot = new G4RotationMatrix();
  EmModuleRot->rotateZ(ModuleRotZ);

//  G4VPhysicalVolume * PhysicalFCALEmModule = 
    new G4PVPlacement(EmModuleRot, 
		      G4ThreeVector(FCALEmModulePosX,FCALEmModulePosY,FCALEmModulePosZ),
		      LogicalFCALEmModule,"FCALEmModulePhysical",LogicalFCALEnvelope,0,0);

       
  //-----------------------------
  // hadronic fcal
  //----------------------------
  HadModule = new FCALHadModule();
  G4LogicalVolume * LogicalFCALHadModule  = HadModule->Construct();

  G4RotationMatrix * HadModuleRot = new G4RotationMatrix();
  HadModuleRot->rotateZ(ModuleRotZ);
  
//  G4VPhysicalVolume * PhysicalFCALHadModule =
    new G4PVPlacement(HadModuleRot, 
		      G4ThreeVector(FCALHadModulePosX,FCALHadModulePosY,FCALHadModulePosZ),
		      LogicalFCALHadModule, "FCALHadModulePhysical",LogicalFCALEnvelope,0,0);
  


  //-------------------------
  // Returning the mother
  //-------------------------

  return LogicalCryostat;

}

