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
//   Author:            Mathieu Fontaine           Rachid Mazini
//                      fontaine@lps.umontreal.ca  Rachid.Mazini@cern.ch
//   Language:          C++
//   Tested on:         g++
//   Prerequisites:     None
//   Purpose:           Header file for FCALFrontVolume.cc, which defines
//                      the  volumes in the testbeam front.
//   Developped:        10-March-2000   M.F.
//
//----------------------------------------------------------------------------


#include "FCALTestbeamSetup.hh"
#include "FCALTestbeamSetupSD.hh"

#include "FCALMaterialConsultant.hh"
#include "FCALCryostatVolumes.hh"

#include "FCALTestbeamSetupSD.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4Threading.hh"

FCALTestbeamSetup::FCALTestbeamSetup() {
#include "FCALTestbeamSetupParameters.input"
}

FCALTestbeamSetup::~FCALTestbeamSetup() {}

G4VPhysicalVolume * FCALTestbeamSetup::Construct()
{
  G4int i=0;

  //-----------------------------
  // construction of materials
  //-----------------------------
  G4cout << "Constructing materials...";
  FCALMaterialConsultant* 
    FCALMaterials = FCALMaterialConsultant::GetInstance();
  G4cout << "... done" << G4endl;

  //-------------------
  // Experimental Hall 
  //-------------------
  G4Box * SolidMother = new G4Box("Mother",MotherSizeX,MotherSizeY,MotherSizeZ);
  G4LogicalVolume * LogicalMother = 
    new G4LogicalVolume(SolidMother,FCALMaterials->Material("Air"),"Mother");
  G4VPhysicalVolume * PhysicalMother =
    new G4PVPlacement(0, G4ThreeVector(),"Mother", LogicalMother, NULL, 0,0);
  
  LogicalMother->SetVisAttributes(G4VisAttributes::GetInvisible());


  //-------------------------------
  //  Scintillators S1, S2 and S3
  //-------------------------------
  G4Box * SolidScintS1andS3 = 
    new G4Box("ScintS1andS3Solid",ScintS1andS3SizeX,  ScintS1andS3SizeY, ScintS1andS3SizeZ);
  G4Box * SolidScintS2 = 
    new G4Box("ScintS2Solid", ScintS2SizeX, ScintS2SizeY, ScintS2SizeZ); 

  G4LogicalVolume * LogicalScintS1andS3 = 
    new G4LogicalVolume(SolidScintS1andS3,FCALMaterials->Material("Polystyrene"),
			"ScintS1andS3Logical");  
  G4LogicalVolume * LogicalScintS2 = 
    new G4LogicalVolume(SolidScintS2, FCALMaterials->Material("Polystyrene"),
			"ScintS2Logical");

  // G4VPhysicalVolume * PhysicalScintS1 = 
    new G4PVPlacement(0, G4ThreeVector(ScintS1_S3PosX, ScintS1_S3PosY, ScintS1PosZ),
   		      "ScintS1Physical",LogicalScintS1andS3,PhysicalMother,0,0);
  // G4VPhysicalVolume * PhysicalScintS3 = 
    new G4PVPlacement(0, G4ThreeVector(ScintS1_S3PosX, ScintS1_S3PosY, ScintS3PosZ),
   		      "ScintS3Physical",LogicalScintS1andS3,PhysicalMother,0,0);  
  // G4VPhysicalVolume * PhysicalScintS2 = 
    new G4PVPlacement(0, G4ThreeVector(ScintS1_S3PosX, ScintS1_S3PosY, ScintS2PosZ),
		      "ScintS2Physical", LogicalScintS2, PhysicalMother,0,0);

  G4VisAttributes * ColorOfScintillator = new G4VisAttributes(G4Colour(0.,0.8,0.));
  LogicalScintS1andS3->SetVisAttributes(ColorOfScintillator);
  LogicalScintS2->SetVisAttributes(ColorOfScintillator);
  

  //-------------------
  //      MWPC's
  //-------------------
  G4Box* SolidMWPC = new G4Box("MWPCSolid",MWPCSizeX,MWPCSizeY,MWPCSizeZ);
  G4LogicalVolume * LogicalMWPC = 
    new G4LogicalVolume(SolidMWPC,FCALMaterials->Material("MWPCArCO2"),"MWPCLogical");
  for(i=0; i<5; i++) 
    {
      // G4VPhysicalVolume * PhysicalMWPC = 
	new G4PVPlacement(0,G4ThreeVector(MWPCPosX,MWPCPosY,MWPCPosZ[i]),
			  "MWPCPhysical", LogicalMWPC, PhysicalMother,0,i+1);
    }
  G4VisAttributes * ColorOfMWPC = new G4VisAttributes(G4Colour(0.,0.,0.5));
  LogicalMWPC->SetVisAttributes(ColorOfMWPC);

  //---------------------------------------
  //  Hole Counter (scintillator + Pb + Al
  //---------------------------------------
  // Scintillator Counter
  G4Box *  SolidHoleCntrScint = 
    new G4Box("ScintSolid", HoleCntrSizeX, HoleCntrSizeY, HoleCntrScintSizeZ);
  G4LogicalVolume * LogicalHoleCntrScint = 
    new G4LogicalVolume(SolidHoleCntrScint, FCALMaterials->Material("Polystyrene"),
			"HoleCntrScintLogical");
  // Hole in scintillator Counter
  G4Tubs * SolidHole = 
    new G4Tubs("HoleSolid", ScintHoleRmin, ScintHoleRmax, ScintHoleLenght,  
	       HoleStartPhi, HoleDPhi);
  G4LogicalVolume * LogicalHole = 
    new G4LogicalVolume(SolidHole, FCALMaterials->Material("Air"), "HoleLogical");
  // G4VPhysicalVolume * PhysicalHoleScint = 
    new G4PVPlacement(0, G4ThreeVector(HolePosX, HolePosY, HolePosZ), LogicalHole, 
		      "HolePhysicalScint", LogicalHoleCntrScint, 0, 0);
  // Scintillator Hole counter placement
  // G4VPhysicalVolume * PhysicalHoleCntrScint =
    new G4PVPlacement(0, 
            G4ThreeVector(HoleCntrScintPosX, HoleCntrScintPosY, HoleCntrScintPosZ), 
            "HoleCntrScintPhysical", LogicalHoleCntrScint, PhysicalMother, 0, 0);

  // Absorber Lead
  G4Box * SolidHoleCntrAbsrb = 
    new G4Box("AbsrbSolid", HoleCntrSizeX, HoleCntrSizeY, HoleCntrAbsrbSizeZ);
  G4LogicalVolume * LogicalHoleCntrPb = 
    new G4LogicalVolume(SolidHoleCntrAbsrb, FCALMaterials->Material("Lead"),
			"HoleCntrPbLoghical");

  //hole in ABsorber, both Lead and Al.
  G4Tubs * SolidHoleAbs = 
    new G4Tubs("HoleSolidAbs", AbsrbHoleRmin, AbsrbHoleRmax, AbsrbHoleLenght, 
	       HoleStartPhi, HoleDPhi);
  G4LogicalVolume * LogicalHoleAbs = 
    new G4LogicalVolume(SolidHoleAbs, FCALMaterials->Material("Air"),"HoleAbsLogical");
  // G4VPhysicalVolume * PhysicalHolePb =
    new G4PVPlacement(0, G4ThreeVector(HolePosX, HolePosY, HolePosZ), LogicalHoleAbs,
		      "HolePbPhysical", LogicalHoleCntrPb, 0, 0);
 
  // Lead Placement
  // G4VPhysicalVolume * PhysicalHoleCntrPb = 
    new G4PVPlacement(0, G4ThreeVector(HoleCntrPbPosX, HoleCntrPbPosY, HoleCntrPbPosZ),
		      "HoleCntrPbPhysical", LogicalHoleCntrPb, PhysicalMother, 0, 0);

  // Absorber Al. 
  G4LogicalVolume * LogicalHoleCntrAl = 
    new G4LogicalVolume(SolidHoleCntrAbsrb,  FCALMaterials->Material("Aluminium"),
			"HoleCntrAlLogical");
  // G4VPhysicalVolume * PhysicalHoleAl =
    new G4PVPlacement(0, G4ThreeVector(HolePosX, HolePosY, HolePosZ), LogicalHoleAbs,
		      "HoleAlPhysical", LogicalHoleCntrAl, 0, 0);
  // G4VPhysicalVolume * PhysicalHoleCntrAl = 
    new G4PVPlacement(0, G4ThreeVector(HoleCntrAlPosX, HoleCntrAlPosY, HoleCntrAlPosZ),
		      "HoleCntrAlPhysical", LogicalHoleCntrAl, PhysicalMother, 0, 0);
  
   LogicalHoleCntrScint->SetVisAttributes(ColorOfScintillator);

   G4VisAttributes * ColorOfLead = new G4VisAttributes(G4Colour(0.5,0.5,0.8));
   G4VisAttributes * ColorOfAlu = new G4VisAttributes(G4Colour(0.5,0.5,0.3));
   LogicalHoleCntrPb->SetVisAttributes(ColorOfLead);
   LogicalHoleCntrAl->SetVisAttributes(ColorOfAlu);

   G4VisAttributes * ColorOfAir = new G4VisAttributes(G4Colour(1.,1.,1.));
   LogicalHole->SetVisAttributes(ColorOfAir);
   LogicalHoleAbs->SetVisAttributes(ColorOfAir);


   //-------------------
   //   Lead Wall
   //-------------------
   G4Box * SolidLeadWall = 
     new G4Box("LeadWallSolid", LeadWallSizeX, LeadWallSizeY, LeadWallSizeZ);
   G4LogicalVolume * LogicalLeadWall = 
     new G4LogicalVolume(SolidLeadWall, FCALMaterials->Material("Lead"), 
			 "LeadWallLogical");    

   G4Box * SolidSlitPb = new G4Box("SlitPb", LeadWallSlitSizeX, LeadWallSlitSizeY, 
				 LeadWallSlitSizeZ);
   G4LogicalVolume * LogicalSlitPb = 
     new G4LogicalVolume(SolidSlitPb, FCALMaterials->Material("Air"), "SlitPbLogical");
   // G4VPhysicalVolume * PhysicalSlitPb = 
     new G4PVPlacement(0, G4ThreeVector(), LogicalSlitPb, "SlitPbPhysical", LogicalLeadWall, 0, 0);

   // G4VPhysicalVolume * PhysicalLeadWall = 
     new G4PVPlacement(0, G4ThreeVector(LeadWallPosX,LeadWallPosY,LeadWallPosZ),
		       "LeadWallPhysical", LogicalLeadWall, PhysicalMother, 0, 0);

   LogicalLeadWall->SetVisAttributes(ColorOfLead);
   LogicalSlitPb->SetVisAttributes(ColorOfAir);


    //-------------------
   //   Iron Wall
   //-------------------
   G4Box * SolidIronWall = 
     new G4Box("IronWallSolid", IronWallSizeX, IronWallSizeY, IronWallSizeZ);
   G4LogicalVolume * LogicalIronWall = 
     new G4LogicalVolume(SolidIronWall, FCALMaterials->Material("Iron"), 
			 "IronWallLogical");    

   G4Box * SolidSlitFe = new G4Box("SlitFe", IronWallSlitSizeX, IronWallSlitSizeY, 
				 IronWallSlitSizeZ);
   G4LogicalVolume * LogicalSlitFe = 
     new G4LogicalVolume(SolidSlitFe, FCALMaterials->Material("Air"), "SlitFeLogical");
   // G4VPhysicalVolume * PhysicalSlitFe = 
     new G4PVPlacement(0, G4ThreeVector(), LogicalSlitFe, "SlitFePhysical", LogicalIronWall, 0, 0);

   // G4VPhysicalVolume * PhysicalIronWall = 
     new G4PVPlacement(0, G4ThreeVector(IronWallPosX,IronWallPosY,IronWallPosZ),
		       "IronWallPhysical", LogicalIronWall, PhysicalMother, 0, 0);

   G4VisAttributes * ColorOfIron = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
   LogicalIronWall->SetVisAttributes(ColorOfIron);
   LogicalSlitFe->SetVisAttributes(ColorOfAir);

   //----------------
   // Tail Catcher
   //----------------
   G4Box * SolidBigScint = 
     new G4Box("BigSolidScint",BigScintSizeX, BigScintSizeY, ScintSizeZ);
   G4LogicalVolume * LogicalBigScint =
     new G4LogicalVolume(SolidBigScint,  FCALMaterials->Material("Polystyrene"),
			 "BigScintLogical");
   
   G4Box * SolidSmallScint = 
     new G4Box("SmallSolidScint",SmallScintSizeX, SmallScintSizeY, ScintSizeZ);
   G4LogicalVolume * LogicalSmallScint =
     new G4LogicalVolume(SolidSmallScint,  FCALMaterials->Material("Polystyrene"),
			 "SmallScintLogical");

   for( i=0; i<(NBigScint+NSmallScint); i++)
     { 
       if(i<NBigScint)  
	 { // G4VPhysicalVolume * PhysicalBigScint =
	     new G4PVPlacement(0, G4ThreeVector(ScintPosX, ScintPosY, ScintPosZ[i]),
			       "BigScintPhysical", LogicalBigScint, PhysicalMother, 
			       0, i+1); 
	 }
       else
	 { // G4VPhysicalVolume * PhysicalSmallScint =
	     new G4PVPlacement(0, G4ThreeVector(ScintPosX, ScintPosY, ScintPosZ[i]),
			       "SmallScintPhysical", LogicalSmallScint, PhysicalMother,
			       0, i+1);
	 }
     }
   LogicalBigScint->SetVisAttributes(ColorOfScintillator);
   LogicalSmallScint->SetVisAttributes(ColorOfScintillator);

   
   G4Box * SolidBigIron = 
     new G4Box("BigSolidIron",BigIronSizeX, BigIronSizeY, IronSizeZ);
   G4LogicalVolume * LogicalBigIron =
     new G4LogicalVolume(SolidBigIron,  FCALMaterials->Material("Polystyrene"),
			 "BigIronLogical");

   G4Box * SolidSmallIron = 
     new G4Box("SmallSolidIron",SmallIronSizeX, SmallIronSizeY, IronSizeZ);
   G4LogicalVolume * LogicalSmallIron =
     new G4LogicalVolume(SolidSmallIron,  FCALMaterials->Material("Iron"),
			 "SmallIronLogical");

   for( i=0; i<(NBigIron+NSmallIron); i++)
     { 
       if(i<NBigIron)  
	 { // G4VPhysicalVolume * PhysicalBigIron =
	     new G4PVPlacement(0, G4ThreeVector(IronPosX, IronPosY, IronPosZ[i]),
			       "BigIronPhysical", LogicalBigIron, PhysicalMother, 
			       0, i+1); 
	 }
       else
	 { // G4VPhysicalVolume * PhysicalSmallIron =
	     new G4PVPlacement(0, G4ThreeVector(IronPosX, IronPosY, IronPosZ[i]),
			       "SmallIronPhysical", LogicalSmallIron, PhysicalMother,
			       0, i+1);
	 }
     }
   LogicalBigIron->SetVisAttributes(ColorOfIron);
   LogicalSmallIron->SetVisAttributes(ColorOfIron);
  
   //-------------------------
   //  Concrete Walls A and B
   //-------------------------
   G4Box * SolidConcWall = 
     new G4Box("ConcWallSolid", ConcWallSizeX, ConcWallSizeY, ConcWallSizeZ);
   G4LogicalVolume * LogicalConcWallA =
     new G4LogicalVolume(SolidConcWall, FCALMaterials->Material("ShieldingConcrete"),
			 "ConcWallALogical");
   G4VPhysicalVolume * PhysicalConcWallA = 
     new G4PVPlacement(0, G4ThreeVector(ConcWallPosX, ConcWallPosY, ConcWallAPosZ),
		       "ConcWallAPhysical", LogicalConcWallA, PhysicalMother, 0, 0);

   G4LogicalVolume * LogicalConcWallB =
     new G4LogicalVolume(SolidConcWall, FCALMaterials->Material("ShieldingConcrete"),
			 "ConcWallBLogical");
   // G4VPhysicalVolume * PhysicalConcWallB = 
     new G4PVPlacement(0, G4ThreeVector(ConcWallPosX, ConcWallPosY, ConcWallBPosZ),
		       "ConcWallBPhysical", LogicalConcWallB, PhysicalMother, 0, 0);

    G4Box * SolidConcWallIns = 
     new G4Box("ConcWallInsSolid", ConcWallInsSizeX,ConcWallInsSizeY,ConcWallInsSizeZ);
    G4LogicalVolume * LogicalConcWallIns =
      new G4LogicalVolume(SolidConcWallIns, FCALMaterials->Material("Iron"),
			  "LogicalConcWallIns");
    // G4VPhysicalVolume * PhysicalConcWallIns =
      new G4PVPlacement(0, G4ThreeVector(), "ConcWallInsPhysical", LogicalConcWallIns, PhysicalConcWallA, 0, 0);

   G4VisAttributes * ColorOfConcrete = new G4VisAttributes(G4Colour(0.,0.,0.));
   LogicalConcWallA->SetVisAttributes(ColorOfConcrete);
   LogicalConcWallB->SetVisAttributes(ColorOfConcrete);
   LogicalConcWallIns->SetVisAttributes(ColorOfIron);

   //------------------
   //  Muon Counter
   //-------------------
   G4Box * SolidMuContr = 
     new G4Box("MuContrSolid", MuCntrSIzeX, MuCntrSIzeY, MuCntrSIzeZ);
   G4LogicalVolume * LogicalMuContr =
     new G4LogicalVolume(SolidMuContr, FCALMaterials->Material("Polystyrene"),
			 "MuContrLogical");
   // G4VPhysicalVolume * PhysicalMuContr =
     new G4PVPlacement(0, G4ThreeVector(MuCntrPosX, MuCntrPosY, MuCntrPosZ),
		       "MuContrPhyiscal", LogicalMuContr, PhysicalMother, 0, 0);

   LogicalMuContr->SetVisAttributes(ColorOfScintillator);

  //-----------------
  // cryostat
  //-----------------


  G4RotationMatrix*  CryostatRotationMatrix = 
    new G4RotationMatrix();

  //    new G4RotationMatrix(1.,0.,0.,0.,0.,-1.,0.,1.,0.);                       

  //  Theta(...)    90.0000  180.0000   90.0000
  //  Phi(...)       0.0000    0.0000   90.0000
  //
  //  Matrix(...) |  1.0000  0.0000  0.0000 |
  //              |  0.0000  0.0000 -1.0000 |
  //              |  0.0000  1.0000  0.0000 |
  //
  // How to input?
  		
  CryostatRotationMatrix->rotateX(90*deg);				     

  FCALCryostatVolumes * CryostatVolumes = new FCALCryostatVolumes();

  G4LogicalVolume * theCryostatVolumes = CryostatVolumes->Construct();

  // G4VPhysicalVolume * PhysiCryostatVolumes = 
    new G4PVPlacement(CryostatRotationMatrix, 
		      G4ThreeVector(CryostatPosX,CryostatPosY,CryostatPosZ),"CryostatVolumes"
		      , theCryostatVolumes, PhysicalMother, 0,0);


  return PhysicalMother;

}

void FCALTestbeamSetup::ConstructSDandField()
{
    //-----------------------
    // Senstive detectors
    //-----------------------
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    const G4String detName = "FCALTB/TBSetupSD";
    FCALTestbeamSetupSD* FCALTBSetupSD = static_cast<FCALTestbeamSetupSD*>(SDman->FindSensitiveDetector(detName));
    if(!FCALTBSetupSD)
    {
        FCALTBSetupSD = new FCALTestbeamSetupSD(detName);
        SDman->AddNewDetector(FCALTBSetupSD);
    }
    SetSensitiveDetector("ScintS1andS3Logical",FCALTBSetupSD);
    SetSensitiveDetector("ScintS2Logical",FCALTBSetupSD);
    SetSensitiveDetector("MWPCLogical",FCALTBSetupSD);
    SetSensitiveDetector("HoleCntrScintLogical",FCALTBSetupSD);
    SetSensitiveDetector("HoleCntrPbLoghical",FCALTBSetupSD);
    SetSensitiveDetector("HoleCntrAlLogical",FCALTBSetupSD);
    
    SetSensitiveDetector("LeadWallLogical",FCALTBSetupSD);
    SetSensitiveDetector("IronWallLogical",FCALTBSetupSD);
    SetSensitiveDetector("BigScintLogical",FCALTBSetupSD);
    SetSensitiveDetector("SmallScintLogical",FCALTBSetupSD);
    
    SetSensitiveDetector("BigIronLogical",FCALTBSetupSD);
    SetSensitiveDetector("SmallIronLogical",FCALTBSetupSD);
    SetSensitiveDetector("ConcWallALogical",FCALTBSetupSD);
    SetSensitiveDetector("ConcWallBLogical",FCALTBSetupSD);
    SetSensitiveDetector("LogicalConcWallIns",FCALTBSetupSD);
    SetSensitiveDetector("MuContrLogical",FCALTBSetupSD);
}

