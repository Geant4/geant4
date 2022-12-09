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
// 
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

#include "TETDetectorConstruction.hh"

#include "G4VisAttributes.hh"

TETDetectorConstruction::TETDetectorConstruction(TETModelImport* _tetData)
:fWorldPhysical(nullptr), fContainer_logic(nullptr), fTetData(_tetData), fTetLogic(nullptr)
{
 // initialisation of the variables for phantom information
 fPhantomSize     = fTetData -> GetPhantomSize();
 fPhantomBoxMin   = fTetData -> GetPhantomBoxMin();
 fPhantomBoxMax   = fTetData -> GetPhantomBoxMax();
 fNOfTetrahedrons = fTetData -> GetNumTetrahedron();
}

TETDetectorConstruction::~TETDetectorConstruction()
{
  delete fTetData;
}

G4VPhysicalVolume* TETDetectorConstruction::Construct()
{
 SetupWorldGeometry();
 ConstructPhantom();
 PrintPhantomInformation();
 return fWorldPhysical;
}

void TETDetectorConstruction::SetupWorldGeometry()
{
 // Define the world box (size: 10*10*10 m3)
 //
 G4double worldXYZ = 10. * m;
 G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

 G4VSolid* worldSolid = new G4Box("worldSolid", worldXYZ/2, worldXYZ/2, worldXYZ/2);

 auto* worldLogical = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");

 fWorldPhysical = new G4PVPlacement(nullptr,G4ThreeVector(), worldLogical,"worldPhysical", nullptr, false,0,false);

 // Define the phantom container (10-cm margins from the bounding box of phantom)
 //
 auto* containerSolid = new G4Box("phantomBox", fPhantomSize.x()/2 + 10.*cm,
					           fPhantomSize.y()/2 + 10.*cm,
						   fPhantomSize.z()/2 + 10.*cm);

 fContainer_logic = new G4LogicalVolume(containerSolid, vacuum, "phantomLogical");

 new G4PVPlacement(nullptr, G4ThreeVector(), fContainer_logic, "PhantomPhysical",
	           worldLogical, false, 0);

 fContainer_logic->SetOptimisation(TRUE);
 fContainer_logic->SetSmartless( 0.5 ); // for optimization (default=2)
}

void TETDetectorConstruction::ConstructPhantom()
{
 // Define the tetrahedral mesh phantom as a parameterised geometry
 //
 // solid and logical volume to be used for parameterised geometry
 G4VSolid* tetraSolid = new G4Tet("TetSolid", G4ThreeVector(),
			           G4ThreeVector(1.*cm,0,0),
			           G4ThreeVector(0,1.*cm,0),
			           G4ThreeVector(0,0,1.*cm));

  G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  fTetLogic = new G4LogicalVolume(tetraSolid, vacuum, "TetLogic");

  // physical volume (phantom) constructed as parameterised geometry
  new G4PVParameterised("wholePhantom",fTetLogic,fContainer_logic,
			  kUndefined, fTetData->GetNumTetrahedron(),
			 new TETParameterisation(fTetData));
}

void TETDetectorConstruction::ConstructSDandField()
{
 // Define detector (Phantom SD) and scorer (eDep)
 //
 G4SDManager* pSDman = G4SDManager::GetSDMpointer();
 G4String phantomSDname = "PhantomSD";

 // MultiFunctional detector
 auto* MFDet = new G4MultiFunctionalDetector(phantomSDname);
 pSDman->AddNewDetector( MFDet );

 // scorer for energy depositon in each organ
 MFDet->RegisterPrimitive(new TETPSEnergyDeposit("eDep", fTetData));

 // attach the detector to logical volume for parameterised geometry (phantom geometry)
 SetSensitiveDetector(fTetLogic, MFDet);
}

void TETDetectorConstruction::PrintPhantomInformation()
{
 // print brief information on the imported phantom
 G4cout<< G4endl;
 G4cout.precision(3);
 G4cout<<"   Phantom name               "<<fTetData->GetPhantomName() << " TET phantom"<<G4endl;
 G4cout<<"   Phantom size               "<<fPhantomSize.x()<<" * "<<fPhantomSize.y()<<" * "<<fPhantomSize.z()<<" mm3"<<G4endl;
 G4cout<<"   Phantom box position (min) "<<fPhantomBoxMin.x()<<" mm, "<<fPhantomBoxMin.y()<<" mm, "<<fPhantomBoxMin.z()<<" mm"<<G4endl;
 G4cout<<"   Phantom box position (max) "<<fPhantomBoxMax.x()<<" mm, "<<fPhantomBoxMax.y()<<" mm, "<<fPhantomBoxMax.z()<<" mm"<<G4endl;
 G4cout<<"   Number of tetrahedrons     "<<fNOfTetrahedrons<<G4endl<<G4endl;
}
