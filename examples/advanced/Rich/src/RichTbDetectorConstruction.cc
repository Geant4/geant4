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
// Rich advanced example for Geant4
// RichTbDetectorConstruction.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "RichTbMaterialParameters.hh"
#include "RichTbGeometryParameters.hh"
#include "RichTbDetectorConstruction.hh"
#include "RichTbHall.hh"
#include "RichTbMaterial.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "RichTbSD.hh"



RichTbDetectorConstruction::RichTbDetectorConstruction(){; }
RichTbDetectorConstruction::RichTbDetectorConstruction(RichTbRunConfig* RConfig){

  runConfiguration=RConfig;
 }
RichTbDetectorConstruction::~RichTbDetectorConstruction(){;}

G4VPhysicalVolume* RichTbDetectorConstruction::Construct(){

  RichTbRunConfig* rConfig= runConfiguration;

  InitializeRichTbMaterial();
  HistoRichTbMaterialProperties(rConfig);

  rMaterial = new RichTbMaterial(rConfig);

  InitializeRichTbGeometry();


  rTbHall = new RichTbHall(rMaterial);
   
  rTbComponent = new RichTbComponent(rMaterial,rTbHall,rConfig,true);
   
  rTbPhotoDetector = new RichTbPhotoDetector(rMaterial,rTbComponent,rConfig,true);


   //Now for the sensitive Detector
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String HpdSDname = "RichTbSD";
  G4String RichTbHColname="RichTbHitsCollection";
  
   
  RichTbSD* HpdSD = new RichTbSD( HpdSDname,
				  NumberOfHpds,NumberOfSiDetSectors,
				  NumberOfPadHpdSiPixels,
				  RichTbHColname );
  
  
  
   G4String ROgeometryName = "RichTbROGeom";
   
   rTbROGeom = new RichTbROGeometry(ROgeometryName,rMaterial, rConfig);
   
   rTbROGeom -> BuildROGeometry();
   HpdSD ->SetROgeometry( rTbROGeom);
   SDman->AddNewDetector( HpdSD );
   G4int NumHpd = NumberOfHpds;

   for (G4int ihpd=0 ; ihpd < NumHpd ; ihpd++ ) { 
     G4int NumSiSect=rTbPhotoDetector-> getRichHPD(ihpd) ->
       getNumSectInHpd(); 
     for (G4int isec=0; isec < NumSiSect; isec++) {
       G4LogicalVolume* RichTbSiDetSect_LV =rTbPhotoDetector -> 
	 getRichHPD(ihpd)-> getRichHpdSiDetSect(isec)->
	 getHpdSiSectLogicalVolume();
       RichTbSiDetSect_LV->SetSensitiveDetector( HpdSD );  
	}
      }


   rTbGraphics = new RichTbGraphics(rTbHall,rTbComponent,
                                rTbPhotoDetector,rTbROGeom, rConfig );


  return rTbHall->getRichTbHallPhysicalVolume();
}






