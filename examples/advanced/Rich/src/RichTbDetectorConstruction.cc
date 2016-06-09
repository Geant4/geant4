//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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






