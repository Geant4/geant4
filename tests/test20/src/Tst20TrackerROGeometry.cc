// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20TrackerROGeometry.cc,v 1.1 2001-05-24 19:49:30 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20TrackerROGeometry class ------
// ************************************************************

#include "Tst20TrackerROGeometry.hh"
#include "Tst20DummySD.hh"
#include "Tst20DetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
 

Tst20TrackerROGeometry::Tst20TrackerROGeometry()
  : G4VReadOutGeometry()
{
}
Tst20TrackerROGeometry::Tst20TrackerROGeometry(G4String aString,Tst20DetectorConstruction* Tst20DC)
  :Tst20Detector(Tst20DC), G4VReadOutGeometry(aString)
{
}

Tst20TrackerROGeometry::Tst20TrackerROGeometry(G4String aString)
  : G4VReadOutGeometry(aString)
{
}

Tst20TrackerROGeometry::~Tst20TrackerROGeometry()
{
}

G4VPhysicalVolume* Tst20TrackerROGeometry::Build()
{
  // A dummy material is used to fill the volumes of the readout geometry.
  // ( It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)
  G4Material* dummyMat  = new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);
  
  //Builds the ReadOut World:

  G4double WorldSizeXY = Tst20Detector->GetWorldSizeXY();
  G4double WorldSizeZ  = Tst20Detector->GetWorldSizeZ();
  
  G4Box* ROWorldBox = new 
    G4Box("ROWorldBox",WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2); 
  
  G4LogicalVolume* ROWorldLog = new G4LogicalVolume(ROWorldBox, dummyMat,
						    "ROWorldLogical");
  G4PVPlacement* ROWorldPhys = 
    new G4PVPlacement(0,G4ThreeVector(),"ROWorldPhysical",
		      ROWorldLog,0,false,0);
  
  // Payload RO volume:
  
  G4double PayloadSizeXY = Tst20Detector->GetPayloadSizeXY();
  G4double PayloadSizeZ  = Tst20Detector->GetPayloadSizeZ();
  
  G4VSolid* solidPayloadRO
    = new G4Box("Payload RO",		
		PayloadSizeXY/2,
		PayloadSizeXY/2,
		PayloadSizeZ/2);
  
  G4LogicalVolume* logicPayloadRO = new 
    G4LogicalVolume(solidPayloadRO,dummyMat,"Payload RO",0,0,0);  
  
  G4VPhysicalVolume* physiPayloadRO =
    new G4PVPlacement(0, G4ThreeVector(), 
		      "Payload RO", logicPayloadRO,ROWorldPhys,false, 0);
  
  // -------------------------------
  // Tracker readout division:
  // -------------------------------
  // TRK Layers of Silicon MicroStrips
  
  
  G4double TKRSizeXY = Tst20Detector->GetTKRSizeXY();
  G4double TKRSizeZ  = Tst20Detector->GetTKRSizeZ();
  G4double ACDTKRDistance  = Tst20Detector->GetACDTKRDistance();
  G4double ACDSizeZ  = Tst20Detector->GetACDSizeZ();
  
  G4VSolid*   ROsolidTKR =
    new G4Box("ReadOutTKR", TKRSizeXY/2,TKRSizeXY/2,TKRSizeZ/2); 
  
  G4LogicalVolume*  ROlogicTKR =
    new G4LogicalVolume(ROsolidTKR,dummyMat, "ReadOutTKR",0,0,0);	
  
  
  G4VPhysicalVolume* ROphysiTKR = 
    new G4PVPlacement(0, G4ThreeVector(0,0,-PayloadSizeZ/2+ACDSizeZ+
				       ACDTKRDistance+TKRSizeZ/2),
		      "ReadOutTKR",ROlogicTKR,physiPayloadRO,
		      false, 0);		
  
  // TKR Layers
  
  
  G4double TKRSiliconThickness  =
    Tst20Detector->GetTKRSiliconThickness();
  G4int NbOfTKRLayers =  Tst20Detector->GetNbOfTKRLayers();
  G4double TKRLayerDistance =  Tst20Detector->GetTKRLayerDistance();
  G4double TKRSupportThickness =  Tst20Detector->GetTKRSupportThickness();
  G4double ConverterThickness =  Tst20Detector->GetConverterThickness();
  
  G4VSolid*  solidTKRDetectorRO = new G4Box
    ("TKRDetectorRO",TKRSizeXY/2,TKRSizeXY/2,TKRSiliconThickness/2); 
  
  G4LogicalVolume* logicTKRDetectorRO =
    new G4LogicalVolume(solidTKRDetectorRO,dummyMat, "TKRDetectorRO",0,0,0);
  
  G4int i=0;
  G4VPhysicalVolume* physiTKRDetectorRO = 0;
  
  for (i = 0; i < NbOfTKRLayers; i++)
    {
      physiTKRDetectorRO = 
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-TKRSizeZ/2
					  +TKRSupportThickness+
                                          +ConverterThickness+
                                          +TKRSiliconThickness/2 
					  +(i)*TKRLayerDistance),
			  "TKRDetectorRO",		
			  logicTKRDetectorRO,
			  ROphysiTKR,
			  false,	
			  i);	
      
    }
  

  // Silicon Tiles 
  // some problems with the RO tree
  
  G4double TKRActiveTileXY =  Tst20Detector->GetTKRActiveTileXY();
  G4double TKRActiveTileZ  =  Tst20Detector->GetTKRActiveTileZ();  
  

  
  G4VSolid * solidTKRActiveTileRO = new
    G4Box("Active Tile", TKRActiveTileXY/2,TKRActiveTileXY/2,TKRActiveTileZ/2);

  G4LogicalVolume* logicTKRActiveTileRO = 
    new G4LogicalVolume(solidTKRActiveTileRO, dummyMat,"Active Tile",0,0,0);
  
  G4int j=0;
  G4int k=0;
  
  G4int NbOfTKRTiles = Tst20Detector->GetNbOfTKRTiles();
  G4double SiliconGuardRing = Tst20Detector->GetSiliconGuardRing();
  G4double TKRTileDistance =  Tst20Detector->GetTKRTileDistance();  
  
  G4VPhysicalVolume* physiTKRActiveTileRO = 0;
  
  G4double x=0.;
  G4double y=0.;
  G4double z=0.;
  
  for (i=0;i< NbOfTKRTiles; i++)
    { 
      for (j=0;j< NbOfTKRTiles; j++)
	{
	  k = i*NbOfTKRTiles + j;
	  
	  x = -TKRSizeXY/2+TKRTileDistance+SiliconGuardRing+TKRActiveTileXY/2+
	    (j)*((2*SiliconGuardRing)+TKRTileDistance+TKRActiveTileXY);
          y = -TKRSizeXY/2+TKRTileDistance+SiliconGuardRing+TKRActiveTileXY/2+
	    (i)*((2*SiliconGuardRing)+TKRTileDistance+TKRActiveTileXY);
          z = 0.;

	  physiTKRActiveTileRO =
	    new G4PVPlacement(0,
			      G4ThreeVector(x,y,z),
			      "Active Tile",		
			      logicTKRActiveTileRO,
			      physiTKRDetectorRO,
			      false,	
			      k);	

	  
	}
    }
  
  
  // Silicon Pixels 

  G4double TKRPixelXY=0.;
  G4double TKRPixelZ=0.;
  
  TKRPixelXY = Tst20Detector->GetTKRSiliconPitch();
  TKRPixelZ  = Tst20Detector->GetTKRSiliconThickness();
  
  G4int NbOfTKRPixels  = Tst20Detector->GetNbOfTKRPixels();
  
  G4VSolid* solidTKRPixel = new G4Box("Pixel",			
				      TKRPixelXY/2,TKRPixelXY/2,
				      TKRPixelZ/2); 
  
  G4LogicalVolume* logicTKRPixel = 
    new G4LogicalVolume(solidTKRPixel,dummyMat,"Pixel",0,0,0);	 
  
  
  G4VPhysicalVolume* physiTKRPixel = 0;

  for (i=0;i< NbOfTKRPixels; i++)
    {
      for (j=0;j< NbOfTKRPixels; j++)
	{  
	  k = i*NbOfTKRPixels + j;
	  physiTKRPixel = new G4PVPlacement
	    (0,G4ThreeVector(-TKRActiveTileXY/2 +TKRPixelXY/2 +
			     (i)*TKRPixelXY, 
			     -TKRActiveTileXY/2 +TKRPixelXY/2 +
			     (j)*TKRPixelXY, 0.),
	     "Pixel",		
	     logicTKRPixel,
	     physiTKRActiveTileRO,
	     false,	
	     i);	
	  
	} 
    }
  

  
  //Flags the strip as sensitive .The pointer here serves
  // as a flag only to check for sensitivity.
  // (Could we make it by a simple cast of a non-NULL value ?)
  
  
  Tst20DummySD * dummySensi = new Tst20DummySD;
  
  logicTKRPixel->SetSensitiveDetector(dummySensi);
  
  return ROWorldPhys;
}










