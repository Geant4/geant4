// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelTrackerROGeometry.cc,v 1.1 2001-03-05 13:58:23 flongo Exp $
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
//      ------------ GammaRayTelTrackerROGeometry class ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelTrackerROGeometry.hh"
#include "GammaRayTelDummySD.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
 

GammaRayTelTrackerROGeometry::GammaRayTelTrackerROGeometry()
  : G4VReadOutGeometry()
{
}
GammaRayTelTrackerROGeometry::GammaRayTelTrackerROGeometry(G4String aString,GammaRayTelDetectorConstruction* GammaRayTelDC)
  :GammaRayTelDetector(GammaRayTelDC), G4VReadOutGeometry(aString)
{
}

GammaRayTelTrackerROGeometry::GammaRayTelTrackerROGeometry(G4String aString)
  : G4VReadOutGeometry(aString)
{
}

GammaRayTelTrackerROGeometry::~GammaRayTelTrackerROGeometry()
{
}

G4VPhysicalVolume* GammaRayTelTrackerROGeometry::Build()
{
  // A dummy material is used to fill the volumes of the readout geometry.
  // ( It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)
  G4Material* dummyMat  = new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);
  
  //Builds the ReadOut World:

  G4double WorldSizeXY = GammaRayTelDetector->GetWorldSizeXY();
  G4double WorldSizeZ  = GammaRayTelDetector->GetWorldSizeZ();
  
  G4Box* ROWorldBox = new 
    G4Box("ROWorldBox",WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2); 
  
  G4LogicalVolume* ROWorldLog = new G4LogicalVolume(ROWorldBox, dummyMat,
						    "ROWorldLogical");
  G4PVPlacement* ROWorldPhys = 
    new G4PVPlacement(0,G4ThreeVector(),"ROWorldPhysical",
		      ROWorldLog,0,false,0);
  
  // Payload RO volume:
  
  G4double PayloadSizeXY = GammaRayTelDetector->GetPayloadSizeXY();
  G4double PayloadSizeZ  = GammaRayTelDetector->GetPayloadSizeZ();
  
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
 
  
  G4double TKRSizeXY = GammaRayTelDetector->GetTKRSizeXY();
  G4double TKRSizeZ  = GammaRayTelDetector->GetTKRSizeZ();
  G4double CALSizeZ  = GammaRayTelDetector->GetCALSizeZ();
  G4double CALTKRDistance  = GammaRayTelDetector->GetCALTKRDistance();
  
  G4VSolid*   ROsolidTKR =
    new G4Box("ReadOutTKR", TKRSizeXY/2,TKRSizeXY/2,TKRSizeZ/2); 
  
  G4LogicalVolume*  ROlogicTKR =
    new G4LogicalVolume(ROsolidTKR,dummyMat, "ReadOutTKR",0,0,0);	
						     
  
  G4VPhysicalVolume* ROphysiTKR = 
    new G4PVPlacement(0, G4ThreeVector(0,0,-PayloadSizeZ/2+CALSizeZ+
				       CALTKRDistance+TKRSizeZ/2),
		      "ReadOutTKR",ROlogicTKR,physiPayloadRO,
		      false, 0);		
  
  // TKR Layers
  
  
  G4double TKRSiliconThickness  =
    GammaRayTelDetector->GetTKRSiliconThickness();
  G4int NbOfTKRLayers =  GammaRayTelDetector->GetNbOfTKRLayers();
  G4double TKRLayerDistance =  GammaRayTelDetector->GetTKRLayerDistance();
  G4double TKRViewsDistance =  GammaRayTelDetector->GetTKRViewsDistance();  
  
  G4VSolid*  solidTKRDetectorYRO = new G4Box
    ("TKRDetectorYRO",TKRSizeXY/2,TKRSizeXY/2,TKRSiliconThickness/2); 
  
  G4LogicalVolume* logicTKRDetectorYRO =
    new G4LogicalVolume(solidTKRDetectorYRO,dummyMat, "TKRDetectorYRO",0,0,0);

  G4VSolid*  solidTKRDetectorXRO = new G4Box
    ("TKRDetectorXRO",TKRSizeXY/2,TKRSizeXY/2,TKRSiliconThickness/2); 
  
  G4LogicalVolume* logicTKRDetectorXRO =
    new G4LogicalVolume(solidTKRDetectorXRO,dummyMat, "TKRDetectorXRO",0,0,0);
  G4int i=0;
  G4VPhysicalVolume* physiTKRDetectorXRO = 0;
  G4VPhysicalVolume* physiTKRDetectorYRO = 0;
  
  
  for (i = 0; i < NbOfTKRLayers; i++)
    {
      
      physiTKRDetectorYRO = 
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-TKRSizeZ/2
					  +TKRSiliconThickness/2 
					  +(i)*TKRLayerDistance),
			  "TKRDetectorYRO",		
			  logicTKRDetectorYRO,
			  ROphysiTKR,
			  false,	
			  i);	
      
      physiTKRDetectorXRO = 
	new G4PVPlacement(0,G4ThreeVector(0.,0.,
					  -TKRSizeZ/2+
					  TKRSiliconThickness/2 +
					  TKRViewsDistance+
					  TKRSiliconThickness+
					  (i)*TKRLayerDistance),
			  "TKRDetectorXRO",		
			  logicTKRDetectorXRO,
			  ROphysiTKR,
			  false,	
			  i);	
    }
  

  // Silicon Tiles 
  // some problems with the RO tree

  G4double TKRActiveTileXY =  GammaRayTelDetector->GetTKRActiveTileXY();
  G4double TKRActiveTileZ  =  GammaRayTelDetector->GetTKRActiveTileZ();  

  G4VSolid * solidTKRActiveTileXRO = new
    G4Box("Active Tile X", TKRActiveTileXY/2,TKRActiveTileXY/2,TKRActiveTileZ/2);

  G4VSolid * solidTKRActiveTileYRO = new
    G4Box("Active Tile Y", TKRActiveTileXY/2,TKRActiveTileXY/2,TKRActiveTileZ/2);
  
  
  G4LogicalVolume* logicTKRActiveTileXRO = 
    new G4LogicalVolume(solidTKRActiveTileXRO, dummyMat,"Active Tile",0,0,0);

  G4LogicalVolume* logicTKRActiveTileYRO = 
    new G4LogicalVolume(solidTKRActiveTileYRO, dummyMat,"Active Tile",0,0,0);
    
  G4int j=0;
  G4int k=0;
  
  G4int NbOfTKRTiles = GammaRayTelDetector->GetNbOfTKRTiles();
  G4double SiliconGuardRing = GammaRayTelDetector->GetSiliconGuardRing();
  G4double TilesSeparation = GammaRayTelDetector->GetTilesSeparation();
  
  G4VPhysicalVolume* physiTKRActiveTileXRO = 0;
  G4VPhysicalVolume* physiTKRActiveTileYRO = 0;

  G4double x=0.;
  G4double y=0.;
  G4double z=0.;

  for (i=0;i< NbOfTKRTiles; i++)
    { 
      for (j=0;j< NbOfTKRTiles; j++)
	{
	  k = i*NbOfTKRTiles + j;
	  
	  x = -TKRSizeXY/2+TilesSeparation+SiliconGuardRing+TKRActiveTileXY/2+
	    (j)*((2*SiliconGuardRing)+TilesSeparation+TKRActiveTileXY);
          y = -TKRSizeXY/2+TilesSeparation+SiliconGuardRing+TKRActiveTileXY/2+
	    (i)*((2*SiliconGuardRing)+TilesSeparation+TKRActiveTileXY);
          z = 0.;

	  physiTKRActiveTileXRO =
	    new G4PVPlacement(0,
			      G4ThreeVector(x,y,z),
			      "Active Tile X",		
			      logicTKRActiveTileXRO,
			      physiTKRDetectorXRO,
			      false,	
			      k);	

	  x = -TKRSizeXY/2+TilesSeparation+SiliconGuardRing+TKRActiveTileXY/2+
	    (i)*((2*SiliconGuardRing)+TilesSeparation+TKRActiveTileXY);
          y = -TKRSizeXY/2+TilesSeparation+SiliconGuardRing+TKRActiveTileXY/2+
	    (j)*((2*SiliconGuardRing)+TilesSeparation+TKRActiveTileXY);
          z = 0.;

      	  physiTKRActiveTileYRO =
	    new G4PVPlacement(0,
			      G4ThreeVector(x,y,z),
			      "Active Tile Y",		
			      logicTKRActiveTileYRO,
			      physiTKRDetectorYRO,
			      false,	
			      k);
	  	
	}
    }
  
  
  // Silicon Strips 
  // some problems with the RO tree
  
  G4double TKRXStripX=0.;
  G4double TKRYStripY=0.;
  G4double TKRYStripX=0.; 
  G4double TKRXStripY=0.;
  
  TKRXStripX = TKRYStripY = GammaRayTelDetector->GetTKRSiliconPitch();
  TKRYStripX = TKRXStripY= GammaRayTelDetector->GetTKRActiveTileXY();
  G4double TKRZStrip  = GammaRayTelDetector->GetTKRSiliconThickness();
  
  G4int NbOfTKRStrips  = GammaRayTelDetector->GetNbOfTKRStrips();
  
  
  G4VSolid* solidTKRStripX = new G4Box("Strip X",			
				       TKRXStripX/2,TKRYStripX/2,
				       TKRZStrip/2); 
  
  G4LogicalVolume* logicTKRStripX = 
    new G4LogicalVolume(solidTKRStripX,dummyMat,"Strip X",0,0,0);	 
  
		
  G4VSolid* solidTKRStripY = new G4Box("Strip Y",			
				       TKRXStripY/2,TKRYStripY/2,
				       TKRZStrip/2); 
  

  G4LogicalVolume* logicTKRStripY = 
    new G4LogicalVolume(solidTKRStripY,dummyMat,"Strip Y",0,0,0);	 
							
						      
  G4VPhysicalVolume* physiTKRStripX = 0;
  G4VPhysicalVolume* physiTKRStripY = 0;
  G4double TKRSiliconPitch = GammaRayTelDetector->GetTKRSiliconPitch();

  for (i=0;i< NbOfTKRStrips; i++)
    {  
      physiTKRStripX = new 
	G4PVPlacement(0,G4ThreeVector(-TKRActiveTileXY/2 +TKRSiliconPitch/2 +
				      (i)*TKRSiliconPitch, 0., 0.),
		      "Strip X",		
		      logicTKRStripX,
		      physiTKRActiveTileXRO,
		      false,	
		      i);	

	
      physiTKRStripY = new 
	G4PVPlacement(0,G4ThreeVector(0.,-TKRActiveTileXY/2 
				      +TKRSiliconPitch/2 +
				      (i)*TKRSiliconPitch, 0.),
		      "Strip Y",		
		      logicTKRStripY,
		      physiTKRActiveTileYRO,
		      false,	
		      i);	
      
      



    }
  


  //Flags the strip as sensitive .The pointer here serves
  // as a flag only to check for sensitivity.
  // (Could we make it by a simple cast of a non-NULL value ?)
  
  
  GammaRayTelDummySD * dummySensi = new GammaRayTelDummySD;
  
  logicTKRStripX->SetSensitiveDetector(dummySensi);
  logicTKRStripY->SetSensitiveDetector(dummySensi);
  
  //logicTKRActiveTileXRO->SetSensitiveDetector(dummySensi);
  //logicTKRActiveTileYRO->SetSensitiveDetector(dummySensi);
    //logicTKRDetectorRO->SetSensitiveDetector(dummySensi);
  
  return ROWorldPhys;
}










