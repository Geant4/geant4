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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyDetectorROGeometry.hh"
#include "IAEADetectorMessenger.hh"
#include "HadrontherapyDetectorSD.hh"
#include "IAEADetectorConstruction.hh"
#include "HadrontherapyModulator.hh"

/////////////////////////////////////////////////////////////////////////////
IAEADetectorConstruction::IAEADetectorConstruction()
  : detectorSD(0), detectorROGeometry(0), 
    passiveProtonBeamLine(0), modulator(0),
    physicalTreatmentRoom(0),
    phantomPhysicalVolume(0), 
    detectorLogicalVolume(0), 
    detectorPhysicalVolume(0)
{
  // Messenger to change parameters of the geometry
  detectorMessenger = new IAEADetectorMessenger(this);

  // Detector sizes
  detectorSizeX = 20.*cm;
  detectorSizeY = 20.*cm;
  detectorSizeZ = 20.*cm;

  // Number of the detector voxels
  // This does notthing, since a primitive score is used  
  numberOfVoxelsAlongX = 400;
  numberOfVoxelsAlongY = 1;
  numberOfVoxelsAlongZ = 1;
  
  startDetectorThickness = 5.*cm; // approximation, exact value not given by Haettner
  phantomCenter = startDetectorThickness + 64.*cm;
  phantomDepth = 27.9 *cm;
  plexiThickness = 0.2 *cm;
  aluWindowThickness = 0.01 *cm;
  endDetectorThickness = 3.7 *cm;
  moveEndDetectorForward = 0 *cm; //if detector needs to be moved
  endDetectorPosition =  startDetectorThickness + 358 *cm + (endDetectorThickness / 2) - moveEndDetectorForward;
  
  noPhantom = false;
  
 }

/////////////////////////////////////////////////////////////////////////////
IAEADetectorConstruction::~IAEADetectorConstruction()
{ 
  if (detectorROGeometry) delete detectorROGeometry;  
  delete detectorMessenger;
}

/////////////////////////////////////////////////////////////////////////////
//void IAEADetectorConstruction::ChangeTheBeamLine(const G4String& name)
//if (name == emName) return;
//if (name == "ProtonPassiveBeamLine")
//{
//ConstructPassiveProtonBeamLine();
//}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* IAEADetectorConstruction::Construct()
{ 

  ConstructPassiveProtonBeamLine();
  ConstructDetector();
#ifdef G4ANALYSIS_USE_ROOT
  //write the metadata for analysis
  HadrontherapyAnalysisManager::GetInstance()->setGeometryMetaData((endDetectorPosition - endDetectorThickness/2)/10, phantomDepth/10, phantomCenter/10); //FIXME! unit correction hardcoded
#endif
  // Set the sensitive detector where the energy deposit is collected
  ConstructSensitiveDetector();
  return physicalTreatmentRoom;
}

/////////////////////////////////////////////////////////////////////////////
void IAEADetectorConstruction::ConstructPassiveProtonBeamLine()
{ 
  // -----------------------------
  // Treatment room - World volume
  //------------------------------
  // Treatment room sizes
  const G4double worldX = 400.0 *cm;
  const G4double worldY = 200.0 *cm; //to fit, new bigger detector
  const G4double worldZ = 200.0 *cm;
  G4bool isotopes = false;
 
  G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  G4Box* treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);
  G4LogicalVolume* logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, 
                                                            airNist, 
                                                            "logicTreatmentRoom", 
															0,0,0);
  physicalTreatmentRoom = new G4PVPlacement(0,
					    G4ThreeVector(),
					    "physicalTreatmentRoom", 
					    logicTreatmentRoom, 
					    0,false,0);
 

  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::Invisible);
 

 
  //----------------------------------------
  // Phantom:
  // A box used to approximate tissues. Is surrounded by plexi-glas.
  //----------------------------------------

  G4Material* waterNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", isotopes);
  G4Material* plexiGlas = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
  //G4Box* phantom = new G4Box("Phantom",10 *cm, 20 *cm, 20 *cm);
  //the below for integrated angular distribution plot
  G4Box* phantom = new G4Box("Phantom",phantomDepth/2, 32 *cm, 32 *cm); //< E.haettner has 32x8, but simulation takes scattering symmetrically
  G4Box* plexiSheet = new G4Box("phantomEdge",plexiThickness/2, 32 *cm, 32 *cm);
  G4LogicalVolume* phantomLogicalVolume = new G4LogicalVolume(phantom,	
							      waterNist, 
							      "phantomLog", 0, 0, 0);
								  
  G4LogicalVolume* phantomEdgeLogicalVolume = new G4LogicalVolume(plexiSheet,	
							      plexiGlas, 
							      "phantomEdgeLog", 0, 0, 0);
					//5.5cm for veto and start detector, see fig 5.1 and 4.1 ref Haettner 2006 
  if(!noPhantom){  
		  phantomPhysicalVolume = new G4PVPlacement(0,G4ThreeVector(phantomCenter, 0.*cm, 0.*cm),
								"phantomPhys",
								phantomLogicalVolume,
								physicalTreatmentRoom,
								false,0);

		  phantomEdge1PhysicalVolume = new G4PVPlacement(0,G4ThreeVector(phantomCenter - phantomDepth/2 - plexiThickness/2, 0.*cm, 0.*cm),
								"phantomEdgePhys",
								phantomEdgeLogicalVolume,
								physicalTreatmentRoom,
								false,0);
		  phantomEdge2PhysicalVolume = new G4PVPlacement(0,G4ThreeVector(phantomCenter + phantomDepth/2 + plexiThickness/2, 0.*cm, 0.*cm),
								"phantomEdgePhys",
								phantomEdgeLogicalVolume,
								physicalTreatmentRoom,
								false,0);
	}
  //----------------------------------------
  // Beamwindow:
  // The aluminium-window of the beam-source
  //----------------------------------------
  G4Material* aluNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
  G4Box* beamWindow = new G4Box("beamwindow",aluWindowThickness/2, 10 *cm, 10 *mm);
  G4LogicalVolume* beamWindowLogicalVolume = new G4LogicalVolume(beamWindow,	
							      aluNist, 
							      "beamWindowLog", 0, 0, 0);
  beamWindowPhysicalVolume = new G4PVPlacement(0,G4ThreeVector(1.*mm, 0.*mm, 0.*mm),
					    "beamPhys",
					    beamWindowLogicalVolume,
					    physicalTreatmentRoom,
					    false,0); //just ahead of beam that starts at origo
  //----------------------------------------
  // NewDetector:
  // A box used to simulate the end detector
  //----------------------------------------
  //G4Material* NewDetectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  //G4Box* NewDetector = new G4Box("NewDetector",1 *cm, 20 *cm, 27.9 *cm);
  
  // Visualisation attributes of the patient
  red = new G4VisAttributes(G4Colour(255/255., 0/255. ,0/255.));
  red -> SetVisibility(true);
  red -> SetForceSolid(true);
  //red -> SetForceWireframe(true);
  phantomLogicalVolume -> SetVisAttributes(red);
}

void IAEADetectorConstruction::setWaterThickness(G4double newWaterThickness){
	//This has to be run before the elements are constructed
	//does not support multiple thicknesses per run.
	if(newWaterThickness > 0){
	this->phantomDepth = newWaterThickness;
	}else{
	this->noPhantom = true;
	}
#ifdef G4ANALYSIS_USE_ROOT
	  //update the geometry metadata
	  HadrontherapyAnalysisManager::GetInstance()->setGeometryMetaData((this->endDetectorPosition - this->endDetectorThickness/2)/10, this->phantomDepth/10, this->phantomCenter/10); //FIXME! unit correction hardcoded
#endif
	}


/////////////////////////////////////////////////////////////////////////////
void IAEADetectorConstruction::ConstructDetector()
{
  //-----------
  // Braggcurve Detector
  //-----------
  //Is currently not used
  
  G4bool isotopes =  false; 
  G4Material* waterNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", isotopes);
  G4Box* detector = new G4Box("Detector",detectorSizeX,detectorSizeY,detectorSizeZ);

  detectorLogicalVolume = new G4LogicalVolume(detector,
					      waterNist,
					      "DetectorLog",
					      0,0,0);
  
  //G4double detectorXtranslation = -180.*mm;
  G4double detectorXtranslation = -180.*cm;
  detectorPhysicalVolume = new G4PVPlacement(0,
					     G4ThreeVector(detectorXtranslation, 0.0 *mm, 0.0 *mm),
					     "DetectorPhys",
					     detectorLogicalVolume,
					     phantomPhysicalVolume,
					     false,0);
  
  // Visualisation attributes of the phantom
  skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
  skyBlue -> SetVisibility(true);
  skyBlue -> SetForceSolid(true);

  detectorLogicalVolume -> SetVisAttributes(skyBlue);
  
  //Visualization attributes for the beamwindow
  
  //-----------
  // NewDetector (mwpc etc. type behind the phantom)
  //-----------
  G4Material* NewDetectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", false);
  G4Box* NewDetector = new G4Box("NewDetector",endDetectorThickness/2,190.*cm,190.*cm); //huge detector, will be scaled in root
  //For integrated angular distribution below
  NewDetectorLogicalVolume = new G4LogicalVolume(NewDetector,
					      NewDetectorMaterial,
					      "NewDetectorLog",
					      0,0,0);
  NewDetectorPhysicalVolume = new G4PVPlacement(0,
					     G4ThreeVector(endDetectorPosition, 0.0 *cm, 0.0 *cm),
					     "NewDetectorPhys",
					     NewDetectorLogicalVolume,
					     physicalTreatmentRoom,
					     false,0);
  
}

/////////////////////////////////////////////////////////////////////////////
void  IAEADetectorConstruction::ConstructSensitiveDetector()
{  
  // Sensitive Detector and ReadOut geometry definition
  G4SDManager* sensitiveDetectorManager = G4SDManager::GetSDMpointer();

  G4String sensitiveDetectorName = "Detector"; 
  if(!detectorSD)
    {
      // The sensitive detector is instantiated
      detectorSD = new HadrontherapyDetectorSD(sensitiveDetectorName);
      // The Read Out Geometry is instantiated
      G4String ROGeometryName = "DetectorROGeometry";
      detectorROGeometry = new HadrontherapyDetectorROGeometry(ROGeometryName,
							     G4ThreeVector(20. *mm, // detector to world position!
								               0. *mm, 
								               0. *mm),
							     phantomDepth/2,
							     detectorSizeY,
							     detectorSizeZ,
							     numberOfVoxelsAlongX,
							     numberOfVoxelsAlongY,
							     numberOfVoxelsAlongZ);
      detectorROGeometry -> BuildROGeometry();
      detectorSD -> SetROgeometry(detectorROGeometry);
      sensitiveDetectorManager -> AddNewDetector(detectorSD);
      detectorLogicalVolume -> SetSensitiveDetector(detectorSD);
    }
}
