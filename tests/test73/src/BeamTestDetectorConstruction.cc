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
// $Id:$
// GEANT4 tag $Name:$
//
// BeamTestDetector Construction class used to make the geometry setup of the simulation
// includes the parametrisation class defined in BeamTestCellPatameterisation 
//

#include "BeamTestDetectorConstruction.hh"
#include "BeamTestSiliconMonitor.hh"
#include "BeamTestCellParameterisation.hh"

#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"     
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
// Defining sensitive detector
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4StepLimiter.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
// Scoring components
#include "G4MultiFunctionalDetector.hh"
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4SDManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

// Filter
#include "G4SDParticleWithEnergyFilter.hh"
#include "BeamTestDetectorMessenger.hh"

////////////////////////////////////////////////////////////////////////

BeamTestDetectorConstruction::BeamTestDetectorConstruction(/*Parameters* parameter*/)
:solidWorld(0),  logicWorld(0),  physiWorld(0),
	solidTracker(0),logicTracker(0),physiTracker(0), 
	solidChamber(0),logicChamber(0),physChamber(0),physiChamber(0), 
	ChamberMater(0),chamberParam(0),stepLimit(0),
	fWorldLength(0.),  fTrackerLength(0.),
	NbOfChambers(20) ,  ChamberWidth(1*mm),  ChamberSpacing(1*um), normalise(0.),
    messenger(new BeamTestDetectorMessenger(this)),
    monitor(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamTestDetectorConstruction::~BeamTestDetectorConstruction()
{
	delete chamberParam;
    delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* BeamTestDetectorConstruction::Construct()
{
    //This method can be called several times, remeber if I've called several times
    static G4bool onlyonce = true;
    //Open and clean old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    
    //****************************************************************************************************
    //Materials defitnion
    //****************************************************************************************************
    
    G4NistManager* matman = G4NistManager::Instance();
    G4Material* vacuum = matman->FindOrBuildMaterial("G4_Galactic");
    G4Material* Al = matman->FindOrBuildMaterial("G4_Al");

    if ( onlyonce ) {
        G4cout<<"Defined materials:"<<G4endl;
        G4cout<<*(G4Material::GetMaterialTable())<<G4endl;
    }
//********************************************************************************************************
	//--------- Sizes of the principal geometrical components (solids)  ---------
//********************************************************************************************************
	fTrackerLength = 1.0*m; // Full length of Tracker

	ChamberMater = Al;

	fWorldLength= 2.0*m;

	G4double trackerSize = 0.5*fTrackerLength;   // Half length of the Tracker


//********************************************************************************************************
	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
//********************************************************************************************************

	//------------------------------ 
	// World
	//------------------------------ 
    G4double HalfWorldLength = 0.5*fWorldLength;
    if ( onlyonce ) {
        G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
    }

	G4cout << "Computed tolerance = "
		<< G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
		<< " mm" << G4endl;
    solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
    logicWorld= new G4LogicalVolume( solidWorld, vacuum, "World", 0, 0, 0);
        
    //  Must place the World Physical volume unrotated at (0,0,0).
    // 
    physiWorld = new G4PVPlacement(0,               // no rotation
                                   G4ThreeVector(), // at (0,0,0)
                                   logicWorld,      // its logical volume
                                   "World",         // its name
                                   0,               // its mother  volume
                                   false,           // no boolean operations
                                   0);              // copy number
	
	//------------------------------ 
	// Tracker
	//------------------------------
	G4ThreeVector positionTracker = G4ThreeVector(0,0,0);

    solidTracker = new G4Box("tracker",trackerSize,trackerSize,trackerSize);
    logicTracker = new G4LogicalVolume(solidTracker ,vacuum, "Tracker",0,0,0);  
    physiTracker = new G4PVPlacement(0,              // no rotation
                                     positionTracker, // at (x,y,z)
                                     logicTracker,    // its logical volume				  
                                     "Tracker",       // its name
                                     logicWorld,      // its mother  volume
                                     false,           // no boolean operations
                                     0);              // copy number 
	//------------------------------ 
	// Tracker segments
	//------------------------------
	// 
	// An example of Parameterised volumes
	// dummy values for G4Box -- modified by parameterised volume
    G4double firstPosition = 0.5*ChamberWidth;
    solidChamber = new G4Box("chamber", 25*cm, 25*cm, ChamberWidth*0.5); 
    logicChamber = new G4LogicalVolume(solidChamber,ChamberMater,"Chamber",0,0,0);

    chamberParam = new BeamTestCellParameterisation(  
                                                    NbOfChambers,          // NoChambers 
                                                    firstPosition,         // Z of center of first 
                                                    ChamberSpacing,        // Z spacing of centers
                                                    ChamberWidth);          // Width Chamber
    physiChamber = new G4PVParameterised(
                                         "Chamber",       // their name
                                         logicChamber,    // their logical volume
                                         logicTracker,    // Mother logical volume
                                         kZAxis,          // Are placed along this axis 
                                         NbOfChambers,    // Number of chambers
                                         chamberParam);   // The parametrisation
	G4cout << "There are " << NbOfChambers << " chambers in the tracker region. "
		<< "The chambers are " << ChamberWidth/mm << " mm of " 
		<< ChamberMater->GetName() << "\n The distance between chamber is "
		<< ChamberSpacing/um << " um" << G4endl;



	//------------------------------------------------ 
	// Sensitive detectors
	//------------------------------------------------ 
    G4String trackerChamberSDname = "/BeamTest/TrackerChamberSD";
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if ( monitor == 0 ) {
      monitor = new BeamTestSiliconMonitor(trackerChamberSDname);
      SDman->AddNewDetector( monitor );
    }
    logicChamber->SetSensitiveDetector( monitor );
    //SDman->ListTree();
	
    //--------- Visualization attributes -------------------------------
    G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    logicWorld  ->SetVisAttributes(BoxVisAtt);  
    logicTracker->SetVisAttributes(BoxVisAtt);

    G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    logicChamber->SetVisAttributes(ChamberVisAtt);
	//--------- example of User Limits -------------------------------

	// below is an example of how to set tracking constraints in a given
	// logical volume(see also in N02PhysicsList how to setup the processes
	// G4StepLimiter or G4UserSpecialCuts).

	// Sets a max Step length in the tracker region, with G4StepLimiter
	//
	//G4double maxStep = ChamberWidth/normalise;
	//stepLimit = new G4UserLimits(50*um);
	//logicChamber->SetUserLimits(stepLimit);

	// Set additional contraints on the track, with G4UserSpecialCuts
	//
	// G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
	// logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
	//                                               minEkin));
    if ( onlyonce ) onlyonce = false;
	return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamTestDetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

