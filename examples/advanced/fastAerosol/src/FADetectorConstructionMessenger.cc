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

// (adapted from B2aDetectorMessenger)
// A.Knaian, N.MacFadden

#include "FADetectorConstructionMessenger.hh"
#include "FADetectorConstruction.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* detectorIn)
 : G4UImessenger()
{
	detector = detectorIn;

	// Directory
	//
	// /geometry
	fGeometryDirectory = new G4UIdirectory("/geometry/");
	fGeometryDirectory->SetGuidance("Geometry setup.");

	// Physics
	//
	// /geometry/stepLim
	fStepLimCmd = new G4UIcmdWithADoubleAndUnit("/geometry/stepLim",this);
	fStepLimCmd->SetGuidance("Maximum step length.");
	fStepLimCmd->SetParameterName("stepLim",false);
	fStepLimCmd->SetRange("stepLim>=0.");
	fStepLimCmd->SetDefaultValue(DBL_MAX);
	fStepLimCmd->SetDefaultUnit("mm");
	fStepLimCmd->AvailableForStates(G4State_PreInit);

	// Cloud droplet settings
	// /geometry/dropletR
	fDropletRCmd = new G4UIcmdWithADoubleAndUnit("/geometry/dropletR",this);
	fDropletRCmd->SetGuidance("Minimal bounding radius of droplet.");
	fDropletRCmd->SetParameterName("dropletR",false);
	fDropletRCmd->SetRange("dropletR>0.");
	fDropletRCmd->SetDefaultValue(1.0);
	fDropletRCmd->SetDefaultUnit("mm");
	fDropletRCmd->AvailableForStates(G4State_PreInit);

	// /geometry/dropletNumDens 		
	fDropletNumDensCmd = new G4UIcmdWithADouble("/geometry/dropletNumDens", this);
	fDropletNumDensCmd->SetGuidance("Number of droplets per mm^3.");				// would be nice to have official number density units
	fDropletNumDensCmd->SetParameterName("dropletCOunt",false);
	fDropletNumDensCmd->SetDefaultValue(0);
	fDropletNumDensCmd->AvailableForStates(G4State_PreInit);

	// Cloud build type
	//
	// /geometry/fastAerosol
	fFastAerosolCloudCmd = new G4UIcmdWithABool("/geometry/fastAerosolCloud",this);
	fFastAerosolCloudCmd->SetGuidance("Whether or not to build the fastAerosol cloud.");
	fFastAerosolCloudCmd->SetParameterName("fastAerosol",false);
	fFastAerosolCloudCmd->SetDefaultValue(false);
	fFastAerosolCloudCmd->AvailableForStates(G4State_PreInit);

	// /geometry/parameterisedCloud
	fParameterisedCloudCmd = new G4UIcmdWithABool("/geometry/parameterisedCloud",this);
	fParameterisedCloudCmd->SetGuidance("Whether or not to build the parameterised cloud.");
	fParameterisedCloudCmd->SetParameterName("parameterisedCloud",false);
	fParameterisedCloudCmd->SetDefaultValue(false);
	fParameterisedCloudCmd->AvailableForStates(G4State_PreInit);

	// /geometry/smoothCloud
	fSmoothCloudCmd = new G4UIcmdWithABool("/geometry/smoothCloud",this);
	fSmoothCloudCmd->SetGuidance("Whether or not to build the smooth cloud.");
	fSmoothCloudCmd->SetParameterName("parameterisedCloud",false);
	fSmoothCloudCmd->SetDefaultValue(false);
	fSmoothCloudCmd->AvailableForStates(G4State_PreInit);

	// fastAerosol cloud details
	//
	// /geometry/cloudShape
	fCloudShapeCmd = new G4UIcmdWithAString("/geometry/cloudShape",this);
	fCloudShapeCmd->SetGuidance("Cloud bulk shape");
	fCloudShapeCmd->SetParameterName("cloudShapeStr",false);
	fCloudShapeCmd->AvailableForStates(G4State_PreInit);

	// /geometry/dropletShape
	fDropletShapeCmd = new G4UIcmdWithAString("/geometry/dropletShape",this);
	fDropletShapeCmd->SetGuidance("Cloud droplet shape");
	fDropletShapeCmd->SetParameterName("dropletShapeStr",false);
	fDropletShapeCmd->AvailableForStates(G4State_PreInit);

	// /geometry/prePopulate
	fPrePopulateCmd = new G4UIcmdWithABool("/geometry/prePopulate",this);
	fPrePopulateCmd->SetGuidance("Whether or not to populate the cloud at the beginning.");
	fPrePopulateCmd->SetParameterName("prePopulate",false);
	fPrePopulateCmd->SetDefaultValue(false);
	fPrePopulateCmd->AvailableForStates(G4State_PreInit);

	// /geometry/minSpacing
	fMinSpacingCmd = new G4UIcmdWithADoubleAndUnit("/geometry/minSpacing",this);
	fMinSpacingCmd->SetGuidance("Minimum spacing between surfaces of spheres when generating random cloud of spheres.");
	fMinSpacingCmd->SetParameterName("minSpacing",false);
	fMinSpacingCmd->SetRange("minSpacing>0.");
	fMinSpacingCmd->SetDefaultValue(10.);
	fMinSpacingCmd->SetDefaultUnit("micrometer");
	fMinSpacingCmd->AvailableForStates(G4State_PreInit);

	/*
	// /geometry/gridPitch
	fGridPitchCmd = new G4UIcmdWithADoubleAndUnit("/geometry/gridPitch",this);
	fGridPitchCmd->SetGuidance("Pitch of the collision detection grid, must be larger than the sphere size.");
	fGridPitchCmd->SetParameterName("gridPitch",false);
	fGridPitchCmd->SetRange("gridPitch>0.");
	fGridPitchCmd->SetDefaultValue(0.5);
	fGridPitchCmd->SetDefaultUnit("mm");
	fGridPitchCmd->AvailableForStates(G4State_PreInit);
	*/

	// /geometry/setSmartless
	fSmartlessCmd = new G4UIcmdWithADouble("/geometry/smartless", this);
	fSmartlessCmd->SetGuidance("Set the 'smartless' parameter for the parameterised cloud.");
	fSmartlessCmd->SetParameterName("smartless",false);
	fSmartlessCmd->SetRange("smartless>0.");
	fSmartlessCmd->SetDefaultValue(2.0);
	fSmartlessCmd->AvailableForStates(G4State_PreInit);

	// /geometry/cloudSeed
	fCloudSeedCmd = new G4UIcmdWithAnInteger("/geometry/cloudSeed", this);
	fCloudSeedCmd->SetGuidance("Base of the random seed for the cloud sphere positions.");
	fCloudSeedCmd->SetParameterName("cloudSeed",false);
	fCloudSeedCmd->SetDefaultValue(0);
	fCloudSeedCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{
	delete fGeometryDirectory;

	delete fStepLimCmd;

	delete fDropletRCmd;
	delete fDropletNumDensCmd;

	delete fFastAerosolCloudCmd;
	delete fParameterisedCloudCmd;
	delete fSmoothCloudCmd;

	delete fCloudShapeCmd;
	delete fDropletShapeCmd;
	delete fPrePopulateCmd;
	delete fMinSpacingCmd;
	//delete fGridPitchCmd;

	delete fSmartlessCmd;

	delete fCloudSeedCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstructionMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{
	// Geometry Commands

	if( command == fFastAerosolCloudCmd ) {
		detector->fastAerosolCloud = (fFastAerosolCloudCmd->GetNewBoolValue(newValue));
	}

	if( command == fStepLimCmd ) {
		detector->stepLim = (fStepLimCmd->GetNewDoubleValue(newValue));
	}

	if( command == fDropletRCmd ) {
		detector->dropletR = (fDropletRCmd->GetNewDoubleValue(newValue));
	}
	if( command == fDropletNumDensCmd ) {
		detector->dropletNumDens = (fDropletNumDensCmd->GetNewDoubleValue(newValue));
	}

	if( command == fParameterisedCloudCmd ) {
		detector->parameterisedCloud = (fParameterisedCloudCmd->GetNewBoolValue(newValue));
	}
	if( command == fSmoothCloudCmd ) {
		detector->smoothCloud = (fSmoothCloudCmd->GetNewBoolValue(newValue));
	}
	if( command == fPrePopulateCmd ) {
		detector->prePopulate = (fPrePopulateCmd->GetNewBoolValue(newValue));
	}

	if( command == fCloudShapeCmd ) {
		detector->cloudShapeStr = newValue;
	}
	if( command == fDropletShapeCmd ) {
		detector->dropletShapeStr = newValue;
	}
	if( command == fMinSpacingCmd ) {
		detector->minSpacing = (fMinSpacingCmd->GetNewDoubleValue(newValue));
	}
	//if( command == fGridPitchCmd ) {
	//	detector->gridPitch = (fGridPitchCmd->GetNewDoubleValue(newValue));
	//}
	if( command == fSmartlessCmd ) {
		detector->smartless = (fSmartlessCmd->GetNewDoubleValue(newValue));
	}
	if( command == fCloudSeedCmd ) {
		detector->cloudSeed = (fCloudSeedCmd->GetNewIntValue(newValue));
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
