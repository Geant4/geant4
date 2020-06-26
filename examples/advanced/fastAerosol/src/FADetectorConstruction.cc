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

// (adapted from B1DetectorConstruction)
// Author: A.Knaian (ara@nklabs.com), N.MacFadden (natemacfadden@gmail.com)

#include "FADetectorConstruction.hh"
#include "FADetectorConstructionMessenger.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

// shapes
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"

// to build FastAerosol cloud
#include "FastAerosolSolid.hh"

// to build parameterised cloud
#include "FACloudParameterisation.hh"
#include "G4PVParameterised.hh"
#include <fstream>

// step limits
#include "G4UserLimits.hh"

// visualization
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// to save distribution
#include <sys/stat.h>
#include <ctime> 		// for measuring FastAerosol droplet center population time


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
fScoringVolume(0)
{ 
fMessenger = new DetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
	delete fMessenger;

	delete fStepLimits;

	delete fCloudShape;
	delete fDropletShape;

	delete fCloud;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	//
	// Check cloud build settings
	//
	if (fFastAerosolCloud + fParameterisedCloud + fSmoothCloud > 1)
	{
		std::ostringstream message;
		message << "Must select at most one build type! Selections:" << G4endl
				<< "     fFastAerosolCloud = " << fFastAerosolCloud << G4endl
				<< "     fParameterisedCloud = " << fParameterisedCloud << G4endl
				<< "     fSmoothCloud = " << fSmoothCloud << G4endl;
		G4Exception("DetectorConstruction::Construct()", "GeomSolids0002",
					FatalException, message);
	}


	//
	// Get nist material manager
	//
	G4NistManager* nist = G4NistManager::Instance();


	//
	// Option to switch on/off checking of volumes overlaps
	//
	G4bool checkOverlaps = false;


	//
	// Large scale geometry dimensions
	//
	G4double cloud_sizeXY = 0.5*m;
	G4double cloud_sizeZ = 5.0*m;

	G4double world_sizeXY = 1.1*(cloud_sizeXY);
	G4double world_sizeZ= 1.1*(cloud_sizeZ);


	//
	// Cloud shape
	//
	if (fCloudShapeStr == "box")
	{
		G4cout << "Cloud shape = box" << G4endl;
		fCloudShape = new G4Box("cloudShape", 0.5*cloud_sizeXY, 0.5*cloud_sizeXY, 0.5*cloud_sizeZ);
	}
	else if (fCloudShapeStr == "ellipsoid")
	{
		G4cout << "Cloud shape = ellipsoid" << G4endl;
		fCloudShape = new G4Ellipsoid("cloudShape", 0.5*cloud_sizeXY, 0.5*cloud_sizeXY, 0.5*cloud_sizeZ, 0, 0);
	}
	else if (fCloudShapeStr == "cylinder")
	{
		G4cout << "Cloud shape = cylinder" << G4endl;
		fCloudShape = new G4Tubs("cloudShape", 0.0, 0.5*cloud_sizeXY, 0.5*cloud_sizeZ, 0, 360*deg);
	}
	else if (fCloudShapeStr == "pipe")
	{
		G4cout << "Cloud shape = pipe" << G4endl;
		fCloudShape = new G4Tubs("cloudShape", 0.25*cloud_sizeXY, 0.5*cloud_sizeXY, 0.5*cloud_sizeZ, 0, 360.*deg);
	}
	else
	{
		std::ostringstream message;
		message << "Invalid cloud shape = " << fCloudShapeStr << "!";
		G4Exception("DetectorConstruction::Construct()", "GeomSolids0002",
					FatalException, message);
	}

	//
	// Droplet Shape
	//

	// The difference in radii of the maximal sphere (centered at the origin) contained in the droplet and the minimal sphere (centered at the origin) containing the droplet
	G4double sphericalUncertainty = 0.0;

	if (fDropletShapeStr == "sphere")
	{
		G4cout << "Droplet shape = sphere" << G4endl;
		fDropletShape = new G4Orb("dropletSV", fDropletR);
		sphericalUncertainty = 0.0;
	}
	else if (fDropletShapeStr == "halfSphere")
	{
		G4cout << "Droplet shape = halfSphere" << G4endl;
		fDropletShape = new G4Sphere("dropletSV", 0.0, fDropletR,
										0.0, 180.*deg,
										0.0, 180.*deg);
		sphericalUncertainty = fDropletR;
	}
	else if (fDropletShapeStr == "cylinder")
	{
		G4cout << "Droplet shape = cylinder" << G4endl;
		fDropletShape = new G4Tubs("dropletSV", 0, fDropletR/std::sqrt(3), fDropletR/std::sqrt(3), 0, 360.*deg);
		sphericalUncertainty = fDropletR*(1-1/std::sqrt(3));
	}
	else if (fDropletShapeStr == "box")
	{
		G4cout << "Droplet shape = box" << G4endl;
		fDropletShape = new G4Box("dropletSV", fDropletR/std::sqrt(3), fDropletR/std::sqrt(3), fDropletR/std::sqrt(3));
		sphericalUncertainty = fDropletR*(1-1/std::sqrt(3));
	}
	else
	{
		std::ostringstream message;
		message << "Invalid droplet shape = " << fCloudShapeStr << "!";
		G4Exception("DetectorConstruction::Construct()", "GeomSolids0002",
					FatalException, message);
	}


	//
	// Materials
	//

	// Compute the density of air at 14 km using the Barometric formula
	// see, e.g., https://en.wikipedia.org/wiki/Density_of_air
	G4double h = 14.0*km;
	G4double p0 = 101325*hep_pascal;
	G4double T0 = 288.15*kelvin;
	G4double grav = 9.80665*m/(s*s);
	G4double La = 0.0065*kelvin/m;
	G4double R = 8.31447*joule/(mole*kelvin);
	G4double M = 0.0289644*kg/mole;
	G4double T = T0 - La*h;
	G4double p = p0*std::pow(1-La*h/T0,grav*M/(R*La));
	G4double air_density = p*M/(R*T);

	// make materials and set densities
	G4Material* air_mat = nist->BuildMaterialWithNewDensity("Atmosphere","G4_AIR",air_density);
	G4Material* water_mat = nist->FindOrBuildMaterial("G4_WATER");

	G4double water_density = water_mat->GetDensity();
	G4double ice_density = 0.9168*g/cm3;
	
	G4Material* ice_mat = new G4Material("Water ice ", ice_density, 1, kStateSolid, T, p);
	ice_mat->AddMaterial(water_mat, 1.);


	//
	// Droplets
	//
	G4double droplet_density = water_density;
	G4Material* droplet_mat = water_mat;

	G4double droplet_count = fDropletNumDens*(fCloudShape->GetCubicVolume());

	G4double droplet_volume = fDropletShape->GetCubicVolume();
	G4double droplet_total_volume = droplet_count*droplet_volume;

	G4double droplet_total_mass =  droplet_total_volume*droplet_density;


	//
	// Cloud macroscopic quantities
	//
	G4double cloud_volume = fCloudShape->GetCubicVolume();
	G4double cloud_air_volume = cloud_volume - droplet_total_volume;
	G4double cloud_air_mass = air_density*cloud_air_volume;


	//
	// Step limit
	//
	fStepLimits = new G4UserLimits(fStepLim);


	//
	// Build world
	//
	G4Box* solidWorld =	
		new G4Box("World",					//its name
				  0.5*world_sizeXY,			//half x-span
				  0.5*world_sizeXY,			//half y-span
				  0.5*world_sizeZ);			//half z-span

	G4LogicalVolume* logicWorld =						 
		new G4LogicalVolume(solidWorld,		//its solid
							air_mat,		//its material
							"World");		//its name

	logicWorld->SetUserLimits(fStepLimits);
									 
	G4VPhysicalVolume* physWorld = 
		new G4PVPlacement(0,				//no rotation
						  G4ThreeVector(),	//at (0,0,0)
						  logicWorld,		//its logical volume
						  "World",			//its name
						  0,				//its mothervolume
						  false,			//no boolean operation
						  0,				//copy number
						  checkOverlaps);	//overlaps checking


	//
	// Build cloud
	//
	G4LogicalVolume* logicCloud;

	// **********************************************************
	// 
	// Build the cloud using the FastAerosol geometry class
	// 
	// ***********************************************************
	if (fFastAerosolCloud) {
		G4cout << "\nFastAerosol geometry with n=" << fDropletNumDens*mm3 << "/mm3, r=" << fDropletR/mm << "mm spheres.\n" << G4endl;
		
		fCloud = new FastAerosol("cloud",
								   fCloudShape,					//cloud shape
								   fDropletR,					//bounding radius of droplets
								   fMinSpacing,					//minimum spacing between droplets
								   fDropletNumDens,				//approximate number of droplets in cloud
								   sphericalUncertainty);		//uncertainty in distance to droplet surface from outside using just droplet's origin as info
		fCloud->SetDropletsPerVoxel(4);
		
		/*
		fCloud = new FastAerosol("fCloud",
								   fCloudShape,					//cloud shape
								   fDropletR,					//bounding radius of droplets
								   fMinSpacing,					//minimum spacing between droplets
								   fDropletNumDens,				//approximate number of droplets in cloud
								   sphericalUncertainty,		//uncertainty in distance to droplet surface from outside using just droplet's origin as info
								   [](G4ThreeVector pos) {return pos.x();});					//number density distribution function
		*/

		FastAerosolSolid* solidCloud =
			new FastAerosolSolid("cloudSV",						//its name
									fCloud,						//its shape
									fDropletShape);				//its droplets

		/*
		FastAerosolSolid* solidCloud =
			new FastAerosolSolid("cloudSV",						//its name
									fCloud,						//its shape
									fDropletShape,				//its droplets
									[](G4ThreeVector) {G4RotationMatrix rotm = G4RotationMatrix(); rotm.rotateY(90.0*deg); return rotm;});	//droplet rotation function
		*/

		solidCloud->SetStepLim(fStepLim);						//FastAerosol can use step limit to speed calculations

		logicCloud =
			new G4LogicalVolume(solidCloud,						//its solid
								droplet_mat,					//its material
								"cloudLV");						//its name
		logicCloud->SetUserLimits(fStepLimits);
		logicCloud->SetVisAttributes(G4VisAttributes(G4Colour(0.0,0.0,1.0,0.4)));

		new G4PVPlacement(0,									//no rotation
						  G4ThreeVector(),						//at position
						  logicCloud,							//its logical volume
						  "cloudPV",							//its name
						  logicWorld,							//its mother volume
						  false,								//no boolean operation
						  0,									//copy number 
						  checkOverlaps);						//overlaps checking


		fCloud->SetSeed(fCloudSeed);

		// fPrePopulate = whether to populate all voxels at the beginning or on the fly
		if (fPrePopulate) {
			// populate (proving it to the user by printing population reports)
			clock_t t;
			t = clock();

			G4cout << "\nBefore populating" << G4endl;
			G4cout <<   "=================" << G4endl;
			fCloud->PrintPopulationReport();
			G4cout << "\nPopulating..." << G4endl;
			fCloud->PopulateAllGrids();
			G4cout << "\nAfter populating" << G4endl;
			G4cout <<   "================" << G4endl;
			fCloud->PrintPopulationReport();
			G4cout << G4endl;

			t = clock() - t;

			G4cout << "\nThis took " << ((float)t)/CLOCKS_PER_SEC << "s\n" << G4endl;

			// make filename variables to save data
			G4String rStr = std::to_string(fDropletR/mm);
			rStr.erase ( rStr.find_last_not_of('0') + 1, std::string::npos );	// drop trailing 0
			replace( rStr.begin(), rStr.end(), '.', 'p');
			if (rStr.back() == 'p') { rStr.pop_back(); }	// don't write "3p" for 3.0, just write "3"

			// want to represent the number density as 1E-ApB for some A, B
			G4int order10 = (G4int) -round(10*std::log10(fDropletNumDens*mm3));	// gives 10x the exponent rounded to the int (10x so we get two decimals)
			G4int leading = order10 / 10;	// first number
			G4int trailing = order10 % 10;	// second number
			G4String nStr = "1E-" + std::to_string(leading) + "p" + std::to_string(trailing);		


			// save population time
			std::ofstream file;
			file.open("popTime_r" + rStr + "mm_n" + nStr + "mm-3.csv");
			file << ((float)t)/CLOCKS_PER_SEC;
			file.close();

			

			// save distribution
			G4String fName = "distribution_r" + rStr + "mm_n" + nStr + "mm-3.csv";
			fCloud->SaveToFile(fName);
		}
	}
	// **********************************************************
	// 
	// (For comparision/benchmarking) Build the cloud using G4VParameterized (does not use FastAerosol)
	// 
	// ***********************************************************

	// the droplet positions for this cloud are those saved in the "distribution" folder of our data
	// this is to make comparable simulations between FastAerosol and parameterised clouds
	// this requires that we first simulate FastAerosol (pre-populated) to generate the positions
	
	else if (fParameterisedCloud)
	{
		G4cout << "\nParameterised geometry with n=" << fDropletNumDens*mm3 << "/mm3 and r=" << fDropletR/mm << "mm spheres.\n" << G4endl;
		std::vector<G4ThreeVector> positions;
		G4double x,y,z;

		// load distribution file
		G4String fName;

		G4String rStr = std::to_string(fDropletR/mm);
		rStr.erase ( rStr.find_last_not_of('0') + 1, std::string::npos );	// drop trailing 0
		replace( rStr.begin(), rStr.end(), '.', 'p');
		if (rStr.back() == 'p') { rStr.pop_back(); }	// don't write "3p" for 3.0, just write "3"

		// want to represent the number density as 1E-ApB for some A, B
		G4int order10 = (G4int) -round(10*std::log10(fDropletNumDens*mm3));	// gives 10x the exponent rounded to the int (10x so we get two decimals)
		G4int leading = order10 / 10;	// first number
		G4int trailing = order10 % 10;	// second number
		G4String nStr = "1E-" + std::to_string(leading) + "p" + std::to_string(trailing);

		fName = "distribution_r" + rStr + "mm_n" + nStr + "mm-3.csv";

		std::ifstream infile(fName);
		std::string line;
		
		while (getline(infile,line)) {
			std::istringstream stream(line);
			std::string field;

			getline(stream,field,','); x = stod(field)*mm;
			getline(stream,field,','); y = stod(field)*mm;
			getline(stream,field,','); z = stod(field)*mm;

			positions.push_back(G4ThreeVector(x,y,z));
		}

		G4VPVParameterisation* cloudParam =
			new CloudParameterisation(positions);

		G4Box* cloudBounding = 
			new G4Box("cloudBounding",				//its name
				  0.5*cloud_sizeXY,					//half x-span
				  0.5*cloud_sizeXY,					//half y-span
				  0.5*cloud_sizeZ);					//half z-span
							
		logicCloud =						 
			new G4LogicalVolume(cloudBounding,		//its solid
								air_mat,			//its material
								"cloudLV");			//its name

		logicCloud->SetSmartless(fSmartless);
		logicCloud->SetUserLimits(fStepLimits);
		logicCloud->SetVisAttributes(G4VisAttributes(false));

		new G4PVPlacement(0,						//no rotation
						  G4ThreeVector(),			//at position
						  logicCloud,				//its logical volume
						  "cloudPV",				//its name
						  logicWorld,				//its mothervolume
						  false,					//no boolean operation
						  0,						//copy number
						  checkOverlaps);			//overlaps checking

		G4LogicalVolume* logicDroplet = 
			new G4LogicalVolume(fDropletShape,		//its solid
								droplet_mat,		//its material
								"dropletLV");		//its name

		logicDroplet->SetUserLimits(fStepLimits);

		/*G4PVParameterised* paramDroplet =*/
			new G4PVParameterised("droplets",		//its name
								  logicDroplet,		//droplet logical volume
								  logicCloud,		//mother logical volume
								  kUndefined,		//droplets placed along this axis
								  positions.size(),	//number of droplets
								  cloudParam);		//the parametrisation 
	}
	// **********************************************************
	// 
	// (For comparision/benchmarking) Simulate the cloud by smearing droplets out into a single solid (does not use FastAerosol)
	// 
	// ***********************************************************
	else if (fSmoothCloud)
	{
		G4cout << "\nSmooth geometry based on a cloud of n=" << fDropletNumDens*mm3 << "/mm3 and r=" << fDropletR/mm << "mm spheres.\n" << G4endl;
		// build cloud by smearing the droplets uniformly across the cloud volume, for comparison/benchmarking purposes (does not use FastAerosol)
		G4Material* cloud_mat = new G4Material("Cloud", (droplet_total_mass+cloud_air_mass)/cloud_volume, 2);
		cloud_mat->AddMaterial(droplet_mat, droplet_total_mass/(cloud_air_mass+droplet_total_mass));
		cloud_mat->AddMaterial(air_mat, cloud_air_mass/(cloud_air_mass+droplet_total_mass));
							
		logicCloud =						 
			new G4LogicalVolume(fCloudShape,	//its solid
								cloud_mat,	//its material
								"cloudLV");	//its name
		logicCloud->SetUserLimits(fStepLimits);
		logicCloud->SetVisAttributes(G4VisAttributes(G4Colour(0.0,0.0,1.0,0.4)));
					 
		new G4PVPlacement(0,				//no rotation
						G4ThreeVector(),	//at position
						logicCloud,			//its logical volume
						"cloudPV",			//its name
						logicWorld,			//its mothervolume
						false,				//no boolean operation
						0,					//copy number
						checkOverlaps);		//overlaps checking
	}
	else
	{
		G4cout << "\nNo cloud.\n" << G4endl;
	}	

	//
	// Build detector
	//
	G4double detector_sizeXY = cloud_sizeXY;
	G4double detector_sizeZ = 0.05*m;
	G4Material* detector_mat = nist->FindOrBuildMaterial("G4_Al");
	G4ThreeVector detector_pos = G4ThreeVector(0, 0, 0.5*1.05*cloud_sizeZ);

	G4Box* soldDetector =	
		new G4Box("detectorSV",				//its name
				0.5*detector_sizeXY,		//half x-span
				0.5*detector_sizeXY,		//half y-span
				0.5*detector_sizeZ);		//half z-span
						
	G4LogicalVolume* logicDetector =						 
		new G4LogicalVolume(soldDetector,	//its solid
							detector_mat,	//its material
							"detectorLV");	//its name

	logicDetector->SetUserLimits(fStepLimits);
				 
	new G4PVPlacement(0,					//no rotation
						detector_pos,		//at position
						logicDetector,		//its logical volume
						"detectorPV",		//its name
						logicWorld,			//its mothervolume
						false,				//no boolean operation
						0,					//copy number
						checkOverlaps);		//overlaps checking


	//
	// Scoring Volume
	//
	fScoringVolume = logicDetector;

	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
