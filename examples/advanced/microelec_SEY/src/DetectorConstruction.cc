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
/// \file electromagnetic/TestEm14/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 92996 2015-09-28 08:03:50Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4NistManager.hh"
#include <vector>
#include "G4ProductionCuts.hh"
#include "G4Region.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//vector <G4VPhysicalVolume*> All_Detectors_info;

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),World_Material(NULL), fSolidWorld(NULL), fLogicWorld(NULL), fPhysiWorld(NULL)
, fLBox(0),fPBox(0), fMaterial(0),detectorMat(NULL), fDetectorMessenger(0)
{
  fBoxSize = 1*mm;
  fBoxWidth = 1*mm;
  fBoxSizeSurface = 0.1*fBoxSize;
  fBoxSizeLayer1 = 0.1*fBoxSize;
  fBoxSizeLayer2 = 0.1*fBoxSize;
  fBoxSizeLayer3 = 0.1*fBoxSize;
  fBoxSizeLayer4 = 0.1*fBoxSize;
  fWorldSizeX = fBoxSize *5;
  fWorldSizeY = fBoxSize *5;
  fWorldSizeZ = (fBoxWidth + fBoxSizeLayer4 + fBoxSizeLayer3 + fBoxSizeLayer2 + fBoxSizeLayer1+ fBoxSizeSurface) * 5;
  WorldDim = std::max(fWorldSizeX, fWorldSizeZ);
  WorldRay = WorldDim / 2;
  DetectorRay = WorldRay*0.8;

  DefineMaterials();
  SetMaterial("G4_Si");
  SetMaterialSurface("G4_Si");
  SetMaterialLayer1("G4_Si");
  SetMaterialLayer2("G4_Si");
  SetMaterialLayer3("G4_Si");
  SetMaterialLayer4("G4_Si");

  fDetectorMessenger = new DetectorMessenger(this);

  //pSDTrajectoires = NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
	delete fDetectorMessenger;
	}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	DefineMaterials();
	return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4double a,z,density,fractionmass;
  G4double pressure, temperature;
  G4String name,symbol;
  G4int nel;

  a=1.01*g/mole;
  G4Element *elH=new G4Element(name="Hydrogen",symbol="H",z=1.,a);
  a=12.01*g/mole;
  G4Element* elB = new G4Element(name = "Boron", symbol = "B", z = 5., a);
  a = 10.811 * g / mole;
  G4Element *elC=new G4Element(name="Carbon",symbol="C",z=6.,a);
  a=14.01*g/mole;
  G4Element *elN=new G4Element(name="Nitrogen",symbol="N",z=7.,a);
  a=16.*g/mole;
  G4Element *elO=new G4Element(name="Oxygen",symbol="O",z=8.,a);
  a=28.0855*g/mole;
  G4Element* elSi = new G4Element("Silicon", "Si", z=14., a);
  a = 47.9 * g / mole;
  G4Element* elTi = new G4Element(name = "Titanium", symbol = "Ti", z = 22., a);

  
  // G4_SILICON_DIOXIDE and G4_KAPTON forced declaration and present in NIST database
  // G4TITANIUM_NITRIDE and G4_BORON_NITRIDE do not exist in the G4 NIST database
  My_SiO2 = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  My_Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
  My_BN = G4NistManager::Instance()->FindOrBuildMaterial("G4_BORON_NITRIDE");
  My_TiN = G4NistManager::Instance()->FindOrBuildMaterial("G4_TITANIUM_NITRIDE");
  
  if (!My_SiO2) {
        density = 2.32*g/cm3;
        My_SiO2 = new G4Material("G4_SILICON_DIOXIDE", density, 2);
        My_SiO2->AddElement(elSi, 1);
        My_SiO2->AddElement(elO , 2);
    }
  if (!My_Kapton) {
        density = 1.42*g/cm3;
		My_Kapton = new G4Material(name="G4_KAPTON",density, nel=4);
		My_Kapton->AddElement(elH, fractionmass = 0.0273);
		My_Kapton->AddElement(elC, fractionmass = 0.7213);
		My_Kapton->AddElement(elN, fractionmass = 0.0765);
		My_Kapton->AddElement(elO, fractionmass = 0.1749);
    }
  if (!My_BN) {
		density = 2.1 * g / cm3;
		My_BN = new G4Material(name = "G4_BORON_NITRIDE", density, nel = 2);
		My_BN->AddElement(elB, fractionmass = 0.436);
		My_BN->AddElement(elN, fractionmass = 0.564);
     } 
  if (!My_TiN) {
		density = 5.4*g / cm3;
		My_TiN = new G4Material(name = "G4_TITANIUM_NITRIDE", density, nel = 2);
		My_TiN->AddElement(elTi, fractionmass = 0.7737);
		My_TiN->AddElement(elN, fractionmass = 0.2263);
     } 

  // Vaccum must be created for MicroElec, should use G4_GALACTIC and changing the name to Vacuum
  a = 1.01*g / mole;
  pressure = 1.e-19*pascal;
  temperature = 0.1*kelvin;
  density = universe_mean_density;

  G4Material* Vide = new G4Material(name = "Vacuum", z = 1., a, density, kStateGas, temperature, pressure);

  // LayerX materials, Surface material and Substrate materials are initialized in the detector construction
  World_Material = Vide;
  detectorMat = Vide;
  DefaultMaterial = My_TiN;
  fMaterial = My_BN;
  fMaterialSurface = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
  fMaterialLayer1 = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
  fMaterialLayer2 = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
  fMaterialLayer3 = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
  fMaterialLayer4 = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry and spectrum variable
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  G4cout << "MicroELec - SEY modellling - Creating geometry" << G4endl;
  // world volume                         
  fWorldSizeX = fBoxSize * 5;
  fWorldSizeY = fBoxSize * 5;
  fWorldSizeZ = (fBoxWidth + fBoxSizeLayer4 + fBoxSizeLayer3 + fBoxSizeLayer2 + fBoxSizeLayer1+ fBoxSizeSurface) * 5;
  WorldDim = std::max(fWorldSizeX, fWorldSizeZ);
  WorldRay = WorldDim / 2;
  DetectorRay = WorldRay * 0.8;

  fSolidWorld = new G4Box("World", WorldDim/2.0, WorldDim / 2.0, WorldDim / 2.0);
  //fLogicWorld = new G4LogicalVolume(fSolidWorld, World_Material, "World", 0, 0, 0);
  fLogicWorld = new G4LogicalVolume(fSolidWorld, World_Material, "World");
  fPhysiWorld = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), "World", fLogicWorld, 0, false, 0);
  

  // construct geometry substrate volume 
  
  fSBox = new G4Box("Substrate",                        //its name
	  fBoxSize , fBoxSize ,fBoxWidth/2);  				//its dimensions

  fLBox = new G4LogicalVolume(fSBox,                    //its shape
								fMaterial,             //its material
								"Substrate");          //fMaterial->GetName());    //its name

  fPBox = new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0.0, (0*fBoxSize), fBoxSizeSurface + fBoxSizeLayer1 + fBoxSizeLayer2 + fBoxSizeLayer3 + fBoxSizeLayer4 + fBoxWidth / 2 ),//(fBoxSize / 2 + fBoxSizeSurface)),            //at (0,0,fBoxSize/2+fBoxSizeSurface)
							"Substrate", 				//fMaterial->GetName() its name
							fLBox,                      //its logical volume
							fPhysiWorld,                //its mother  volume
							false,                      //no boolean operation
							0);                         //copy number//*/
  
// construct geometry surface volume 
  
  fSBoxSurface = new G4Box("Surface",                   //its name
		  fBoxSize , fBoxSize , fBoxSizeSurface / 2);  //its dimensions X & Y leteral dimensions = sBox

  fLBoxSurface = new G4LogicalVolume(fSBoxSurface,                  //its shape
										fMaterialSurface,           //its material
										"Surface");					// fMaterialSurface->GetName());    //its name

  fPBoxSurface = new G4PVPlacement(0,                         //no rotation
									G4ThreeVector(0.0, 0.0, fBoxSizeSurface/2),            //at (0,0,0) on the top of the sbox surface at Z= 0.
									"Surface", 				    //	fMaterialSurface->GetName(), its name
									fLBoxSurface,               //its logical volume
									fPhysiWorld,                //its mother  volume
									false,                      //no boolean operation
									0);                         //copy number

  
	  // construct geometry Layer1 volume 
  
  fSBoxLayer1 = new G4Box("Layer1",                 //its name
		  fBoxSize, fBoxSize, fBoxSizeLayer1 / 2);  //its dimensions X & Y leteral dimensions = sBox

  fLBoxLayer1 = new G4LogicalVolume(fSBoxLayer1,               //its shape
									fMaterialLayer1,           //its material
									"Layer1");				   // fMaterialSurface->GetName());    //its name

  fPBoxLayer1 = new G4PVPlacement(0,                         //no rotation
									G4ThreeVector(0.0, 0.0, fBoxSizeSurface + fBoxSizeLayer1 / 2),            //at (0,0,0) on the top of the sbox surface at Z= 0.
									"Layer1",                   //	fMaterialSurface->GetName(),       //its name
									fLBoxLayer1,                //its logical volume
									fPhysiWorld,                //its mother  volume
									false,                      //no boolean operation
									0);                         //copy number
	// construct geometry Layer2 volume 

  fSBoxLayer2 = new G4Box("Layer2",                 //its name
		  fBoxSize, fBoxSize, fBoxSizeLayer2 / 2);  //its dimensions X & Y leteral dimensions = sBox

  fLBoxLayer2 = new G4LogicalVolume(fSBoxLayer2,               //its shape
									fMaterialLayer2,           //its material
									"Layer2");				// fMaterialSurface->GetName());    //its name

  fPBoxLayer2 = new G4PVPlacement(0,                         //no rotation
									G4ThreeVector(0.0, 0.0, fBoxSizeSurface + fBoxSizeLayer1 +fBoxSizeLayer2 / 2),            //at (0,0,0) on the top of the sbox surface at Z= 0.
									"Layer2",                   //fMaterialSurface->GetName(),       //its name
									fLBoxLayer2,                //its logical volume
									fPhysiWorld,                //its mother  volume
									false,                      //no boolean operation
									0);                         //copy number

	// construct geometry Layer3 volume 
  
  fSBoxLayer3 = new G4Box("Layer3",                        //its name
		  fBoxSize, fBoxSize, fBoxSizeLayer3 / 2);         //its dimensions X & Y leteral dimensions = sBox

  fLBoxLayer3 = new G4LogicalVolume(fSBoxLayer3,               //its shape
									fMaterialLayer3,           //its material
									"Layer3");				// fMaterialSurface->GetName());    //its name

  fPBoxLayer3 = new G4PVPlacement(0,                         //no rotation
									G4ThreeVector(0.0, 0.0, fBoxSizeSurface + fBoxSizeLayer1 + fBoxSizeLayer2 + fBoxSizeLayer3 / 2),            //at (0,0,0) on the top of the sbox surface at Z= 0.
									"Layer3",                   //fMaterialSurface->GetName(),       //its name
									fLBoxLayer3,                //its logical volume
									fPhysiWorld,                //its mother  volume
									false,                      //no boolean operation
									0);                         //copy number

	// construct geometry Layer4 volume 
  
  fSBoxLayer4 = new G4Box("Layer4",                        //its name
		  fBoxSize, fBoxSize, fBoxSizeLayer4 / 2);         //its dimensions X & Y leteral dimensions = sBox

  fLBoxLayer4 = new G4LogicalVolume(fSBoxLayer4,               //its shape
									fMaterialLayer4,           //its material
									"Layer4");				   // fMaterialSurface->GetName());    //its name

  fPBoxLayer4 = new G4PVPlacement(0,                         //no rotation
									G4ThreeVector(0.0, 0.0, fBoxSizeSurface + fBoxSizeLayer1 + fBoxSizeLayer2 + fBoxSizeLayer3 + fBoxSizeLayer4 / 2),            //at (0,0,0) on the top of the sbox surface at Z= 0.
									"Layer4",                //fMaterialSurface->GetName(),       //its name
									fLBoxLayer4,             //its logical volume
									fPhysiWorld,             //its mother  volume
									false,                   //no boolean operation
									0);                      //copy number



	 // ---------------------------- - detector 1
	 // ----------------------------------------

  Detector_1_box = new G4Sphere("Detector_1", DetectorRay, DetectorRay + 1*angstrom, 0.0, 2.0*3.1415, 0.0, 3.1415);
  Detector_1_log = new G4LogicalVolume(Detector_1_box, detectorMat, "Detector_1");
  Detector_1_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), "Detector_1", Detector_1_log, fPhysiWorld, false, 0);//*/


  fRegion = new G4Region("Target");
  G4ProductionCuts* cuts = new G4ProductionCuts();


  // for detailed transport 
  G4double defCut = 1 * nanometer;
  cuts->SetProductionCut(defCut, "gamma");
  cuts->SetProductionCut(defCut, "e-");
  cuts->SetProductionCut(defCut, "e+");
  cuts->SetProductionCut(defCut, "proton");

  fRegion->SetProductionCuts(cuts);
  fRegion->AddRootLogicalVolume(fLBoxSurface);
  fRegion->AddRootLogicalVolume(fLBox);
  fRegion->AddRootLogicalVolume(fLBoxLayer1);
  fRegion->AddRootLogicalVolume(fLBoxLayer2);
  fRegion->AddRootLogicalVolume(fLBoxLayer3);
  fRegion->AddRootLogicalVolume(fLBoxLayer4);


  G4VisAttributes* subVisAtt = new G4VisAttributes(G4Colour(1, 0, 0)); //Red, //subVisAtt->SetForceSolid(true);											  
  fLBoxSurface->SetVisAttributes(subVisAtt);
  fLBoxLayer1->SetVisAttributes(subVisAtt);
  fLBoxLayer2->SetVisAttributes(subVisAtt);
  fLBoxLayer3->SetVisAttributes(subVisAtt);
  fLBoxLayer4->SetVisAttributes(subVisAtt);
  

  G4VisAttributes* subVisAtt2 = new G4VisAttributes(G4Colour(0, 1, 0)); //green		//subVisAtt->SetForceSolid(true);								    
  fLBox->SetVisAttributes(subVisAtt2); 
  
  G4VisAttributes* subVisAtt3 = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0)); //White
  subVisAtt3->SetVisibility(true);
  Detector_1_log->SetVisAttributes(subVisAtt3);


  G4cout << "MicroELec - SEY modellling - End geometry creation " << G4endl;
  PrintParameters();
  G4cout << fPhysiWorld->GetLogicalVolume()->GetMaterial()->GetName() << G4endl;
  G4cout << fPhysiWorld->GetLogicalVolume()->GetDaughter(0)->GetLogicalVolume()->GetMaterial()->GetName() << G4endl;
  return fPhysiWorld;

  //always return the root volume
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
	// sensitive detectors -----------------------------------------------------
	auto sdManager = G4SDManager::GetSDMpointer();
	G4String SDname;
	
	auto MypMicroElecSdSey = new MicroElecSdSey(SDname = "/SdSey", "DetecteurMicroElecSdSey");
	sdManager->AddNewDetector(MypMicroElecSdSey);
	Detector_1_log->SetSensitiveDetector(MypMicroElecSdSey);
	G4cout<<"Sensitive Detector Name : "<<MypMicroElecSdSey->GetName()<<G4endl;
}







//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(fBoxSize,"Length")
         << " of " << fMaterial->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
	G4Material* pttoMaterial = nullptr;
  // search the material by its name, or build it from nist data base
	if (materialChoice == "G4_SiO2"|| materialChoice == "G4_SILICON_DIOXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	}
	else if (materialChoice == "G4_KAPTON"|| materialChoice == "G4_KAPTON_FILM"|| materialChoice == "Kapton") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
	}
	else if (materialChoice == "G4_Al2O3"|| materialChoice == "G4_ALUMINUM_OXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	}
	else if (materialChoice == "G4_TiN" || materialChoice == "G4_TITANIUM_NITRIDE") {
		pttoMaterial = My_TiN;
		pttoMaterial->SetName("G4_TITANIUM_NITRIDE");
	}
	else if (materialChoice == "G4_BN" || materialChoice == "G4_BORON_NITRIDE") {
		pttoMaterial = My_BN;
		pttoMaterial->SetName("G4_BORON_NITRIDE");
	}
	else { pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice); }

  if (pttoMaterial) {
   fMaterial = pttoMaterial;
   G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
	  G4String str = "Material ";
	  str += materialChoice + " does not exist";
	  G4Exception("DetectorConstruction::SetMaterial", "em0002", FatalException, str);
  } 
}

void DetectorConstruction::SetMaterialSurface(G4String materialChoice)
{
	G4Material* pttoMaterial = nullptr;
	// search the material by its name, or build it from nist data base
	if (materialChoice == "G4_SiO2"|| materialChoice == "G4_SILICON_DIOXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	}
	else if (materialChoice == "G4_KAPTON" || materialChoice == "G4_KAPTON_FILM" || materialChoice == "Kapton") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
	}
	else if (materialChoice == "G4_Al2O3"|| materialChoice == "G4_ALUMINUM_OXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	}
	else if (materialChoice == "G4_TiN" || materialChoice == "G4_TITANIUM_NITRIDE") {
		pttoMaterial = My_TiN;
		pttoMaterial->SetName("G4_TITANIUM_NITRIDE");
	}
	else if (materialChoice == "G4_BN" || materialChoice == "G4_BORON_NITRIDE") {
		pttoMaterial = My_BN;
		pttoMaterial->SetName("G4_BORON_NITRIDE");
	}
	else { pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice); }

	if (pttoMaterial) {
		fMaterialSurface = pttoMaterial;
		G4RunManager::GetRunManager()->PhysicsHasBeenModified();
	}
	else {
		G4String str = "Material ";
		str += materialChoice + " does not exist";
		G4Exception("DetectorConstruction::SetMaterialSurface", "em0002", FatalException, str);
	}
}

void DetectorConstruction::SetMaterialLayer1(G4String materialChoice)
{
	G4Material* pttoMaterial = nullptr;
	// search the material by its name, or build it from nist data base
	if (materialChoice == "G4_SiO2"|| materialChoice == "G4_SILICON_DIOXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	}
	else if (materialChoice == "G4_KAPTON" || materialChoice == "G4_KAPTON_FILM" || materialChoice == "Kapton") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
	}
	else if (materialChoice == "G4_Al2O3"|| materialChoice == "G4_ALUMINUM_OXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	}
	else if (materialChoice == "G4_TiN" || materialChoice == "G4_TITANIUM_NITRIDE") {
		pttoMaterial = My_TiN;
		pttoMaterial->SetName("G4_TITANIUM_NITRIDE");
	}
	else if (materialChoice == "G4_BN" || materialChoice == "G4_BORON_NITRIDE") {
		pttoMaterial = My_BN;
		pttoMaterial->SetName("G4_BORON_NITRIDE");
	}
	else { pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice); }

	if (pttoMaterial) {
		fMaterialLayer1 = pttoMaterial;
		G4RunManager::GetRunManager()->PhysicsHasBeenModified();
	}
	else {
		G4String str = "Material ";
		str += materialChoice + " does not exist";
		G4Exception("DetectorConstruction::SetMaterialSurface", "em0002", FatalException, str);
	}
}

void DetectorConstruction::SetMaterialLayer2(G4String materialChoice)
{
	G4Material* pttoMaterial = nullptr;
	// search the material by its name, or build it from nist data base
	if (materialChoice == "G4_SiO2"|| materialChoice == "G4_SILICON_DIOXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	}
	else if (materialChoice == "G4_KAPTON" || materialChoice == "G4_KAPTON_FILM" || materialChoice == "Kapton") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
	}
	else if (materialChoice == "G4_Al2O3"|| materialChoice == "G4_ALUMINUM_OXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	}
	else if (materialChoice == "G4_TiN" || materialChoice == "G4_TITANIUM_NITRIDE") {
		pttoMaterial = My_TiN;
		pttoMaterial->SetName("G4_TITANIUM_NITRIDE");
	}
	else if (materialChoice == "G4_BN" || materialChoice == "G4_BORON_NITRIDE") {
		pttoMaterial = My_BN;
		pttoMaterial->SetName("G4_BORON_NITRIDE");
	}
	else { pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice); }

	if (pttoMaterial) {
		fMaterialLayer2 = pttoMaterial;
		G4RunManager::GetRunManager()->PhysicsHasBeenModified();
	}
	else {
		G4String str = "Material ";
		str += materialChoice + " does not exist";
		G4Exception("DetectorConstruction::SetMaterialSurface", "em0002", FatalException, str);
	}
}

void DetectorConstruction::SetMaterialLayer3(G4String materialChoice)
{
	G4Material* pttoMaterial = nullptr;
	// search the material by its name, or build it from nist data base
	if (materialChoice == "G4_SiO2"|| materialChoice == "G4_SILICON_DIOXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	}
	else if (materialChoice == "G4_KAPTON" || materialChoice == "G4_KAPTON_FILM" || materialChoice == "Kapton") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
	}
	else if (materialChoice == "G4_Al2O3"|| materialChoice == "G4_ALUMINUM_OXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	}
	else if (materialChoice == "G4_TiN" || materialChoice == "G4_TITANIUM_NITRIDE") {
		pttoMaterial = My_TiN;
		pttoMaterial->SetName("G4_TITANIUM_NITRIDE");
	}
	else if (materialChoice == "G4_BN" || materialChoice == "G4_BORON_NITRIDE") {
		pttoMaterial = My_BN;
		pttoMaterial->SetName("G4_BORON_NITRIDE");
	}
	else { pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice); }

	if (pttoMaterial) {
		fMaterialLayer3 = pttoMaterial;
		G4RunManager::GetRunManager()->PhysicsHasBeenModified();
	}
	else {
		G4String str = "Material ";
		str += materialChoice + " does not exist";
		G4Exception("DetectorConstruction::SetMaterialSurface", "em0002", FatalException, str);
	}
}

void DetectorConstruction::SetMaterialLayer4(G4String materialChoice)
{
	G4Material* pttoMaterial = nullptr;
	// search the material by its name, or build it from nist data base
	if (materialChoice == "G4_SiO2"|| materialChoice == "G4_SILICON_DIOXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	}
	else if (materialChoice == "G4_KAPTON" || materialChoice == "G4_KAPTON_FILM" || materialChoice == "Kapton") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
	}
	else if (materialChoice == "G4_Al2O3"|| materialChoice == "G4_ALUMINUM_OXIDE") {
		pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	}
	else if (materialChoice == "G4_TiN" || materialChoice == "G4_TITANIUM_NITRIDE") {
		pttoMaterial = My_TiN;
		pttoMaterial->SetName("G4_TITANIUM_NITRIDE");
	}
	else if (materialChoice == "G4_BN" || materialChoice == "G4_BORON_NITRIDE") {
		pttoMaterial = My_BN;
		pttoMaterial->SetName("G4_BORON_NITRIDE");
	}
	else { pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice); }

	if (pttoMaterial) {
		fMaterialLayer4 = pttoMaterial;
		G4RunManager::GetRunManager()->PhysicsHasBeenModified();
	}
	else {
		G4String str = "Material ";
		str += materialChoice + " does not exist";
		G4Exception("DetectorConstruction::SetMaterialSurface", "em0002", FatalException, str);
	}
}




void DetectorConstruction::UpdateGeometry()
{
	G4cout << "UpdateGeometry Begin" << G4endl;
	G4GeometryManager::GetInstance()->OpenGeometry();

	fWorldSizeX = fBoxSize * 5;
	fWorldSizeY = fBoxSize * 5;
	fWorldSizeZ = (fBoxWidth + fBoxSizeLayer4 + fBoxSizeLayer3 + fBoxSizeLayer2 + fBoxSizeLayer1+ fBoxSizeSurface) * 5;
	WorldDim = std::max(fWorldSizeX, fWorldSizeZ);
	WorldRay = WorldDim / 2;
	DetectorRay = WorldRay * 0.8;
	
	if (fSolidWorld) {
		fSolidWorld->SetXHalfLength(WorldDim / 2);
		fSolidWorld->SetYHalfLength(WorldDim / 2);
		fSolidWorld->SetZHalfLength(WorldDim / 2);
	}

	if (fLBox) { fLBox->SetMaterial(fMaterial);}
	if (fLBoxSurface) { fLBoxSurface->SetMaterial(fMaterialSurface); }
	if (fLBoxLayer1) { fLBoxLayer1->SetMaterial(fMaterialLayer1); }
	if (fLBoxLayer2) { fLBoxLayer2->SetMaterial(fMaterialLayer2); }
	if (fLBoxLayer3) { fLBoxLayer3->SetMaterial(fMaterialLayer3); }
	if (fLBoxLayer4) { fLBoxLayer4->SetMaterial(fMaterialLayer4); }

	if (fSBox) { 
		fSBox->SetZHalfLength(fBoxWidth / 2);
		fSBox->SetXHalfLength(fBoxSize / 2);
		fSBox->SetYHalfLength(fBoxSize / 2);
		fPBox->SetTranslation(G4ThreeVector(0.0, 0.0, fBoxSizeSurface + fBoxSizeLayer1 + fBoxSizeLayer2 + fBoxSizeLayer3 + fBoxSizeLayer4 + fBoxSize / 2));
	}
	if (fSBoxSurface) { 
		fSBoxSurface->SetZHalfLength(fBoxSizeSurface / 2);
		fSBoxSurface->SetXHalfLength(fBoxSize / 2);
		fSBoxSurface->SetYHalfLength(fBoxSize / 2);
		fPBoxSurface->SetTranslation(G4ThreeVector(0.0, 0.0, fBoxSizeSurface / 2));
	}
	if (fSBoxLayer1) { 
		fSBoxLayer1->SetZHalfLength(fBoxSizeLayer1 / 2); 
		fSBoxLayer1->SetXHalfLength(fBoxSize / 2);
		fSBoxLayer1->SetYHalfLength(fBoxSize / 2);
		fPBoxLayer1->SetTranslation(G4ThreeVector(0.0, 0.0, fBoxSizeSurface + fBoxSizeLayer1 / 2));
	}
	if (fSBoxLayer2) { 
		fSBoxLayer2->SetZHalfLength(fBoxSizeLayer2 / 2);
		fSBoxLayer2->SetXHalfLength(fBoxSize / 2);
		fSBoxLayer2->SetYHalfLength(fBoxSize / 2);
		fPBoxLayer2->SetTranslation(G4ThreeVector(0.0, 0.0, fBoxSizeSurface + fBoxSizeLayer1 + fBoxSizeLayer2 / 2));
	}
	if (fSBoxLayer3) { 
		fSBoxLayer3->SetZHalfLength(fBoxSizeLayer3 / 2); 
		fSBoxLayer3->SetXHalfLength(fBoxSize / 2);
		fSBoxLayer3->SetYHalfLength(fBoxSize / 2);
		fPBoxLayer3->SetTranslation(G4ThreeVector(0.0, 0.0, fBoxSizeSurface + fBoxSizeLayer1 + fBoxSizeLayer2 + fBoxSizeLayer3 / 2));
	}
	if (fSBoxLayer4) {
		fSBoxLayer4->SetZHalfLength(fBoxSizeLayer4 / 2);
		fSBoxLayer4->SetXHalfLength(fBoxSize / 2);
		fSBoxLayer4->SetYHalfLength(fBoxSize / 2);
		fPBoxLayer4->SetTranslation(G4ThreeVector(0.0, 0.0, fBoxSizeSurface + fBoxSizeLayer1 + fBoxSizeLayer2 + fBoxSizeLayer3 + fBoxSizeLayer4 / 2));
	}
	
	if (Detector_1_box) {
		Detector_1_box->SetInnerRadius(DetectorRay);
		Detector_1_box->SetOuterRadius(DetectorRay + 1 * angstrom);
	}
	G4RunManager::GetRunManager()->PhysicsHasBeenModified();
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
	

	G4cout << "Substrate:" << this->GetMaterial()->GetName()<< ", 4:" << this->GetMaterialLayer4()->GetName() << ", 3:" << this->GetMaterialLayer3()->GetName() << ", 2:" << this->GetMaterialLayer2()->GetName() << ", 1:" << this->GetMaterialLayer1()->GetName() << ", surf:" << this->GetMaterialSurface()->GetName() << G4endl;
	G4cout << "Substrate:" << this->GetSize() << ", 4:" << this->GetSizeLayer4()<< ", 3:" << this->GetSizeLayer3() << ", 2:" << this->GetSizeLayer2() << ", 1:" << this->GetSizeLayer1() << ", surf:" << this->GetSizeSurface()<< G4endl;
	G4cout << "UpdateGeometry End" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fBoxSize = value;
}

void DetectorConstruction::SetWidth(G4double value)
{
	fBoxWidth = value;
	//G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetSizeSurface(G4double value)
{
	fBoxSizeSurface = value;
	//G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetSizeLayer1(G4double value)
{
	fBoxSizeLayer1 = value;
	//G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetSizeLayer2(G4double value)
{
	fBoxSizeLayer2 = value;
	//G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetSizeLayer3(G4double value)
{
	fBoxSizeLayer3 = value;
	//G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetSizeLayer4(G4double value)
{
	fBoxSizeLayer4 = value;
	//G4RunManager::GetRunManager()->ReinitializeGeometry();
}

