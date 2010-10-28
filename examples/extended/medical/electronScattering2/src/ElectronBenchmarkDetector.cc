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

#include "ElectronBenchmarkDetector.hh"

#include "ElectronBenchmarkDetectorMessenger.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4NistManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4SDParticleFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSCellFlux.hh"
#include "G4PSPopulation.hh"

ElectronBenchmarkDetector::ElectronBenchmarkDetector()
:fLogWorld(NULL),
worldVisAtt(0),
windowVisAtt(0),
primFoilVisAtt(0), 
monVisAtt(0),
bagVisAtt(0),
heliumVisAtt(0),
ringVisAtt(0),
scorerVisAtt(0)
{
	// Exit Window
	fPosWindow0		=   0.000000*cm;
	fPosWindow1		=   0.004120*cm;
	
	// Scattering Foil
	fPosPrimFoil	=   2.650000*cm;
	
	// Monitor Chamber
	fPosMon0		=   5.000000*cm;
	fPosMon1		=   5.011270*cm;
	
	// Helium Bag
	fPosBag0		=   6.497500*cm;
	fPosHelium0		=   6.500000*cm;
	fPosHelium1		= 116.500000*cm;
	fPosBag1		= 116.502500*cm;
	fThicknessRing	=   1.4*cm;
	
	// Scoring Plane
	fPosScorer		= 118.200000*cm;	
	fThicknessScorer=   0.001*cm;
	fWidthScorerRing=   0.1*cm;
	
	// Radii
	fRadOverall		=  23.3*cm;
	fRadRingInner	=  20.0*cm;
	
	// Extra space remaining in world volume around apparatus
	fPosDelta		=   1.*cm;
	fRadDelta		=   0.1*cm;

	fMessenger = new ElectronBenchmarkDetectorMessenger(this);		
}


ElectronBenchmarkDetector::~ElectronBenchmarkDetector()
{
	delete fMessenger;
	
	delete worldVisAtt;
	delete windowVisAtt;
	delete primFoilVisAtt;
	delete monVisAtt;
	delete bagVisAtt;
	delete heliumVisAtt;
	delete ringVisAtt;
	delete scorerVisAtt;
}


G4VPhysicalVolume* ElectronBenchmarkDetector::Construct()
{ 
	DefineMaterials();
	
	G4VPhysicalVolume* physiworld = CreateGeometry();
	return physiworld;
}


void ElectronBenchmarkDetector::DefineMaterials(){ 
	// Use NIST database for elements and materials whereever possible.
	G4NistManager* man = G4NistManager::Instance();
	man->SetVerbose(1);
	
	// Take all elements and materials from NIST 
	man->FindOrBuildMaterial("G4_He");
	man->FindOrBuildMaterial("G4_Be");
	man->FindOrBuildMaterial("G4_Al");
	man->FindOrBuildMaterial("G4_Ti");
	man->FindOrBuildMaterial("G4_Ta");
	man->FindOrBuildMaterial("G4_AIR");
	man->FindOrBuildMaterial("G4_MYLAR");
	
	G4Element* C  = man->FindOrBuildElement("C");
	G4Element* Cu = man->FindOrBuildElement("Cu");
	G4Element* Au = man->FindOrBuildElement("Au");
	G4Element* Ti = man->FindOrBuildElement("Ti");
	G4Element* Al = man->FindOrBuildElement("Al");
	G4Element* V  = man->FindOrBuildElement("V");
	
	// Define materials not in NIST
	// While the NIST database does contain default materials for C, Cu and Au, those defaults have different
	// densities than the ones used in the benchmark specification.
	G4double density;
	G4int ncomponents;
	G4double fractionmass;
	
	G4Material* G4_C = new G4Material("G4_C", density= 2.18*g/cm3, ncomponents=1);
	G4_C->AddElement(C, fractionmass=1.00);

	G4Material* G4_Cu = new G4Material("G4_Cu", density= 8.92*g/cm3, ncomponents=1);
	G4_Cu->AddElement(Cu, fractionmass=1.00);

	G4Material* G4_Au = new G4Material("G4_Au", density= 19.30*g/cm3, ncomponents=1);
	G4_Au->AddElement(Au, fractionmass=1.00);

	G4Material* TiAlloy = new G4Material("TiAlloy", density= 4.42*g/cm3, ncomponents=3);
	TiAlloy->AddElement(Ti, fractionmass=0.90);
	TiAlloy->AddElement(Al, fractionmass=0.06);
	TiAlloy->AddElement(V,  fractionmass=0.04);

	// Print materials table	
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


G4VPhysicalVolume* ElectronBenchmarkDetector::CreateGeometry(){	
	// Clean old geometry, if any
	fLogWorld=NULL;
	fScorerRingLog=NULL;
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();
	
	// Instantiate the world
	G4VPhysicalVolume* physiworld = CreateWorld();
	fLogWorld = physiworld->GetLogicalVolume();
	
	// Instantiate the geometry
	CreateExitWindow(fLogWorld);
	CreatePrimaryFoil(fLogWorld);
	CreateMonitor(fLogWorld);
	CreateHeliumBag(fLogWorld);
	
	// Create and activate the scorers
	CreateScorer(fLogWorld);
	ActivateScorer();
	
	return physiworld;	
}


void ElectronBenchmarkDetector::UpdateGeometry(){
	G4RunManager::GetRunManager()->DefineWorldVolume(CreateGeometry());
}


G4VPhysicalVolume* ElectronBenchmarkDetector::CreateWorld(){
	G4double halfLengthWorld = fPosScorer/2. + fPosDelta;
	G4double radWorld = fRadOverall + fRadDelta;
	G4VSolid* worldSolid = new G4Tubs("WorldSolid", 0.*cm, radWorld, halfLengthWorld, 0.*deg, 360.*deg);
	G4LogicalVolume* worldLog = new G4LogicalVolume(worldSolid, G4Material::GetMaterial("G4_AIR"), "WorldLog");
	
	worldVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	worldLog->SetVisAttributes(worldVisAtt);
	
	G4VPhysicalVolume* worldPhys = 
		new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), worldLog,"World", 0, false, 0);
	
    return worldPhys;
}


void ElectronBenchmarkDetector::CreateExitWindow(G4LogicalVolume* worldLog){
	G4double halfLengthWorld = fPosScorer/2.;
	G4double halfThicknessWindow = fPosWindow1/2.;
	G4VSolid* windowSolid = new G4Tubs("windowSolid", 0.*cm, fRadOverall, halfThicknessWindow, 0.*deg, 360.*deg);
	G4LogicalVolume* windowLog = new G4LogicalVolume(windowSolid, G4Material::GetMaterial("TiAlloy"), "windowLog");
	
	windowVisAtt = new G4VisAttributes(G4Colour(0.5,1.0,0.5));
	windowLog->SetVisAttributes(windowVisAtt);
	
	new G4PVPlacement(0, 
		G4ThreeVector(0.,0., halfThicknessWindow - halfLengthWorld),
		windowLog,"ExitWindow",worldLog,false,0);	
}


void ElectronBenchmarkDetector::CreatePrimaryFoil(G4LogicalVolume* worldLog){
	G4double halfLengthWorld = fPosScorer/2.;

	// For some energies, we have no Primary Foil.
	if (fHalfThicknessPrimFoil==0.) return;
	
	G4VSolid* primFoilSolid = new G4Tubs("PrimFoilSolid", 
		0.*cm, fRadOverall, fHalfThicknessPrimFoil, 0.*deg, 360.*deg);
	G4LogicalVolume* primFoilLog = new G4LogicalVolume(primFoilSolid, fMaterialPrimFoil, "PrimFoilLog");
	
	primFoilVisAtt = new G4VisAttributes(G4Colour(0.5,1.0,0.5));
	primFoilLog->SetVisAttributes(primFoilVisAtt);
	
	new G4PVPlacement(0, 
		G4ThreeVector(0.,0.,fPosPrimFoil + fHalfThicknessPrimFoil - halfLengthWorld),
		primFoilLog,"ScatteringFoil",worldLog,false,0);
}


void ElectronBenchmarkDetector::CreateMonitor(G4LogicalVolume* worldLog){
	G4double halfLengthWorld = fPosScorer/2.;
	G4double halfThicknessMon = (fPosMon1 - fPosMon0) /2.;
	G4VSolid* monSolid = new G4Tubs("monSolid", 0.*cm, fRadOverall, halfThicknessMon, 0.*deg, 360.*deg);
	G4LogicalVolume* monLog = new G4LogicalVolume(monSolid, G4Material::GetMaterial("G4_MYLAR"), "monLog");
	
	monVisAtt = new G4VisAttributes(G4Colour(0.5,1.0,0.5));
	monLog->SetVisAttributes(monVisAtt);
	
	new G4PVPlacement(0, 
		G4ThreeVector(0.,0., fPosMon0 + halfThicknessMon - halfLengthWorld),
		monLog,"MonitorChamber",worldLog,false,0);	
}


void ElectronBenchmarkDetector::CreateHeliumBag(G4LogicalVolume* worldLog){
	G4double halfLengthWorld = fPosScorer/2.;
// Construct cylinder of Mylar
	G4double halfThicknessBag = (fPosBag1 - fPosBag0) /2.;
	G4VSolid* bagSolid = new G4Tubs("bagSolid", 0.*cm, fRadOverall, halfThicknessBag, 0.*deg, 360.*deg);
	G4LogicalVolume* bagLog = new G4LogicalVolume(bagSolid, G4Material::GetMaterial("G4_MYLAR"), "bagLog");
	
	bagVisAtt = new G4VisAttributes(G4Colour(0.5,1.0,0.5));
	bagLog->SetVisAttributes(bagVisAtt);
	
	new G4PVPlacement(0, 
		G4ThreeVector(0.,0., fPosBag0 + halfThicknessBag - halfLengthWorld),
		bagLog,"HeliumBag",worldLog,false,0);	

// Insert cylinder of Helium into the Cylinder of Mylar
	G4double halfThicknessHelium = (fPosHelium1 - fPosHelium0) /2.;
	G4VSolid* heliumSolid = new G4Tubs("heliumSolid", 0.*cm, fRadOverall, halfThicknessHelium, 0.*deg, 360.*deg);
	G4LogicalVolume* heliumLog = new G4LogicalVolume(heliumSolid, G4Material::GetMaterial("G4_He"), "heliumLog");
	
	heliumVisAtt = new G4VisAttributes(G4Colour(0.5,1.0,0.5));
	heliumLog->SetVisAttributes(heliumVisAtt);
	
	new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),heliumLog,"Helium",bagLog,false,0);	

// Insert two rings of Aluminum into the Cylinder of Helium
	G4double halfThicknessRing = fThicknessRing /2.;
	G4VSolid* ringSolid = new G4Tubs("ringSolid", fRadRingInner, fRadOverall, halfThicknessRing, 0.*deg, 360.*deg);
	G4LogicalVolume* ring0Log = new G4LogicalVolume(ringSolid, G4Material::GetMaterial("G4_Al"), "ring0Log");
	G4LogicalVolume* ring1Log = new G4LogicalVolume(ringSolid, G4Material::GetMaterial("G4_Al"), "ring1Log");
	
	ringVisAtt = new G4VisAttributes(G4Colour(0.5,1.0,0.5));
	ring0Log->SetVisAttributes(ringVisAtt);
	ring1Log->SetVisAttributes(ringVisAtt);
	
	new G4PVPlacement(0, 
		G4ThreeVector(0.,0., -halfThicknessHelium + halfThicknessRing),
		ring0Log,"Ring0",heliumLog,false,0);	
	
	new G4PVPlacement(0, 
		G4ThreeVector(0.,0., halfThicknessHelium - halfThicknessRing),
		ring1Log,"Ring1",heliumLog,false,0);	
}


void ElectronBenchmarkDetector::CreateScorer(G4LogicalVolume* worldLog){
	G4double halfLengthWorld = fPosScorer/2.;
	G4double halfThicknessScorer = fThicknessScorer /2.;
	
	G4VSolid* scorerSolid = new G4Tubs("scorerSolid", 0.*cm, fRadOverall, halfThicknessScorer, 0.*deg, 360.*deg);
	G4LogicalVolume* scorerLog = new G4LogicalVolume(scorerSolid, G4Material::GetMaterial("G4_AIR"), "scorerLog");

	scorerVisAtt = new G4VisAttributes(G4Colour(0.5,1.0,0.5));
	scorerLog->SetVisAttributes(scorerVisAtt);
	new G4PVPlacement(0,
		G4ThreeVector(0.,0., halfLengthWorld - halfThicknessScorer),
		scorerLog,"Scorer",worldLog,false,0);	
	
	G4VSolid* scorerRingSolid = new G4Tubs("scorerRingSolid", 0.*cm, fRadOverall, halfThicknessScorer, 0.*deg, 360.*deg);
	fScorerRingLog = new G4LogicalVolume(scorerRingSolid, G4Material::GetMaterial("G4_AIR"), "scorerRingLog");
	new G4PVReplica("ScorerRing",fScorerRingLog,scorerLog,kRho,G4int(fRadOverall/fWidthScorerRing),fWidthScorerRing); 
}


void ElectronBenchmarkDetector::ActivateScorer()
{
	G4String detName = "MyDetector";
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4VSensitiveDetector* sdet = SDman->FindSensitiveDetector(detName,false);
  
	if (sdet)
	{
		G4cout << "Already have the detector" << G4endl;
		fScorerRingLog->SetSensitiveDetector(sdet);
	}
	else
	{
  		G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

		G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector(detName);
	
		G4String fltName,particleName;
		G4SDParticleFilter* electronFilter = 
		new G4SDParticleFilter(fltName="electronFilter", particleName="e-");
	
		G4VPrimitiveScorer* primitive;
	
		primitive = new G4PSCellFlux("cell flux");
		det->RegisterPrimitive(primitive);
	
		primitive = new G4PSCellFlux("e cell flux");
		primitive->SetFilter(electronFilter);
		det->RegisterPrimitive(primitive);
	
		primitive = new G4PSPopulation("population");
		det->RegisterPrimitive(primitive);
	
		primitive = new G4PSPopulation("e population");
		primitive->SetFilter(electronFilter);
		det->RegisterPrimitive(primitive);
	
		G4SDManager::GetSDMpointer()->AddNewDetector(det);
		fScorerRingLog->SetSensitiveDetector(det); 
	
		G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	}
}


void ElectronBenchmarkDetector::SetPrimFoilMaterial(G4String matname){
	fMaterialPrimFoil = G4Material::GetMaterial(matname);
	UpdateGeometry();
}


void ElectronBenchmarkDetector::SetPrimFoilThickness(G4double thicknessPrimFoil){
	fHalfThicknessPrimFoil = thicknessPrimFoil / 2.;
	UpdateGeometry();
}
