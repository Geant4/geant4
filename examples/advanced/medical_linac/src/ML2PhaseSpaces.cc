#include "ML2PhaseSpaces.h"
#include "ML2ReadOutGeometry.h"

CML2PhaseSpaces::CML2PhaseSpaces():sensDetParticle(0), sensDetVoxelized(0)
{}

CML2PhaseSpaces::~CML2PhaseSpaces(void)
{}

bool CML2PhaseSpaces::createPlane(G4VPhysicalVolume  *PVWorld, G4String name, G4ThreeVector centre, G4ThreeVector halfSize)
{
	// constructor for killer plane
	bool bCreated=false;
	G4Material *Vacum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VPhysicalVolume *phVol;
	box = new G4Box("KBox", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, Vacum, name+"KLV", 0, 0, 0);
	phVol= new G4PVPlacement(0, centre, name+"KPV", logVol, PVWorld, false, 0);

	G4VisAttributes* simplePhSpVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simplePhSpVisAtt->SetVisibility(true);
	simplePhSpVisAtt->SetForceSolid(true);
	logVol->SetVisAttributes(simplePhSpVisAtt);


	this->sensDetParticle=new CML2SDWithParticle();
	G4SDManager *SDManager=G4SDManager::GetSDMpointer();
	SDManager->AddNewDetector(this->sensDetParticle);
	logVol->SetSensitiveDetector(this->sensDetParticle);
	bCreated=true;
	return bCreated;
}

bool CML2PhaseSpaces::createPlane(G4int idSD_Type, G4int max_N_particles_in_PhSp_File, G4int seed, G4int nMaxParticlesInRamPhaseSpace, G4VPhysicalVolume  *PVWorld, G4String name, G4String PhaseSpaceOutFile, G4bool bSavePhaseSpace, G4bool bStopAtPhaseSpace, G4ThreeVector centre, G4ThreeVector halfSize, SPrimaryParticle *primaryParticleData)
{
	// constructor for phase space plane
	bool bCreated=false;
	G4Material *Vacum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VPhysicalVolume *phVol;
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, Vacum, name+"LV", 0, 0, 0);
	phVol= new G4PVPlacement(0, centre, name+"PV", logVol, PVWorld, false, 0);

	G4VisAttributes* simplePhSpVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simplePhSpVisAtt->SetVisibility(true);
	simplePhSpVisAtt->SetForceSolid(true);
	logVol->SetVisAttributes(simplePhSpVisAtt);

	this->sensDetParticle=new CML2SDWithParticle(idSD_Type, max_N_particles_in_PhSp_File, seed, nMaxParticlesInRamPhaseSpace, name, PhaseSpaceOutFile, bSavePhaseSpace, bStopAtPhaseSpace, primaryParticleData);
	G4SDManager *SDManager=G4SDManager::GetSDMpointer();
	SDManager->AddNewDetector(this->sensDetParticle);
	logVol->SetSensitiveDetector(this->sensDetParticle);
	bCreated=true;
	return bCreated;
}

