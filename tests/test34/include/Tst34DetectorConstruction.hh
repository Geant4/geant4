#ifndef Tst34DetectorConstruction_h
#define Tst34DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"


class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class GFlashHomoShowerParamterisation;
class GFlashHitMaker;
class GFlashShowerModel;
class GFlashParticleBounds;


class Tst34DetectorConstruction : public G4VUserDetectorConstruction
{
	public:
	Tst34DetectorConstruction();
	~Tst34DetectorConstruction();
	G4VPhysicalVolume* Construct();
	
	
	private:
	G4LogicalVolume* m_experimentalHall_log;
	G4LogicalVolume* m_calo_log;
	
	G4VPhysicalVolume* m_experimentalHall_phys;  	
	G4VPhysicalVolume* m_calo_phys;
	
	G4Box *m_experimentalHall_box;
	
	G4double m_experimentalHall_x;
	G4double m_experimentalHall_y;
	G4double m_experimentalHall_z;
	
	G4double m_calo_xside;
	G4double m_calo_yside;
	G4double m_calo_zside;	

	// Gflash members	
	GFlashHomoShowerParamterisation *m_theParametrisation;
	GFlashHitMaker *m_theHMaker;
	GFlashParticleBounds *m_theParticleBounds;
	GFlashShowerModel* m_theFastShowerModel;     
	
};

#endif




















