//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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




















