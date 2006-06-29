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
#ifndef ExGflashDetectorConstruction_h
#define ExGflashDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"


class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class GFlashHomoShowerParameterisation;
class GFlashHitMaker;
class GFlashShowerModel;
class GFlashParticleBounds;


class ExGflashDetectorConstruction : public G4VUserDetectorConstruction
{
	public:
	ExGflashDetectorConstruction();
	~ExGflashDetectorConstruction();
	G4VPhysicalVolume* Construct();
	
	const G4VPhysicalVolume* GetCristal(int num__crystal)
                                {return m_crystal_phys[num__crystal];};


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
	
	G4int    m_NbOfCrystals;            	// Nb of chambers in the tracker region
	G4double m_CrystalWidth;            	// width of the chambers
	G4double m_CrystalWidht;
	G4double m_CrystalLenght;
	//@@@  ExGflashDetectorConstruction : wie mache ich das am besten ?
	G4Box *m_crystal[100];      
	G4LogicalVolume* m_crystal_log[100];
	G4VPhysicalVolume*  m_crystal_phys[100];
	
	// Gflash members	
	GFlashHomoShowerParameterisation *m_theParameterisation;
	GFlashHitMaker *m_theHMaker;
	GFlashParticleBounds *m_theParticleBounds;
	GFlashShowerModel* m_theFastShowerModel;     
        G4Region* aRegion;
};

#endif




















