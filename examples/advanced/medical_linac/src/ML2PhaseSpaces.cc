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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2PhaseSpaces.hh"

CML2PhaseSpaces::CML2PhaseSpaces()
{}

CML2PhaseSpaces::~CML2PhaseSpaces(void)
{}

/* NOT implemented at the moment
bool CML2PhaseSpaces::createPlane(G4VPhysicalVolume  *, G4String , G4ThreeVector , G4ThreeVector )
{
 G4cout << " Not implemented at the moment" << G4endl;

	// constructor for killer plane
	bool bCreated = false;
	G4Material *Vacum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Box *box;
	G4LogicalVolume *logVol;
	box = new G4Box("KBox", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, Vacum, name+"KLV", 0, 0, 0);
	phVol = new G4PVPlacement(0, centre, name+"KPV", logVol, PVWorld, false, 0);

	G4VisAttributes* simplePhSpVisAtt = new G4VisAttributes(G4Colour::Cyan());
	simplePhSpVisAtt->SetVisibility(true);
	simplePhSpVisAtt->SetForceSolid(true);
	logVol->SetVisAttributes(simplePhSpVisAtt);

	bCreated = true;
	return bCreated;
}

bool CML2PhaseSpaces::createPlane(G4int, //idSD_Type, 
G4int,//max_N_particles_in_PhSp_File,
G4int, //seed,
G4int,// nMaxParticlesInRamPhaseSpace, 
G4VPhysicalVolume*,//  *PVWorld, 
G4String, //name, G4String PhaseSpaceOutFile, G4bool bSavePhaseSpace, G4bool bStopAtPhaseSpace, G4ThreeVector centre, G4ThreeVector halfSize, SPrimaryParticle *primaryParticleData, G4double  accTargetZPosition)
{
	// constructor for phase space plane
	bool bCreated = false;
	G4Material *Vacum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Box *box;
	G4LogicalVolume *logVol;
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, Vacum, name+"LV", 0, 0, 0);
	phVol= new G4PVPlacement(0, centre, name+"PV", logVol, PVWorld, false, 0);

	G4VisAttributes* simplePhSpVisAtt = new G4VisAttributes(G4Colour::Yellow());
	simplePhSpVisAtt -> SetVisibility(true);
	simplePhSpVisAtt -> SetForceSolid(true);
	logVol -> SetVisAttributes(simplePhSpVisAtt);

	sensDetParticle = new CML2SDWithParticle(idSD_Type, max_N_particles_in_PhSp_File, seed, nMaxParticlesInRamPhaseSpace, name, PhaseSpaceOutFile, bSavePhaseSpace, bStopAtPhaseSpace, primaryParticleData, accTargetZPosition);
	G4SDManager *SDManager = G4SDManager::GetSDMpointer();
	SDManager -> AddNewDetector(sensDetParticle);
	logVol -> SetSensitiveDetector(sensDetParticle);
	bCreated = true;
	return bCreated;

}
*/
