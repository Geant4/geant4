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


#ifndef CML2PhaseSpacesH
#define CML2PhaseSpacesH

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "ML2SinputData.hh"

class CML2PhaseSpaces
{
public:
	CML2PhaseSpaces();
	~CML2PhaseSpaces(void);
//	bool createPlane(G4VPhysicalVolume  *PVWorld, G4String name, G4ThreeVector centre, G4ThreeVector halfSize);
//	bool createPlane(G4int idSD_Type, G4int max_N_particles_in_PhSp_File, G4int seed, G4int nMaxParticlesInRamPhaseSpace, G4VPhysicalVolume  *PVWorld, G4String name, G4String PhaseSpaceOutFile, G4bool bSavePhaseSpace, G4bool bStopAtPhaseSpace, G4ThreeVector centre, G4ThreeVector halfSize, SPrimaryParticle *primaryParticleData, G4double  accTargetZPosition);
	//G4int getCML2SensDetNParticle(){return sensDetParticle->getTotalNumberOfParticles();}
//	inline CML2SDWithParticle* getCML2SensitiveDetectorParticle(){return sensDetParticle->getCML2SensitiveDetectorParticle();}
//	inline void save(){sensDetParticle->save();}
private:

  //      G4VPhysicalVolume *phVol;

};


#endif
