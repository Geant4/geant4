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


#ifndef CML2WorldConstructionH
#define CML2WorldConstructionH

#include "ML2SinputData.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "ML2AcceleratorConstruction.hh"
#include "ML2PhantomConstruction.hh"
#include "ML2PhaseSpaces.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RegionStore.hh"

#include "ML2PrimaryGenerationAction.hh"

class G4VPhysicalVolume;
class CML2PhantomConstruction;
class CML2AcceleratorConstruction;
class CML2PhaseSpaces;

class CML2WorldConstruction : public G4VUserDetectorConstruction
{
public:
	CML2WorldConstruction(void);
	~CML2WorldConstruction(void);
	G4VPhysicalVolume* Construct();
	bool create(SInputData *inputData, bool bOnlyVisio);
	static CML2WorldConstruction* GetInstance(void);
	G4int getNParticleBackScattered(){return backScatteredPlane->getCML2SensDetNParticle();}
	G4int getNParticlePhaseSpace(){return phaseSpace->getCML2SensDetNParticle();}
	inline G4int getTotalNumberOfEventsInPhantom(){return phantomEnv->getTotalNumberOfEvents();}
	inline CML2AcceleratorConstruction * getCML2AcceleratorConstruction(){return acceleratorEnv;}

	bool newGeometry();
	bool getWorldCreated(){return bWorldCreated;}
	
	inline void savePhantomData(){if (phantomEnv!=0){phantomEnv->saveData();}}
	inline void savePhaseSpaceData(){if (phaseSpace!=0){phaseSpace->save();}}
	inline CML2PhantomConstruction * getPhantomWorld(){return phantomEnv;}
	inline CML2AcceleratorConstruction * getAcceleratorWorld(){return acceleratorEnv;}
	void checkVolumeOverlap();

	inline G4bool getbOnlyVisio(){return bOnlyVisio;}
private:
	static CML2WorldConstruction * instance;

	bool bWorldCreated;
	CML2AcceleratorConstruction *acceleratorEnv;
	CML2PhantomConstruction *phantomEnv;
	G4VPhysicalVolume* PVWorld;
	CML2PhaseSpaces *phaseSpace;
	CML2PhaseSpaces *backScatteredPlane;
	G4bool bOnlyVisio;
};

#endif

