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
// The code was written by :
//	^Claudio Andenna claudio.andenna@iss.infn.it, claudio.andenna@ispesl.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//
// ^ISPESL and INFN Roma, gruppo collegato Sanità, Italy
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
	void create(SInputData *inputData);
	static CML2WorldConstruction* GetInstance(void);
	G4int getNParticleBackScattered(){return this->backScatteredPlane->getCML2SensDetNParticle();};
	G4int getNParticlePhaseSpace(){return this->phaseSpace->getCML2SensDetNParticle();};
	inline G4int getTotalNumberOfEventsInPhantom(){return this->phantomEnv->getTotalNumberOfEvents();};
	inline bool getBContinueRun(){return this->phaseSpace->getBContinueRun();};
private:
	void checkVolumeOverlap();
	static CML2WorldConstruction * instance;

	CML2AcceleratorConstruction *acceleratorEnv;
	CML2PhantomConstruction *phantomEnv;
	G4VPhysicalVolume* PVWorld;
	CML2PhaseSpaces *phaseSpace;
	CML2PhaseSpaces *backScatteredPlane;
};

#endif

