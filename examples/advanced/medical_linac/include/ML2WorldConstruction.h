#ifndef CML2WorldConstructionH
#define CML2WorldConstructionH

#include "ML2SinputData.h"

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "ML2AcceleratorConstruction.h"
#include "ML2PhantomConstruction.h"
#include "ML2PhaseSpaces.h"


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

