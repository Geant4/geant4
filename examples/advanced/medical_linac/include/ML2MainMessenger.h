#ifndef CML2MainMessengerH
#define CML2MainMessengerH

#include "ML2Main.h"
#include "G4UImessenger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "ML2CInputData.h"


class CML2Main;
class CML2CInputData;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAInteger;

class CML2MainMessenger : public G4UImessenger 
{
public:
	CML2MainMessenger(CML2CInputData *CInputData);
	~CML2MainMessenger(void);
	void SetNewValueOLD(G4UIcommand* cmd, G4String newValue);
	void SetNewValue(G4UIcommand* cmd, G4String newValue);

private:

	CML2CInputData *CInputData;
	CML2Main *ML2Main;

	G4UIcmdWith3VectorAndUnit *phaseSpaceCentre, *phaseSpaceHalfSize;
	G4UIcmdWithAString *phaseSPaceOutFile, *ROGOutFile;
	G4UIcmdWithABool *bSavePhaseSpace, *bStopAtPhaseSpace, *bSaveROG;

	G4UIcmdWithAnInteger *nBeam, *nMaxParticlesInRamPlanePhaseSpace, *minNumberOfEvents;
	G4UIcmdWithABool *bCompareExp, *bOnlyVisio;
	G4UIcmdWithAString * fileExperimentalData;
	G4UIcmdWithAnInteger *saving_in_Selected_Voxels_every_events;
	G4UIcmdWithAnInteger *saving_in_ROG_Voxels_every_events;
	G4UIcmdWithAnInteger *max_N_particles_in_PhSp_File; 
};

#endif
