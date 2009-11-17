#ifndef CML2InputDataH
#define CML2InputDataH

#include "G4UImessenger.hh"
#include "ML2MainMessenger.h"

class CML2MainMessenger;

class CML2CInputData
{
public:
	CML2CInputData(void);
	~CML2CInputData(void);

	inline void setbOnlyVisio(G4bool val){this->bOnlyVisio=val;};
	inline void setPhaseSpaceCentre(G4ThreeVector val){this->inputData.generalData.centrePhaseSpace.set(val.getX(), val.getY(), val.getZ());};
	inline void setPhaseSpaceHalfSize(G4ThreeVector val){this->inputData.generalData.halfSizePhaseSpace.set(val.getX(), val.getY(), val.getZ());};
	inline void setbSavePhaseSPace(G4bool val){this->inputData.generalData.bSavePhaseSpace=val;};
	inline void setbStopAtPhaseSpace(G4bool val){this->inputData.generalData.bStopAtPhaseSpace=val;};
	inline void setPhaseSpaceOutFile(G4String val){this->inputData.generalData.PhaseSpaceOutFile=val;};

	inline void setbSaveROG(G4bool val){this->inputData.generalData.bSaveROG=val;};
	inline void setROGOutFile(G4String val){this->inputData.generalData.ROGOutFile=val;};

	inline void setMinNumberOfEvents(G4double val){this->inputData.generalData.minNumberOfEvents=val;};

	inline void setBCompareExp(G4bool val){this->inputData.generalData.bCompareExp=val;};
	inline void setFileExperimentalData(G4String val){this->inputData.generalData.fileExperimentalData=val;};
	inline void setNBeams(G4int val){this->inputData.generalData.nBeam=val;};
	inline void setNMaxParticlesInRamPlanePhaseSpace(G4int val){this->inputData.generalData.nMaxParticlesInRamPlanePhaseSpace=val;};

	inline void setSaving_in_Selected_Voxels_every_events(G4int val){this->inputData.generalData.saving_in_Selected_Voxels_every_events=val;};
	inline void setSaving_in_ROG_Voxels_every_events(G4int val){this->inputData.generalData.saving_in_ROG_Voxels_every_events=val;};
	inline void setMax_N_particles_in_PhSp_File(G4int val){this->inputData.generalData.max_N_particles_in_PhSp_File=val;};

	G4bool bOnlyVisio;
	SInputData inputData;
private:
	CML2MainMessenger *ML2MainMessenger;
};


#endif

