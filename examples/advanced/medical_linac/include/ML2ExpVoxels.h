#ifndef CML2ReadOutGeometryVoxelsH
#define CML2ReadOutGeometryVoxelsH

#include "G4ThreeVector.hh"
#include "ML2SinputData.h"
#include "G4Step.hh"

class CML2ExpVoxels 
{
public:
	CML2ExpVoxels(G4bool bHasExperimentalData, G4int saving_in_Selected_Voxels_every_events, G4int seed, G4String FileExperimentalData);
	~CML2ExpVoxels(void);
	void add(G4ThreeVector pos, G4double depEnergy, G4double density);

	inline std::vector <Svoxel> getVoxels(){return this->voxels;}

	G4int getMinNumberOfEvents();

	G4bool loadData();
private:
	void saveHeader(G4String fullOutFileName);
	void saveResults(G4String fullOutFileName, std::vector <Svoxel> voxels);
	void calculateNormalizedEd(std::vector <Svoxel> &voxels);
	std::vector <Svoxel> voxels;
	G4ThreeVector minZone, maxZone;
	G4int nCurves;
	G4int *startCurve, *stopCurve;
	G4double *chi2Factor;
	G4String headerText1, headerText2, fullFileIn, fullFileOut;
	G4String seedName, loopName;
	SGeneralData *generalData;
	G4int nParticle;
	G4int nTotalEvents, saving_in_Selected_Voxels_every_events;
	G4bool bHasExperimentalData;
};

#endif

