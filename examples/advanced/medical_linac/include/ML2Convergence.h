#ifndef CML2ConvergenceH
#define CML2ConvergenceH

#include "G4Step.hh"

#include "ML2SinputData.h"
#include "ML2ExpVoxels.h"

class CML2Convergence
{
public:
	CML2Convergence(G4int seed, G4int saving_in_Selected_Voxels_every_events, G4String FileExperimentalData, G4bool bCompareExp, G4int minNumberOfEvents);
	~CML2Convergence(void);
	void add(const G4Step* aStep);
	G4bool runAgain();
private:
	G4bool convergenceCriteria();

	std::vector <Svoxel> voxels;
	CML2ExpVoxels *ML2ExpVoxels;

	G4String fileExperimentalData;

	G4bool bCompareExp;
	G4int minNumberOfEvents;
};

#endif

