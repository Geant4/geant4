#ifndef CML2SteppingActionH
#define CML2SteppingActionH

#include "G4UserSteppingAction.hh"
#include "ML2Convergence.h"

class CML2ReadOutGeometryVoxels;

class CML2SteppingAction : public G4UserSteppingAction
{
public:
	CML2SteppingAction(CML2Convergence *convergence);
	~CML2SteppingAction(void);
	void UserSteppingAction(const G4Step* aStep);
private:
	CML2Convergence *convergence;
};

#endif

