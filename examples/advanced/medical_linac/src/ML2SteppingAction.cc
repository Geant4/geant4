#include "ML2SteppingAction.h"
#include "globals.hh"
#include "G4Step.hh"
#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"

#include "G4VReadOutGeometry.hh"

CML2SteppingAction::CML2SteppingAction(CML2Convergence *convergence)
{
	this->convergence=convergence;
}

CML2SteppingAction::~CML2SteppingAction(void)
{}
void CML2SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	this->convergence->add(aStep);
}

