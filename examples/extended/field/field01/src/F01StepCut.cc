// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01StepCut.cc,v 1.1 2001-03-27 16:22:23 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "F01StepCut.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

F01StepCut::F01StepCut(const G4String& aName)
  : G4VDiscreteProcess(aName),MaxChargedStep(DBL_MAX)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

F01StepCut::~F01StepCut()
{
}

F01StepCut::F01StepCut(F01StepCut& right)
    :G4VDiscreteProcess(right)
{}

void F01StepCut::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
}


