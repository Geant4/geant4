// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F03StepCut.cc,v 1.1 2001-06-08 11:55:58 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "F03StepCut.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

F03StepCut::F03StepCut(const G4String& aName)
  : G4VDiscreteProcess(aName),MaxChargedStep(DBL_MAX)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

F03StepCut::~F03StepCut()
{
}

F03StepCut::F03StepCut(F03StepCut& right)
    :G4VDiscreteProcess(right)
{}

void F03StepCut::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
}


