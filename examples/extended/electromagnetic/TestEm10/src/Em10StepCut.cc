// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em10StepCut.cc,v 1.1 2000-07-14 15:51:37 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Em10StepCut.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

Em10StepCut::Em10StepCut(const G4String& aName)
  : G4VDiscreteProcess(aName),MaxChargedStep(DBL_MAX)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

Em10StepCut::~Em10StepCut()
{
}

Em10StepCut::Em10StepCut(Em10StepCut& right)
    :G4VDiscreteProcess(right)
{}

void Em10StepCut::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
}


