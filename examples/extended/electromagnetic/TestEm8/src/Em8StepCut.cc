// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8StepCut.cc,v 1.1 2000-01-07 14:50:47 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//	05 March 1999 L. Urban
// --------------------------------------------------------------
// --------------------------------------------------------------

#include "Em8StepCut.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

Em8StepCut::Em8StepCut(const G4String& aName)
  : G4VDiscreteProcess(aName),MaxChargedStep(DBL_MAX)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
   }
}

Em8StepCut::~Em8StepCut()
{
}

Em8StepCut::Em8StepCut(Em8StepCut& right)
    :G4VDiscreteProcess(right)
{}

void Em8StepCut::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
}

