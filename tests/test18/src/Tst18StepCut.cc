// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst18StepCut.cc,v 1.2 2000-06-06 11:36:19 flei Exp $
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

#include "Tst18StepCut.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

Tst18StepCut::Tst18StepCut(const G4String& aName)
  : G4VDiscreteProcess(aName),MaxChargedStep(DBL_MAX)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
   }
}

Tst18StepCut::~Tst18StepCut()
{
}

Tst18StepCut::Tst18StepCut(Tst18StepCut& right)
    :G4VDiscreteProcess(right)
{}

void Tst18StepCut::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
}

