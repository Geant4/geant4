// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14StepCut.cc,v 1.1 1999-05-29 14:12:13 stesting Exp $
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

#include "Tst14StepCut.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

Tst14StepCut::Tst14StepCut(const G4String& aName)
  : G4VDiscreteProcess(aName),MaxChargedStep(DBL_MAX)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
   }
}

Tst14StepCut::~Tst14StepCut()
{
}

Tst14StepCut::Tst14StepCut(Tst14StepCut& right)
    :G4VDiscreteProcess(right)
{}

void Tst14StepCut::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
}

