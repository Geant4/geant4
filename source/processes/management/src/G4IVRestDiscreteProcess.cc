// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IVRestDiscreteProcess.cc,v 1.4 1999-04-30 08:02:31 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// --------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
// ------------------------------------------------------------

#include "G4IVRestDiscreteProcess.hh"
G4IVRestDiscreteProcess::G4IVRestDiscreteProcess()
                   :G4VProcess("No Name Discrete Process"),BIGSTEP(1.e10)

{
  G4Exception("G4IVRestDiscreteProcess:: default constructor is called");
}

G4IVRestDiscreteProcess::G4IVRestDiscreteProcess(const G4String& aName , G4ProcessType aType)
                  : G4VProcess(aName, aType),BIGSTEP(1.e10)
{
}

G4IVRestDiscreteProcess::~G4IVRestDiscreteProcess()
{
}

G4IVRestDiscreteProcess::G4IVRestDiscreteProcess(G4IVRestDiscreteProcess& right)
                  : G4VProcess(right),BIGSTEP(1.e10)
{
}

G4double G4IVRestDiscreteProcess::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{
  G4double value = DBL_MAX ;

  return value;
}
