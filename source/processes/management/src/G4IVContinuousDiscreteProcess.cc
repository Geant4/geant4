// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IVContinuousDiscreteProcess.cc,v 1.1 1999-01-07 16:13:57 gunter Exp $
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

#include "G4IVContinuousDiscreteProcess.hh"
G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess()
                   :G4VProcess("No Name Discrete Process") 
{
  G4Exception("G4IVContinuousDiscreteProcess:: default constructor is called");
}

G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess(const G4String& aName , G4ProcessType aType)
                  : G4VProcess(aName, aType)
{
}

G4IVContinuousDiscreteProcess::~G4IVContinuousDiscreteProcess()
{
}

G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess(G4IVContinuousDiscreteProcess& right)
                  : G4VProcess(right)
{
}











