// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VEvtBiasMechanism.cc,v 1.2 1999-12-15 14:53:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//	
//	
// ------------------------------------------------------------
//   Implemented for the new scheme                 13 Apr. 1999  H.Kurahige

#include "G4VEvtBiasMechanism.hh"

G4VEvtBiasMechanism::G4VEvtBiasMechanism(const G4String& name):
  theEBName(name)
{
}

G4VEvtBiasMechanism::G4VEvtBiasMechanism(const G4VEvtBiasMechanism& right)
{
  theEBName = right.theEBName;
}

G4VEvtBiasMechanism::~G4VEvtBiasMechanism()
{
}













