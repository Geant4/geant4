// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Run.cc,v 1.2 1999-12-15 14:53:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4Run.hh"

G4Allocator<G4Run> aRunAllocator;

G4Run::G4Run()
:runID(0),numberOfEvent(0),HCtable(NULL),DCtable(NULL)
{;}

G4Run::~G4Run()
{;}

