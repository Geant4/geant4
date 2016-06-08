// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Run.cc,v 1.1.10.1 1999/12/07 20:52:59 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#include "G4Run.hh"

G4Allocator<G4Run> aRunAllocator;

G4Run::G4Run()
:runID(0),numberOfEvent(0),HCtable(NULL),DCtable(NULL)
{;}

G4Run::~G4Run()
{;}

