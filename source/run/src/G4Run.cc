// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Run.cc,v 2.0 1998/07/02 17:28:15 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//

#include "G4Run.hh"

G4Allocator<G4Run> aRunAllocator;

G4Run::G4Run()
:runID(0),numberOfEvent(0),HCtable(NULL),DCtable(NULL)
{;}

G4Run::~G4Run()
{;}

