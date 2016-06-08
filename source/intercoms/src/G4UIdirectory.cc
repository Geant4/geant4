// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIdirectory.cc,v 1.1.10.1 1999/12/07 20:49:03 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
//

#include "G4UIdirectory.hh"

G4UIdirectory::G4UIdirectory(char * theCommandPath)
:G4UIcommand(theCommandPath,NULL)
{;}

G4UIdirectory::G4UIdirectory(const char * theCommandPath)
:G4UIcommand(theCommandPath,NULL)
{;}

