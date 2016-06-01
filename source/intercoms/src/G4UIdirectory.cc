// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIdirectory.cc,v 2.0 1998/07/02 17:08:55 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
//

#include "G4UIdirectory.hh"

G4UIdirectory::G4UIdirectory(char * theCommandPath)
:G4UIcommand(theCommandPath,NULL)
{;}

G4UIdirectory::G4UIdirectory(const char * theCommandPath)
:G4UIcommand(theCommandPath,NULL)
{;}

