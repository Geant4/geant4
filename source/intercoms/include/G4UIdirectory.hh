// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIdirectory.hh,v 2.0 1998/07/02 17:08:07 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
//

#ifndef G4UIdirectory_H
#define G4UIdirectory_H 1

#include "G4UIcommand.hh"

class G4UIdirectory : public G4UIcommand
{
  public:
    G4UIdirectory(char * theCommandPath);
    G4UIdirectory(const char * theCommandPath);
};

#endif
