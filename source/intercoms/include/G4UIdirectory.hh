// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIdirectory.hh,v 1.2 1999-10-29 06:06:45 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4UIdirectory_H
#define G4UIdirectory_H 1

#include "G4UIcommand.hh"

// class description:
//  A concrete class of G4UIcommand. This class defines a command
// directory which can have commands.
//  General information of G4UIcommand is given in G4UIcommand.hh.

class G4UIdirectory : public G4UIcommand
{
  public: // with description
    G4UIdirectory(char * theCommandPath);
    G4UIdirectory(const char * theCommandPath);
    // Constructors. The argument is a full path directory which
    // starts and ends with "/".
};

#endif
