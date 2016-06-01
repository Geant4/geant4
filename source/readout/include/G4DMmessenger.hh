// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DMmessenger.hh,v 2.1 1998/07/12 03:08:24 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4DMmessenger_h
#define G4DMmessenger_h 1

#include "G4UImessenger.hh"

class G4DigiManager;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

class G4DMmessenger: public G4UImessenger
{
  public:
    G4DMmessenger(G4DigiManager * DigiManager);
    ~G4DMmessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
  
  private:
    G4DigiManager * fDMan;
    G4UIdirectory* digiDir;
    G4UIcmdWithoutParameter* listCmd;
    G4UIcmdWithAString* digiCmd;
    G4UIcmdWithAnInteger* verboseCmd;
};




#endif

