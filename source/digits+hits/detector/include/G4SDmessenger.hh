// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SDmessenger.hh,v 1.1 1999/01/07 16:06:25 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#ifndef G4SDmessenger_h
#define G4SDmessenger_h 1

#include "G4UImessenger.hh"

class G4SDManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class G4SDmessenger: public G4UImessenger
{
  public:
    G4SDmessenger(G4SDManager * SDManager);
    ~G4SDmessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
  
  private:
    G4SDManager * fSDMan;
    G4UIdirectory* hitsDir;
    G4UIcmdWithoutParameter* listCmd;
    G4UIcmdWithAString* activeCmd;
    G4UIcmdWithAString* inactiveCmd;
    G4UIcmdWithAnInteger* verboseCmd;
};




#endif

