// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EvManMessenger.hh,v 1.1 1999-01-07 16:06:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4EvManMessenger_h
#define G4EvManMessenger_h 1

#include "G4UImessenger.hh"
class G4EventManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;

class G4EvManMessenger: public G4UImessenger
{
  public:
    G4EvManMessenger(G4EventManager * fEvMan);
    ~G4EvManMessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);
  private:
    G4EventManager * fEvManager;
    G4UIdirectory* eventDirectory;
    G4UIcmdWithoutParameter* abortCmd;
    G4UIcmdWithAnInteger* verboseCmd;
};

#endif

