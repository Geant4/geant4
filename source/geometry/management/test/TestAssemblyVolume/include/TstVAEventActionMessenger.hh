// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVAEventActionMessenger.hh,v 1.3 2001-02-07 17:30:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#ifndef TstVAEventActionMessenger_h 
#define TstVAEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class TstVAEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class TstVAEventActionMessenger: public G4UImessenger
{
  public:
    TstVAEventActionMessenger(TstVAEventAction*);
   ~TstVAEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    TstVAEventAction*     eventAction;   
    G4UIcmdWithAString*   DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
