// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17SteppingMessenger.hh,v 1.2 1999-12-15 14:54:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
#ifndef Tst17SteppingMessenger_h
#define Tst17SteppingMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

class Tst17SteppingAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

class Tst17SteppingMessenger: public G4UImessenger
{
  public:

   Tst17SteppingMessenger(Tst17SteppingAction* );
   ~Tst17SteppingMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   Tst17SteppingAction* steppingAction;

   G4UIdirectory* steppingDir;

};

#endif

