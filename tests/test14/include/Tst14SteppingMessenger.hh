// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14SteppingMessenger.hh,v 1.1 1999-05-29 14:12:07 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
#ifndef Tst14SteppingMessenger_h
#define Tst14SteppingMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

class Tst14SteppingAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

class Tst14SteppingMessenger: public G4UImessenger
{
  public:

   Tst14SteppingMessenger(Tst14SteppingAction* );
   ~Tst14SteppingMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   Tst14SteppingAction* steppingAction;

   G4UIdirectory* steppingDir;

};

#endif

