// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst16DetectorMessenger.hh,v 1.1 1999-11-18 14:58:17 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst16DetectorMessenger_h
#define Tst16DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst16DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst16DetectorMessenger: public G4UImessenger
{
  public:
    Tst16DetectorMessenger(Tst16DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst16DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

