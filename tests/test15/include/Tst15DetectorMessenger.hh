// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst15DetectorMessenger.hh,v 1.1 1999-11-18 14:48:13 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst15DetectorMessenger_h
#define Tst15DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst15DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst15DetectorMessenger: public G4UImessenger
{
  public:
    Tst15DetectorMessenger(Tst15DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst15DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

