// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst06DetectorMessenger.hh,v 1.2 1999-12-15 14:54:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst06DetectorMessenger_h
#define Tst06DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst06DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst06DetectorMessenger: public G4UImessenger
{
  public:
    Tst06DetectorMessenger(Tst06DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst06DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

