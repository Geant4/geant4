// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst02DetectorMessenger.hh,v 1.1 1999-01-08 16:34:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst02DetectorMessenger_h
#define Tst02DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst02DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst02DetectorMessenger: public G4UImessenger
{
  public:
    Tst02DetectorMessenger(Tst02DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst02DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

