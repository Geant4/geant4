// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst12DetectorMessenger.hh,v 1.2 1999-12-15 14:54:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst12DetectorMessenger_h
#define Tst12DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst12DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst12DetectorMessenger: public G4UImessenger
{
  public:
    Tst12DetectorMessenger(Tst12DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst12DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

