// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst09DetectorMessenger.hh,v 1.2 1999-12-15 14:54:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst09DetectorMessenger_h
#define Tst09DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst09DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst09DetectorMessenger: public G4UImessenger
{
  public:
    Tst09DetectorMessenger(Tst09DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst09DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

