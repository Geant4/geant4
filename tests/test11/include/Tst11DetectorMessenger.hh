// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst11DetectorMessenger.hh,v 1.1 1999-01-08 16:35:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst11DetectorMessenger_h
#define Tst11DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst11DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst11DetectorMessenger: public G4UImessenger
{
  public:
    Tst11DetectorMessenger(Tst11DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst11DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

