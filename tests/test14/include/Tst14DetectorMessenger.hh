// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14DetectorMessenger.hh,v 1.2 1999-06-14 14:28:33 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst14DetectorMessenger_h
#define Tst14DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst14DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst14DetectorMessenger: public G4UImessenger
{
  public:
    Tst14DetectorMessenger(Tst14DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst14DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

