// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: DetectorMessenger.hh,v 1.1 2001-05-29 19:17:32 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

