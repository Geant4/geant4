// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVADetectorMessenger.hh,v 1.3 2001-02-07 17:30:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#ifndef TstVADetectorMessenger_h
#define TstVADetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class TstVADetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class TstVADetectorMessenger: public G4UImessenger
{
  public:
    TstVADetectorMessenger(TstVADetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    TstVADetectorConstruction* myDetector;
    G4UIdirectory*             mydetDir;
    G4UIcmdWithAString*        selDetCmd;
    G4UIcmdWithAString*        switchCmd;
    G4UIcmdWithAString*        selMatCmd;

};

#endif

