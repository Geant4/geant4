// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02DetectorMessenger.hh,v 1.2 1999-12-15 14:49:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef ExN02DetectorMessenger_h
#define ExN02DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ExN02DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class ExN02DetectorMessenger: public G4UImessenger
{
  public:
    ExN02DetectorMessenger(ExN02DetectorConstruction * myDet);
    ~ExN02DetectorMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    ExN02DetectorConstruction* myDetector;
    
    G4UIdirectory*             mydetDir;
    G4UIcmdWithAString*        MatCmd;
    G4UIcmdWithADoubleAndUnit* SizeCmd;
    G4UIcmdWithADoubleAndUnit* FieldCmd;
};

#endif

