// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08DetectorMessenger.hh,v 1.1 1999-01-08 16:35:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef T08DetectorMessenger_h
#define T08DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class T08DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class T08DetectorMessenger: public G4UImessenger
{
  public:
    T08DetectorMessenger(T08DetectorConstruction * myDet);
    ~T08DetectorMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    T08DetectorConstruction* myDetector;
    
    G4UIdirectory*             mydetDir;
    G4UIcmdWithAString*        MatCmd;
    G4UIcmdWithADoubleAndUnit* SizeCmd;
    G4UIcmdWithADoubleAndUnit* FieldCmd;
};

#endif

