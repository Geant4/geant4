// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02DetectorMessenger.hh,v 1.3 2000-12-04 16:24:05 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#ifndef ExN02DetectorMessenger_h
#define ExN02DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ExN02DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class ExN02DetectorMessenger: public G4UImessenger
{
  public:
    ExN02DetectorMessenger(ExN02DetectorConstruction*);
   ~ExN02DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    ExN02DetectorConstruction* myDetector;
    
    G4UIdirectory*             mydetDir;
    G4UIcmdWithAString*        TargMatCmd;
    G4UIcmdWithAString*        ChamMatCmd;    
    G4UIcmdWithADoubleAndUnit* FieldCmd;
};

#endif

