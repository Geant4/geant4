// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5DetectorMessenger.hh,v 1.1 1999-10-12 12:23:30 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5DetectorMessenger_h
#define Em5DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em5DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5DetectorMessenger: public G4UImessenger
{
  public:
    Em5DetectorMessenger(Em5DetectorConstruction* );
   ~Em5DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em5DetectorConstruction*   Em5Detector;
    
    G4UIdirectory*             Em5detDir;

    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithADoubleAndUnit* AbsThickCmd;
    G4UIcmdWithADoubleAndUnit* AbsSizYZCmd;

    G4UIcmdWithADoubleAndUnit* AbsXposCmd;

    G4UIcmdWithAString*        WorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* WorldXCmd;
    G4UIcmdWithADoubleAndUnit* WorldYZCmd;

    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;

};

#endif

