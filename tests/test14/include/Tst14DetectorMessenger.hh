// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14DetectorMessenger.hh,v 1.4 1999-12-15 14:54:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst14DetectorMessenger_h
#define Tst14DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst14DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst14DetectorMessenger: public G4UImessenger
{
  public:
    Tst14DetectorMessenger(Tst14DetectorConstruction* );
   ~Tst14DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Tst14DetectorConstruction*   Tst14Detector;
    
    G4UIdirectory*             Tst14detDir;

    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithADoubleAndUnit* AbsThickCmd;
    G4UIcmdWithADoubleAndUnit* AbsRadCmd;

    G4UIcmdWithADoubleAndUnit* AbsZposCmd;

    G4UIcmdWithAString*        WorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* WorldZCmd;
    G4UIcmdWithADoubleAndUnit* WorldRCmd;

    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;

};

#endif

