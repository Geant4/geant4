// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em10DetectorMessenger.hh,v 1.2 2001-03-19 17:59:25 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em10DetectorMessenger_h
#define Em10DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em10DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;


class Em10DetectorMessenger: public G4UImessenger
{
  public:
    Em10DetectorMessenger(Em10DetectorConstruction* );
   ~Em10DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    void SetNewValue(G4UIcommand*, G4int);
    
  private:
    Em10DetectorConstruction*   Em10Detector;
    
    G4UIdirectory*             Em10detDir;

    G4UIcmdWithAnInteger*       ModelCmd;
    G4UIcmdWithAnInteger*       FoilNumCmd;

    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithADoubleAndUnit* AbsThickCmd;
    G4UIcmdWithADoubleAndUnit* AbsRadCmd;

    G4UIcmdWithAString*        RadiatorMaterCmd;
    G4UIcmdWithADoubleAndUnit* RadiatorThickCmd;

    G4UIcmdWithADoubleAndUnit* GasGapThickCmd;

    G4UIcmdWithADoubleAndUnit* AbsZposCmd;

    G4UIcmdWithAString*        WorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* WorldZCmd;
    G4UIcmdWithADoubleAndUnit* WorldRCmd;

    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;

};

#endif

