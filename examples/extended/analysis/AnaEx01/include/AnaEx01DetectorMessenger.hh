// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01DetectorMessenger.hh,v 1.1.1.1 2000-09-14 11:37:21 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef AnaEx01DetectorMessenger_h
#define AnaEx01DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class AnaEx01DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class AnaEx01DetectorMessenger: public G4UImessenger
{
  public:
    AnaEx01DetectorMessenger(AnaEx01DetectorConstruction* );
   ~AnaEx01DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    AnaEx01DetectorConstruction* AnaEx01Detector;
    
    G4UIdirectory*             AnaEx01detDir;
    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithAString*        GapMaterCmd;
    G4UIcmdWithADoubleAndUnit* AbsThickCmd;
    G4UIcmdWithADoubleAndUnit* GapThickCmd;
    G4UIcmdWithADoubleAndUnit* SizeYZCmd;
    G4UIcmdWithAnInteger*      NbLayersCmd;    
    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif

