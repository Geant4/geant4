#ifndef hTestDetectorMessenger_h
#define hTestDetectorMessenger_h 1

// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestDetectorMessenger -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4UImessenger.hh"

class hTestDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestDetectorMessenger: public G4UImessenger
{
public: // Without description

    hTestDetectorMessenger(hTestDetectorConstruction* );
   ~hTestDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:

    hTestDetectorConstruction* hDet;
    
    G4UIdirectory*             hTestdetDir;
    G4UIdirectory*             hTestdetDir1;
    G4UIdirectory*             hTestdetDir2;

    G4UIcmdWithAnInteger*      NumOfAbsCmd;
    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithADoubleAndUnit* AbsThickCmd;
    G4UIcmdWithADoubleAndUnit* AbsSizYZCmd;
    G4UIcmdWithAString*        WorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* WorldXCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
    G4UIcmdWithADoubleAndUnit* XMagFieldCmd;
    G4UIcmdWithADoubleAndUnit* YMagFieldCmd;
    G4UIcmdWithADoubleAndUnit* ZMagFieldCmd;
    G4UIcmdWithAString*        HistoCmd;
    G4UIcmdWithAnInteger*      NumOfEvt;
    G4UIcmdWithAnInteger*      verbCmd;
    G4UIcmdWithAnInteger*      intCmd;
    G4UIcmdWithAnInteger*      nhistCmd;
    G4UIcmdWithAnInteger*      nDebugSCmd;
    G4UIcmdWithAnInteger*      nDebugECmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif





