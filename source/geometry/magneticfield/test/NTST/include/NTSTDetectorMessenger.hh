// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTDetectorMessenger.hh,v 1.1 2003-11-07 21:30:28 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTDetectorMessenger_h
#define NTSTDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class NTSTDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTDetectorMessenger: public G4UImessenger
{
  public:
    NTSTDetectorMessenger(NTSTDetectorConstruction* );
   ~NTSTDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    NTSTDetectorConstruction*   NTSTDetector;
    
    G4UIdirectory*             NTSTdetDir;
    G4UIcmdWithAnInteger*      DebugCmd;
    G4UIcmdWithAnInteger*      NSubLayer;
    G4UIcmdWithADoubleAndUnit* MotherOuterRadius;
    G4UIcmdWithAString*        InputFileNameCmd;
    G4UIcmdWithAString*        DisableDet;
    G4UIcmdWithoutParameter*   fieldStat;
};

#endif

