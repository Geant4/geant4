// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T07DetectorMessenger.hh,v 1.2 1999-12-15 14:54:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef T07DetectorMessenger_h
#define T07DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class T07DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class T07DetectorMessenger: public G4UImessenger
{
  public:
    T07DetectorMessenger(T07DetectorConstruction* );
   ~T07DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    T07DetectorConstruction* T07Detector;
    
    G4UIdirectory*             T07detDir;
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

