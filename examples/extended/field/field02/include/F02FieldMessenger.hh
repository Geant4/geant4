// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F02FieldMessenger.hh,v 1.1 2001-04-04 15:11:15 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F02FieldMessenger_h
#define F02FieldMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class F02ElectroMagneticField;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;


class F02FieldMessenger: public G4UImessenger
{
  public:
    F02FieldMessenger(F02ElectroMagneticField* );
   ~F02FieldMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    void SetNewValue(G4UIcommand*, G4int);
    
  private:

    F02ElectroMagneticField*   fEMfield;
    
    G4UIdirectory*             F02detDir;

    G4UIcmdWithAnInteger*      StepperCmd;
    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithADoubleAndUnit* MinStepCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;

    G4UIcmdWithAString*        AbsMaterCmd;


};

#endif

