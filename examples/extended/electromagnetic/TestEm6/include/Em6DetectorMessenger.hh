//  Em6DetectorMessenger.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em6DetectorMessenger_h
#define Em6DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em6DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWith3Vector;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6DetectorMessenger: public G4UImessenger
{
  public:
    Em6DetectorMessenger(Em6DetectorConstruction* );
   ~Em6DetectorMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:
    Em6DetectorConstruction* Em6Detector;

    G4UIdirectory*             Em6detDir;
    G4UIcmdWithAString*        MaterCmd;
    G4UIcmdWith3Vector*        LBinCmd;
    G4UIcmdWithADoubleAndUnit* FieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

