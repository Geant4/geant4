//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....#ifndef 

#ifndef myDetectorMessenger_h
#define myDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class myDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myDetectorMessenger: public G4UImessenger
{
  public:
    myDetectorMessenger(myDetectorConstruction* );
   ~myDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);

 private:
    myDetectorConstruction*    Detector;
    G4UIdirectory*             detDir;
    G4UIcmdWithAString*        SiMaterCmd;
    G4UIcmdWithAString*        SamMaterCmd;
    G4UIcmdWithAString*        HPGeMaterCmd;
    G4UIcmdWithAString*        Dia1MaterCmd;
    G4UIcmdWithAString*        Dia2MaterCmd;
    G4UIcmdWithAString*        Dia3MaterCmd;

    G4UIcmdWithADoubleAndUnit* SiThickCmd;
    G4UIcmdWithADoubleAndUnit* SamThickCmd;
    G4UIcmdWithADoubleAndUnit* HPGeThickCmd;
    G4UIcmdWithADoubleAndUnit* Dia1ThickCmd;
    G4UIcmdWithADoubleAndUnit* Dia2ThickCmd;
    G4UIcmdWithADoubleAndUnit* Dia3ThickCmd;

    G4UIcmdWithADoubleAndUnit* SiSizeYZCmd;
    G4UIcmdWithADoubleAndUnit* SamSizeYZCmd;
    G4UIcmdWithADoubleAndUnit* HPGeSizeYZCmd;
    G4UIcmdWithADoubleAndUnit* Dia1SizeYZCmd;
    G4UIcmdWithADoubleAndUnit* Dia2SizeYZCmd;
    G4UIcmdWithADoubleAndUnit* Dia3SizeYZCmd;

    G4UIcmdWithADoubleAndUnit* SenDistCmd;
    G4UIcmdWithADoubleAndUnit* DiaDistCmd;

    G4UIcmdWithADoubleAndUnit* SiAngleCmd;
    G4UIcmdWithADoubleAndUnit* HPGeAngleCmd;
    G4UIcmdWithADoubleAndUnit* Dia1AngleCmd;
    G4UIcmdWithADoubleAndUnit* Dia2AngleCmd;
    G4UIcmdWithADoubleAndUnit* Dia3AngleCmd;

    G4UIcmdWithADoubleAndUnit* SiRotCmd;
    G4UIcmdWithADoubleAndUnit* HPGeRotCmd;
    G4UIcmdWithADoubleAndUnit* Dia1RotCmd;
    G4UIcmdWithADoubleAndUnit* Dia2RotCmd;
    G4UIcmdWithADoubleAndUnit* Dia3RotCmd;

    G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif
