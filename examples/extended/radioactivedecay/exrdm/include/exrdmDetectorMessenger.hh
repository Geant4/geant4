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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef exrdmDetectorMessenger_h
#define exrdmDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class exrdmDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class exrdmDetectorMessenger: public G4UImessenger
{
  public:
    exrdmDetectorMessenger(exrdmDetectorConstruction*);
   ~exrdmDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    exrdmDetectorConstruction* myDetector;
    
    G4UIdirectory*             exrdmDir;
    G4UIdirectory*             detDir;
    G4UIcmdWithAString*        TargMatCmd;
    G4UIcmdWithAString*        DetectMatCmd;
    G4UIcmdWithADoubleAndUnit* TargRadiusCmd;
    G4UIcmdWithADoubleAndUnit* DetectThicknessCmd;
    G4UIcmdWithADoubleAndUnit* TargLengthCmd;
    G4UIcmdWithADoubleAndUnit* DetectLengthCmd;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

