// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelDetectorMessenger.hh,v 1.1 2000-10-05 09:47:58 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelDetectorMessenger_h
#define GammaRayTelDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class GammaRayTelDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelDetectorMessenger: public G4UImessenger
{
public:
  GammaRayTelDetectorMessenger(GammaRayTelDetectorConstruction* );
  ~GammaRayTelDetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  GammaRayTelDetectorConstruction* GammaRayTelDetector;
  
  G4UIdirectory*             GammaRayTeldetDir;
  
  // Converter
  
  G4UIcmdWithAString*        ConverterMaterCmd;
  G4UIcmdWithADoubleAndUnit* ConverterThickCmd;
  
  // Silicon Tile
  
  G4UIcmdWithADoubleAndUnit* SiliconThickCmd;
  G4UIcmdWithADoubleAndUnit* SiliconTileXYCmd;
  G4UIcmdWithAnInteger*      NbSiTilesCmd;    
  G4UIcmdWithADoubleAndUnit* SiliconPitchCmd;

  // Tracker
  
  G4UIcmdWithAnInteger*      NbTKRLayersCmd;    
  G4UIcmdWithADoubleAndUnit* LayerDistanceCmd;
  G4UIcmdWithADoubleAndUnit* ViewsDistanceCmd;
  
  // Calorimeter 
  
  G4UIcmdWithADoubleAndUnit* CALThickCmd;
  G4UIcmdWithAnInteger*      NbCALBarsCmd;
  G4UIcmdWithAnInteger*      NbCALLayersCmd;    

  // Anticoincidence

  G4UIcmdWithADoubleAndUnit* ACDThickCmd;

  // Total

  G4UIcmdWithADoubleAndUnit* MagFieldCmd;
  G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif










