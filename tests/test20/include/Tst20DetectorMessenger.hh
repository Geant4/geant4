// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20DetectorMessenger.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20DetectorMessenger  ------
//
// ************************************************************
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst20DetectorMessenger_h
#define Tst20DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst20DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20DetectorMessenger: public G4UImessenger
{
public:
  Tst20DetectorMessenger(Tst20DetectorConstruction* );
  ~Tst20DetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  Tst20DetectorConstruction* Tst20Detector;
  
  G4UIdirectory*             Tst20detDir;
  
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
  G4UIcmdWithADoubleAndUnit* TileDistanceCmd;
  
  // Anticoincidence

  G4UIcmdWithADoubleAndUnit* ACDThickCmd;
  G4UIcmdWithAString*        ACDMaterialCmd;

  // Total
  
  G4UIcmdWithADoubleAndUnit* MagFieldCmd;
  G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif










