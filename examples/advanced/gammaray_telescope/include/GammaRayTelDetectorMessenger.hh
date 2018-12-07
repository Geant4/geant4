//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDetectorMessenger  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************
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










