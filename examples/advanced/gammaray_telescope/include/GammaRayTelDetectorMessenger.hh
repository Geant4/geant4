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
//
// $Id: GammaRayTelDetectorMessenger.hh,v 1.3 2001-07-11 09:56:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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










