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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it
//
// Adapted from Hadrontherapy example
// by Jacopo Magini: j.magini@surrey.ac.uk

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4String.hh"

#include "AnalysisManager.hh"

class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class DetectorMessenger: public G4UImessenger
{
public:
  DetectorMessenger(AnalysisManager* analysis);
  ~DetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
  
  // functions to read parameters
  G4String GetDetectorType() { return detectorType; }
  G4double GetDetectorSizeWidth() { return detectorWidth; }
  G4double GetDetectorSizeThickness() { return detectorThickness; }
  G4double GetSecondStageSizeWidth() { return secondStageWidth; }
  G4double GetSecondStageSizeThickness() { return secondStageThickness; }
  
  G4bool IsPhantomEnabled() { return phantomEnabled; }
  G4double GetDetectorPositionDepth() { return detectorDepth; }
  
  G4bool AreTherePendingChanges() { return geometryHasChanged; }
    
private:
  G4UIdirectory* changeTheGeometryDir;      // UI directory for the geometry control
  // subdirectories
  G4UIdirectory* changeDetectorDimensionDir;
  G4UIdirectory* changeDetectorPositionDir;
  G4UIdirectory* changeSecondStageDimensionDir;
  
  // commands
  G4UIcmdWithAString* changeTheDetectorCmd;
  G4UIcmdWithADoubleAndUnit* changeDetectorSizeWidthCmd;
  G4UIcmdWithADoubleAndUnit* changeDetectorSizeThicknessCmd;
  G4UIcmdWithADoubleAndUnit* changeSecondStageSizeWidthCmd;
  G4UIcmdWithADoubleAndUnit* changeSecondStageSizeThicknessCmd;
  
  G4UIcmdWithABool* enableWaterPhantomCmd;
  
  G4UIcmdWithADoubleAndUnit* changeDetectorPositionDepthCmd;
  
  G4UIcmdWithoutParameter *applyChangesToGeometryCmd;	// applies changes to detector
  
  // parameters to store
  G4String detectorType;
  G4double detectorWidth;
  G4double detectorThickness;
  G4double secondStageWidth;
  G4double secondStageThickness;
  G4bool phantomEnabled;
  G4double detectorDepth;
  
  G4bool geometryHasChanged;
  
  AnalysisManager* analysis;
};
#endif

