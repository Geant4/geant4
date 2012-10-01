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
/// \file optical/wls/include/WLSDetectorMessenger.hh
/// \brief Definition of the WLSDetectorMessenger class
//
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSDetectorMessenger_h
#define WLSDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include "WLSDetectorConstruction.hh"

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithoutParameter;

class WLSDetectorMessenger : public G4UImessenger
{
  public:

    WLSDetectorMessenger(WLSDetectorConstruction* );
    ~WLSDetectorMessenger();
 
    void SetNewValue(G4UIcommand*, G4String);

  private:

    WLSDetectorConstruction*   Detector;
 
    G4UIdirectory*          detDir;

    G4UIcmdWithoutParameter*   UpdateCmd;

    G4UIcmdWithAString*        SetPhotonDetGeometryCmd;
    G4UIcmdWithAnInteger*      SetNumOfCladLayersCmd;
    G4UIcmdWithADoubleAndUnit* SetWLSLengthCmd;
    G4UIcmdWithADoubleAndUnit* SetWLSRadiusCmd;
    G4UIcmdWithADoubleAndUnit* SetClad1RadiusCmd;
    G4UIcmdWithADoubleAndUnit* SetClad2RadiusCmd;
    G4UIcmdWithADoubleAndUnit* SetPhotonDetHalfLengthCmd;
    G4UIcmdWithADoubleAndUnit* SetGapCmd;
    G4UIcmdWithADoubleAndUnit* SetPhotonDetAlignmentCmd;
    G4UIcmdWithADouble*        SetXYRatioCmd;
    G4UIcmdWithADouble*        SetSurfaceRoughnessCmd;
    G4UIcmdWithADouble*        SetMirrorPolishCmd;
    G4UIcmdWithADouble*        SetMirrorReflectivityCmd;
    G4UIcmdWithADouble*        SetPhotonDetPolishCmd;
    G4UIcmdWithADouble*        SetPhotonDetReflectivityCmd;
    G4UIcmdWithABool*          SetMirrorCmd;
    G4UIcmdWithADoubleAndUnit* SetBarLengthCmd;
    G4UIcmdWithADoubleAndUnit* SetBarBaseCmd;
    G4UIcmdWithADoubleAndUnit* SetHoleRadiusCmd;
    G4UIcmdWithADoubleAndUnit* SetCoatingThicknessCmd;
    G4UIcmdWithADoubleAndUnit* SetCoatingRadiusCmd;

};

#endif
