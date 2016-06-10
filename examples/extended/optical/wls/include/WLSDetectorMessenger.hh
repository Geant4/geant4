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
// $Id: WLSDetectorMessenger.hh 84595 2014-10-17 07:33:27Z gcosmo $
//
/// \file optical/wls/include/WLSDetectorMessenger.hh
/// \brief Definition of the WLSDetectorMessenger class
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
    virtual ~WLSDetectorMessenger();
 
    virtual void SetNewValue(G4UIcommand*, G4String);

  private:

    WLSDetectorConstruction*   fDetector;
 
    G4UIdirectory*             fDetDir;

    G4UIcmdWithAString*        fSetPhotonDetGeometryCmd;
    G4UIcmdWithAnInteger*      fSetNumOfCladLayersCmd;
    G4UIcmdWithADoubleAndUnit* fSetWLSLengthCmd;
    G4UIcmdWithADoubleAndUnit* fSetWLSRadiusCmd;
    G4UIcmdWithADoubleAndUnit* fSetClad1RadiusCmd;
    G4UIcmdWithADoubleAndUnit* fSetClad2RadiusCmd;
    G4UIcmdWithADoubleAndUnit* fSetPhotonDetHalfLengthCmd;
    G4UIcmdWithADoubleAndUnit* fSetGapCmd;
    G4UIcmdWithADoubleAndUnit* fSetPhotonDetAlignmentCmd;
    G4UIcmdWithADouble*        fSetXYRatioCmd;
    G4UIcmdWithADouble*        fSetSurfaceRoughnessCmd;
    G4UIcmdWithADouble*        fSetMirrorPolishCmd;
    G4UIcmdWithADouble*        fSetMirrorReflectivityCmd;
    G4UIcmdWithADouble*        fSetPhotonDetPolishCmd;
    G4UIcmdWithADouble*        fSetPhotonDetReflectivityCmd;
    G4UIcmdWithABool*          fSetMirrorCmd;
    G4UIcmdWithADoubleAndUnit* fSetBarLengthCmd;
    G4UIcmdWithADoubleAndUnit* fSetBarBaseCmd;
    G4UIcmdWithADoubleAndUnit* fSetHoleRadiusCmd;
    G4UIcmdWithADoubleAndUnit* fSetCoatingThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fSetCoatingRadiusCmd;

};

#endif
