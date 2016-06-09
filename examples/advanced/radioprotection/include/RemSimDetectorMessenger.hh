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
//    *****************************************
//    *                                       *
//    *      RemSimDetectrorMessenger.hh      *
//    *                                       *
//    *****************************************
//
// $Id$
//
// 
#ifndef RemSimDetectorMessenger_h
#define RemSimDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class RemSimDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class RemSimDetectorMessenger: public G4UImessenger
{
public:
  RemSimDetectorMessenger(RemSimDetectorConstruction* );
  ~RemSimDetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  RemSimDetectorConstruction*  detector;//pointer to detector
  G4UIdirectory*               vehicleDir; // control the geometry configuration
  G4UIcmdWithAString*          vehicleCmd; //change set-up configuration 
  G4UIcmdWithAString*          shieldingCmd;//add the shielding layer
  G4UIcmdWithAString*          SPECmd; //add the SPE shelter
  G4UIcmdWithAString*          roofCmd; // add the roof to the moon habitat
  G4UIdirectory*               shieldingDir;//control the shielding parameters 
  G4UIcmdWithADoubleAndUnit*   thicknessCmd; //change the thickness of the
                                             // shielding

 G4UIdirectory*  roofDir; // control the roof height
 G4UIcmdWithADoubleAndUnit*   thicknessRoofCmd;//change the thickness of the
                                             // roof
};
#endif

