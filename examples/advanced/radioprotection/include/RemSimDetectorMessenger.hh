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
//    *****************************************
//    *                                       *
//    *      RemSimDetectrorMessenger.hh      *
//    *                                       *
//    *****************************************
//
// $Id: RemSimDetectorMessenger.hh,v 1.8 2004/05/22 12:57:04 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

