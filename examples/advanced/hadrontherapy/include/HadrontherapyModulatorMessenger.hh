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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#ifndef HadrontherapyModulatorMessenger_h
#define HadrontherapyModulatorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class HadrontherapyModulator;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;

class HadrontherapyModulatorMessenger: public G4UImessenger
{
  public:
   HadrontherapyModulatorMessenger(HadrontherapyModulator*);
  ~HadrontherapyModulatorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:

  // Pointer to the modulator wheel
  HadrontherapyModulator* Modulator;

  G4UIdirectory* modulatorDir; // Control of the modulator 
  
  G4UIcmdWithADoubleAndUnit* modulatorAngleCmd;
  // UI command to rotate the modulator wheel

  G4UIcmdWithAString*   modulatorMatCmd;
  // UI command to set the material of the modulator wheel 
  
  G4UIcmdWithAString*   modulatorExternalFile;
  // UI command to set an external file for inputing modulator properties
  
  G4UIcmdWith3VectorAndUnit* modulatorPositionCmd;
  // UI command to change the position of the modulator  in the the beam line 


  G4UIcmdWithADoubleAndUnit* modulatorOuterRadiusCmd;
  // UI command to set the outer radius of the modulator wheel

  G4UIcmdWithADoubleAndUnit* modulatorInnerRadiusCmd;
  // UI command to set the inner radius of the modulator wheel
};
#endif

