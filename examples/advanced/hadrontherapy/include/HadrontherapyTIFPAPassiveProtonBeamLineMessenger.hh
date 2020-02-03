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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#ifndef TrentoPassiveProtonBeamLineMessenger_h
#define TrentoPassiveProtonBeamLineMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class TrentoPassiveProtonBeamLine;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

class TrentoPassiveProtonBeamLineMessenger: public G4UImessenger
{
  public:
  TrentoPassiveProtonBeamLineMessenger(TrentoPassiveProtonBeamLine*);
  ~TrentoPassiveProtonBeamLineMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:

  // Pointer to the detector component
  TrentoPassiveProtonBeamLine* TrentoPassiveProton;

  G4UIdirectory* changeTheBeamLineDir;
  G4UIcmdWithAString*        changeTheBeamLineNameCmd; // Control the name of the beam line 

  G4UIdirectory* modulatorDir; // Control of the modulator 
  G4UIdirectory* beamLineDir;  // Control of the beam line
 
  //G4UIdirectory* ScatteringFoilDir;
  // Control of the first scattering foil component of the beam line
  
  G4UIcmdWithADoubleAndUnit* ScatteringFoilXSizeCmd;
  // UI command to set half of the X size of the rangeShifter component of
  // the beam line
    
  G4UIcmdWithAString* scatteringFoilMatCmd;
  // UI command to set the material of the rangeShifter component of
  // the beam line
    
  G4UIcmdWithADoubleAndUnit* preCollimatorXSizeCmd;
    // UI command to set half of the X size of the pre collimator component of
    // the beam line
    
  G4UIcmdWithADoubleAndUnit* preCollimatorXPositionCmd;
  // UI command to set the position X of the pre collimator in
  // the beam line
    
  G4UIcmdWithADoubleAndUnit* AirTubeYSizeCmd;
  // UI command to set the position X of the pre collimator in
  // the beam line

  G4UIcmdWithADoubleAndUnit* AirTubeZSizeCmd;
  // UI command to set outer radius of collimator

};
#endif

