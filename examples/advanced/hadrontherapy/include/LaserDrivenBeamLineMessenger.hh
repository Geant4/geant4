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

#ifndef LaserDrivenBeamLineMessenger_h
#define LaserDrivenBeamLineMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class LaserDrivenBeamLine;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class LaserDrivenBeamLineMessenger: public G4UImessenger
{
public:
  LaserDrivenBeamLineMessenger(LaserDrivenBeamLine*);
	~LaserDrivenBeamLineMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:
	
	// Pointer to the detector component
	LaserDrivenBeamLine *laserDrivenMessengerPointer;
	
  G4UIdirectory *laserDrivenDir;
  G4UIdirectory *energySelectorDir;
  G4UIdirectory *FcollimatorDir;
  G4UIdirectory *ScollimatorDir;
  G4UIdirectory *slitDir;
  G4UIdirectory *quadrupoleDir;
  G4UIdirectory *relativePosDir;

  G4UIcmdWithoutParameter *DisableESSCmd;

  G4UIcmdWithADoubleAndUnit *FcollimatorRadiusCmd;
  G4UIcmdWithADoubleAndUnit *FcollimatorThicknessCmd;
  G4UIcmdWithADoubleAndUnit *FcollimatorZpositionCmd;
  G4UIcmdWithADoubleAndUnit *ScollimatorRadiusCmd;
  G4UIcmdWithADoubleAndUnit *ScollimatorThicknessCmd;
  G4UIcmdWithADoubleAndUnit *ScollimatorZpositionCmd;

  G4UIcmdWithADoubleAndUnit *SlitThicknessCmd;
  G4UIcmdWithADoubleAndUnit *holeSlitDimensionYCmd;
  G4UIcmdWithADoubleAndUnit *holeSlitDimensionZCmd;
  G4UIcmdWithADoubleAndUnit *slitHolePositionZCmd;

  G4UIcmdWithoutParameter *DisableQuadsCmd;

};
#endif
