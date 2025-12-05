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

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

class PrimaryGeneratorAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class PrimaryGeneratorMessenger: public G4UImessenger
{
public:
  
  PrimaryGeneratorMessenger(PrimaryGeneratorAction* gen);
  virtual ~PrimaryGeneratorMessenger() override;
    
  void SetNewValue(G4UIcommand* command, G4String newValue) override;

private:
  
  PrimaryGeneratorAction*  fGen;

  G4UIdirectory*             fBeamDir;
  G4UIcmdWithADoubleAndUnit* fKinECmd;
  G4UIcmdWithADoubleAndUnit* fDECmd;
  G4UIcmdWithADoubleAndUnit* fX0Cmd;
  G4UIcmdWithADoubleAndUnit* fY0Cmd;
  G4UIcmdWithADoubleAndUnit* fZ0Cmd;
  G4UIcmdWithADoubleAndUnit* fDXCmd;
  G4UIcmdWithADoubleAndUnit* fDYCmd;
  G4UIcmdWithADoubleAndUnit* fDZCmd;
  G4UIcmdWithAnInteger*      fVerboseCmd;

};

#endif

