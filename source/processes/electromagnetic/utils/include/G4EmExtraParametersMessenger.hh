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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4EmExtraParametersMessenger
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 07-05-2019
//
// -------------------------------------------------------------------
//

// Class Description:
//  This is a messenger class to interface to exchange information
//  between G4VEm and UI.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4EmExtraParametersMessenger_h
#define G4EmExtraParametersMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;
class G4EmExtraParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EmExtraParametersMessenger: public G4UImessenger
{
public:   // with description
  
  explicit G4EmExtraParametersMessenger(G4EmExtraParameters*);
  ~G4EmExtraParametersMessenger() override;

  void SetNewValue(G4UIcommand*, G4String) override;

  G4EmExtraParametersMessenger & operator=
  (const G4EmExtraParametersMessenger &right) = delete;
  G4EmExtraParametersMessenger(const G4EmExtraParametersMessenger&) = delete;

private:

  G4EmExtraParameters*       theParameters;

  G4UIcmdWithABool*          dirSplitCmd;
  G4UIcmdWithABool*          qeCmd;

  G4UIcmdWithADoubleAndUnit* dirSplitRadiusCmd;

  G4UIcommand*               paiCmd;
  G4UIcommand*               mscoCmd;

  G4UIcmdWithAString*        SubSecCmd;
  G4UIcommand*               bfCmd;
  G4UIcommand*               fiCmd;
  G4UIcommand*               bsCmd;
  G4UIcommand*               StepFuncCmd;
  G4UIcommand*               StepFuncCmd1;
  G4UIcommand*               StepFuncCmd2;
  G4UIcommand*               StepFuncCmd3;

  G4UIcmdWith3VectorAndUnit* dirSplitTargetCmd;
};

#endif

