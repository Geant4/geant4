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
// File name:     G4EmLowEParametersMessenger
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

#ifndef G4EmLowEParametersMessenger_h
#define G4EmLowEParametersMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;
class G4EmLowEParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EmLowEParametersMessenger: public G4UImessenger
{
public:   // with description
  
  explicit G4EmLowEParametersMessenger(G4EmLowEParameters*);
  ~G4EmLowEParametersMessenger() override;

  void SetNewValue(G4UIcommand*, G4String) override;

  G4EmLowEParametersMessenger & operator=
  (const G4EmLowEParametersMessenger &right) = delete;
  G4EmLowEParametersMessenger(const G4EmLowEParametersMessenger&) = delete;

private:

  G4EmLowEParameters*        theParameters;

  G4UIcmdWithABool*          deCmd;
  G4UIcmdWithABool*          dirFluoCmd;
  G4UIcmdWithABool*          dirFluoCmd1;
  G4UIcmdWithABool*          auCmd;
  G4UIcmdWithABool*          auCascadeCmd;
  G4UIcmdWithABool*          pixeCmd;
  G4UIcmdWithABool*          dcutCmd;
  G4UIcmdWithABool*          dnafCmd;
  G4UIcmdWithABool*          dnasCmd;
  G4UIcmdWithABool*          dnamscCmd;

  G4UIcmdWithAString*        pixeXsCmd;
  G4UIcmdWithAString*        pixeeXsCmd;
  G4UIcmdWithAString*        livCmd;
  G4UIcmdWithAString*        dnaSolCmd;
  G4UIcmdWithAString*        direFluoCmd;

  G4UIcmdWithAString*        meCmd;
  G4UIcommand*               dnaCmd;
  G4UIcommand*               deexCmd;
};

#endif

