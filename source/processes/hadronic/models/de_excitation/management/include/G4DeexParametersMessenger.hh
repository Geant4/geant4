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
// $Id: G4DeexParametersMessenger.hh 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4DeexParametersMessenger
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 17-10-2017
//
// Modifications:
//
// -------------------------------------------------------------------
//

// Class Description:
//  This is a messenger class to interface to exchange information
//  between deexcitation module and UI.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4DeexParametersMessenger_h
#define G4DeexParametersMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4DeexPrecoParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DeexParametersMessenger: public G4UImessenger
{
public:   // with description
  
  explicit G4DeexParametersMessenger(G4DeexPrecoParameters*);
  virtual ~G4DeexParametersMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String) override;

private:

  G4DeexPrecoParameters*     theParameters;

  G4UIdirectory*             fDirectory;

  G4UIcmdWithABool*          readCmd;
  G4UIcmdWithABool*          icCmd;
  G4UIcmdWithABool*          corgCmd;

  G4UIcmdWithAnInteger*      maxjCmd;

};

#endif

