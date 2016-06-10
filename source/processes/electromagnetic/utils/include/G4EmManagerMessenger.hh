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
// $Id: G4EmManagerMessenger.hh 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4EmManagerMessenger
//
// Author:        Vladimir Ivanchenko created from G4EnergyLossMessenger
//
// Creation date: 22-05-2013
//
// Modifications:
//
// -------------------------------------------------------------------
//

// Class Description:
//  This is a messenger class to interface to exchange information
//  between G4VEm and UI.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4EmManagerMessenger_h
#define G4EmManagerMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4EmManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EmManagerMessenger: public G4UImessenger
{
public:   // with description
  
  G4EmManagerMessenger(G4EmManager*);
  virtual ~G4EmManagerMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

private:

  G4EmManager*               theManager;

  G4UIdirectory*             eLossDirectory;
  G4UIdirectory*             mscDirectory;
  G4UIdirectory*             emDirectory;
  G4UIcmdWithABool*          RndmStepCmd;
  G4UIcmdWithABool*          EnlossFlucCmd;
  G4UIcmdWithABool*          SubSecCmd;
  G4UIcmdWithADouble*        MinSubSecCmd;
  G4UIcommand*               StepFuncCmd;
  G4UIcmdWithADoubleAndUnit* MinEnCmd;
  G4UIcmdWithADoubleAndUnit* MaxEnCmd;
  G4UIcmdWithABool*          IntegCmd;
  G4UIcmdWithABool*          rangeCmd;
  G4UIcmdWithABool*          lpmCmd;
  G4UIcmdWithABool*          splCmd;
  G4UIcommand*               deexCmd;
  /*
  G4UIcmdWithABool*          aplCmd;
  G4UIcmdWithABool*          deCmd;
  G4UIcmdWithABool*          auCmd;
  G4UIcmdWithABool*          pixeCmd;
  G4UIcmdWithAString*        pixeXsCmd;
  G4UIcmdWithAString*        pixeeXsCmd;
  G4UIcmdWithAnInteger*      verCmd;
  G4UIcmdWithAnInteger*      ver1Cmd;
*/
  G4UIcmdWithAnInteger*      dedxCmd;
  G4UIcmdWithAnInteger*      lamCmd;
  G4UIcmdWithADouble*        lllCmd;
  /*
  G4UIcmdWithAString*        mscCmd;
  G4UIcmdWithADouble*        labCmd;
  G4UIcmdWithABool*          latCmd;
  G4UIcmdWithADouble*        skinCmd;
  G4UIcmdWithADouble*        frCmd;
  G4UIcmdWithADouble*        fgCmd;
  G4UIcmdWithADouble*        mscfCmd;
  G4UIcmdWithADoubleAndUnit* angCmd;
  G4UIcommand*               bfCmd;
  G4UIcommand*               fiCmd;
  G4UIcommand*               brCmd;
*/
};

#endif

