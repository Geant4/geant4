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
// $Id: G4EnergyLossMessenger.hh,v 1.17 2007/03/17 19:24:39 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-03 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4EnergyLossMessenger
//
// Author:        Michel Maire
//
// Creation date: 22-06-2000
//
// Modifications:
// 10-01-06 SetStepLimits -> SetStepFunction (V.Ivantchenko)
// 10-05-06 Add command MscStepLimit (V.Ivantchenko) 
// 10-11-06 Add command binsDEDX (V.Ivantchenko) 
// 12-02-07 Add SetSkin, SetLinearLossLimit (V.Ivanchenko)
// 16-03-07 modify /process/eLoss/minsubsec command (V.Ivanchenko)
//
// -------------------------------------------------------------------
//

//
// Class Description:
//  This is a messenger class to interface to exchange information
//  between G4VEnergyLoss and UI.
//
//  /process/eLoss/   directory
//
//   Commands :
//
//    rndmStep *     Randomize the proposed step by eLoss (false/true)
//    fluct *        Switch true/false the energy loss fluctuations
//    subsec *       Switch true/false the subcutoff generation
//    minsubsec *    Set the min. cut for subcutoff delta in range
//    StepFunction * Set the energy loss step limitation parameters
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4EnergyLossMessenger_h
#define G4EnergyLossMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EnergyLossMessenger: public G4UImessenger
{
public:   // with description
  
  G4EnergyLossMessenger();
  virtual ~G4EnergyLossMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:

  G4UIdirectory*             eLossDirectory;
  G4UIcmdWithABool*          RndmStepCmd;
  G4UIcmdWithABool*          EnlossFlucCmd;
  G4UIcmdWithABool*          SubSecCmd;
  G4UIcmdWithADouble*        MinSubSecCmd;
  G4UIcommand*               StepFuncCmd;
  G4UIcommand*               mscCmd;
  G4UIcmdWithADoubleAndUnit* MinEnCmd;
  G4UIcmdWithADoubleAndUnit* MaxEnCmd;
  G4UIcmdWithABool*          IntegCmd;
  G4UIcmdWithABool*          rangeCmd;
  G4UIcmdWithABool*          lpmCmd;
  G4UIcmdWithABool*          latCmd;
  G4UIcmdWithAnInteger*      verCmd;
  G4UIcmdWithAnInteger*      dedxCmd;
  G4UIcmdWithAnInteger*      lbCmd;
  G4UIcmdWithADouble*        lllCmd;
  G4UIcmdWithADouble*        skinCmd;
};

#endif

