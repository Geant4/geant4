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
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelPrimaryGeneratorMessenger  ------
//           by G.Santin, F.Longo & R.Giannitrapani (13 nov 2000) 
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelPrimaryGeneratorMessenger_h
#define GammaRayTelPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class GammaRayTelPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  GammaRayTelPrimaryGeneratorMessenger(GammaRayTelPrimaryGeneratorAction*);
  ~GammaRayTelPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  GammaRayTelPrimaryGeneratorAction* GammaRayTelAction; 
  G4UIcmdWithAString*          RndmCmd;
  G4UIcmdWithAnInteger*        SourceTypeCmd;
  G4UIcmdWithADoubleAndUnit*   VertexRadiusCmd;
  G4UIcmdWithAnInteger*        SpectrumTypeCmd;
  G4UIcmdWithABool*            SourceGenCmd;

};

#endif


