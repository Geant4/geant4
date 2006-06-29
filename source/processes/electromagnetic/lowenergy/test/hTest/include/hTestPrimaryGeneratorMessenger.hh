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
#ifndef hTestPrimaryGeneratorMessenger_h
#define hTestPrimaryGeneratorMessenger_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestPrimaryGeneratorAction
//  
// Description: Definition of physics list parameters via UI interface
//
// Author:      V.Ivanchenko 26/09/00
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestPrimaryGeneratorAction;

class hTestPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
  
    hTestPrimaryGeneratorMessenger(hTestPrimaryGeneratorAction* gen);
   ~hTestPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand* command, G4String newValue);

  private:
  
    hTestPrimaryGeneratorAction*  theGen;

    G4UIcmdWithADoubleAndUnit* beamXCmd;
    G4UIcmdWithADoubleAndUnit* beamYCmd;
    G4UIcmdWithADoubleAndUnit* beamZCmd;
    G4UIcmdWithADoubleAndUnit* beamECmd;
    G4UIcmdWithADoubleAndUnit* sigmaXCmd;
    G4UIcmdWithADoubleAndUnit* sigmaYCmd;
    G4UIcmdWithADoubleAndUnit* sigmaZCmd;
    G4UIcmdWithADoubleAndUnit* sigmaECmd;
    G4UIcmdWithADoubleAndUnit* maxThetaCmd;
    G4UIcmdWithADouble* beamBetaCmd;
    G4UIcmdWithADouble* sigmaBetaCmd;
    G4UIcmdWithAString* partCmd;
    G4UIcmdWithAString* randCmd;

};

#endif

