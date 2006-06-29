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

//---------------------------------------------------------------------------
//
// ClassName:   PrimaryGeneratorAction
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

class PrimaryGeneratorAction;

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
  
    PrimaryGeneratorMessenger(PrimaryGeneratorAction* gen);
   ~PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand* command, G4String newValue);

  private:
  
    PrimaryGeneratorAction*  theGen;

    G4UIcmdWithADoubleAndUnit* beamXCmd;
    G4UIcmdWithADoubleAndUnit* beamYCmd;
    G4UIcmdWithADoubleAndUnit* beamZCmd;
    G4UIcmdWithADoubleAndUnit* beamECmd;
    G4UIcmdWithADoubleAndUnit* sigmaXCmd;
    G4UIcmdWithADoubleAndUnit* sigmaYCmd;
    G4UIcmdWithADoubleAndUnit* sigmaZCmd;
    G4UIcmdWithADoubleAndUnit* sigmaECmd;
    G4UIcmdWithADoubleAndUnit* maxThetaCmd;
    G4UIcmdWithADoubleAndUnit* sThetaCmd;
    G4UIcmdWithADouble* beamBetaCmd;
    G4UIcmdWithADouble* sigmaBetaCmd;
    G4UIcmdWithAString* partCmd;
    G4UIcmdWithAString* randCmd;

};

#endif

