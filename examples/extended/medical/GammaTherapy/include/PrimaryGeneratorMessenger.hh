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
// $Id: PrimaryGeneratorMessenger.hh 103795 2017-04-27 13:38:36Z gcosmo $
//
/// \file medical/GammaTherapy/include/PrimaryGeneratorMessenger.hh
/// \brief Definition of the PrimaryGeneratorMessenger class
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
  virtual ~PrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  
  PrimaryGeneratorMessenger & operator=(const PrimaryGeneratorMessenger &right);
  PrimaryGeneratorMessenger(const PrimaryGeneratorMessenger&);

  PrimaryGeneratorAction*  fGen;
  G4bool fVerbose;

  G4UIcmdWithADoubleAndUnit* fBeamXCmd;
  G4UIcmdWithADoubleAndUnit* fBeamYCmd;
  G4UIcmdWithADoubleAndUnit* fBeamZCmd;
  G4UIcmdWithADoubleAndUnit* fBeamECmd;
  G4UIcmdWithADoubleAndUnit* fSigmaXCmd;
  G4UIcmdWithADoubleAndUnit* fSigmaYCmd;
  G4UIcmdWithADoubleAndUnit* fSigmaZCmd;
  G4UIcmdWithADoubleAndUnit* fSigmaECmd;

  G4UIcmdWithAString* fRandCmd;
  G4UIcmdWithADoubleAndUnit* fMaxThetaCmd;
  G4UIcmdWithADoubleAndUnit* fThetaCmd;
};

#endif

