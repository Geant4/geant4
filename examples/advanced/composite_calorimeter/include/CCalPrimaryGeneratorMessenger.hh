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
///////////////////////////////////////////////////////////////////////////////
// File: CCalPrimaryGeneratorMessenger.hh
// Description: Adds a new command to (un)select random shooting.
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalPrimaryGeneratorMessenger_h
#define CCalPrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CCalPrimaryGeneratorAction;

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

class CCalPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  CCalPrimaryGeneratorMessenger(CCalPrimaryGeneratorAction* myGun);
  ~CCalPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand * command,G4String newValues);
  
private:
  CCalPrimaryGeneratorAction* myAction;
  
  G4UIcmdWithAnInteger*      verboseCmd;
  G4UIcmdWithAString*        rndmCmd;
  G4UIcmdWithAString*        scanCmd;
  G4UIcmdWithADoubleAndUnit* minEnergyCmd;
  G4UIcmdWithADoubleAndUnit* maxEnergyCmd;
  G4UIcmdWithADoubleAndUnit* minPhiCmd;
  G4UIcmdWithADoubleAndUnit* maxPhiCmd;
  G4UIcmdWithADouble*        minEtaCmd; 
  G4UIcmdWithADouble*        maxEtaCmd; 
  G4UIcmdWithAnInteger*      stepsPhiCmd;
  G4UIcmdWithAnInteger*      stepsEtaCmd;
  G4UIcmdWithAnInteger*      runNoCmd;
};

#endif
