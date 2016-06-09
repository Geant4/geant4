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
/// \file optical/wls/src/WLSPrimaryGeneratorMessenger.cc
/// \brief Implementation of the WLSPrimaryGeneratorMessenger class
//
//
//

#include "G4UIdirectory.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"

#include "WLSPrimaryGeneratorAction.hh"
#include "WLSPrimaryGeneratorMessenger.hh"

WLSPrimaryGeneratorMessenger::
                   WLSPrimaryGeneratorMessenger(WLSPrimaryGeneratorAction* gun)
  : action(gun)
{
  gunDir = new G4UIdirectory("/WLS/gun/");
  gunDir->SetGuidance("WLSPrimaryGenerator control");

  SetPolarizationCmd =
                 new G4UIcmdWithADoubleAndUnit("/WLS/gun/optPhotonPolar",this);
  SetPolarizationCmd->SetGuidance("Set linear polarization");
  SetPolarizationCmd->SetGuidance("  angle w.r.t. (k,n) plane");
  SetPolarizationCmd->SetParameterName("angle",true);
  SetPolarizationCmd->SetUnitCategory("Angle");
  SetPolarizationCmd->SetDefaultValue(0.);
  SetPolarizationCmd->AvailableForStates(G4State_Idle);
 
  SetDecayTimeConstantCmd =
           new G4UIcmdWithADoubleAndUnit("/WLS/gun/setDecayTimeConstant",this);
  SetDecayTimeConstantCmd->SetGuidance("Set the decay time constant");
  SetDecayTimeConstantCmd->SetGuidance("for the starting time of each photon");
  SetDecayTimeConstantCmd->SetParameterName("time_const",false);
  SetDecayTimeConstantCmd->SetUnitCategory("Time");
  SetDecayTimeConstantCmd->SetRange("time_const>=0");
  SetDecayTimeConstantCmd->AvailableForStates(G4State_Idle);
}

WLSPrimaryGeneratorMessenger::~WLSPrimaryGeneratorMessenger()
{
  delete gunDir;
  delete SetPolarizationCmd;
  delete SetDecayTimeConstantCmd;
}

void WLSPrimaryGeneratorMessenger::
                           SetNewValue(G4UIcommand * command,G4String val)
{
  if ( command == SetPolarizationCmd )
     action->
          SetOptPhotonPolar(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
  else if ( command == SetDecayTimeConstantCmd )
     action->
       SetDecayTimeConstant(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
}
