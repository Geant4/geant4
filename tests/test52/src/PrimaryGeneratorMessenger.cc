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

#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGenerator.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGenerator* primGen) :
    primaryGenerator(primGen) {

  sourceDirectory = new G4UIdirectory("/source/");
  sourceDirectory -> SetGuidance("Particle source commands");

  primEnergyCmd = 
             new G4UIcmdWithADoubleAndUnit("/source/primEnergy", this);
  primEnergyCmd -> SetGuidance("Specification of primary electron energy");
  primEnergyCmd -> SetParameterName("primEnergy", true);
  primEnergyCmd -> SetDefaultValue(1.0);
  primEnergyCmd -> SetDefaultUnit("MeV");
  primEnergyCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);  

  sigmaEnergyCmd = 
             new G4UIcmdWithADoubleAndUnit("/source/sigmaEnergy", this);
  sigmaEnergyCmd -> SetGuidance("Specification of sigma (energy distr.)");
  sigmaEnergyCmd -> SetParameterName("sigmaEnergy", true);
  sigmaEnergyCmd -> SetDefaultValue(0.0);
  sigmaEnergyCmd -> SetDefaultUnit("keV");
  sigmaEnergyCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);  

  sigmaSpatialCmd = 
              new G4UIcmdWithADoubleAndUnit("/source/sigmaSpatial", this);
  sigmaSpatialCmd -> SetGuidance("Specification of sigma (energy distr.)");
  sigmaSpatialCmd -> SetParameterName("sigmaSpatial", true);
  sigmaSpatialCmd -> SetDefaultValue(0.0);
  sigmaSpatialCmd -> SetDefaultUnit("mm");
  sigmaSpatialCmd -> AvailableForStates(G4State_PreInit, G4State_Idle); 

  incidAngleCmd = 
              new G4UIcmdWithADoubleAndUnit("/source/incidentAngle", this);
  incidAngleCmd -> SetGuidance("Specification of incident beam angle");
  incidAngleCmd -> SetParameterName("incidentAngle", true);
  incidAngleCmd -> SetDefaultValue(0.0);
  incidAngleCmd -> SetDefaultUnit("deg");
  incidAngleCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);  
}


PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {

  delete incidAngleCmd;
  delete sigmaSpatialCmd;   
  delete sigmaEnergyCmd;
  delete primEnergyCmd;
  delete sourceDirectory;
}


void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* cmd, G4String val) {

  if(cmd == primEnergyCmd) 
     primaryGenerator -> SetPrimaryKineticEnergy(
                            primEnergyCmd -> GetNewDoubleValue(val));
  if(cmd == sigmaEnergyCmd) 
     primaryGenerator -> SetSigmaKineticEnergy(
                            sigmaEnergyCmd -> GetNewDoubleValue(val));
  if(cmd == sigmaSpatialCmd) 
     primaryGenerator -> SetSigmaSpatialPlacement(
                            sigmaSpatialCmd -> GetNewDoubleValue(val));
  if(cmd == incidAngleCmd) 
     primaryGenerator -> SetIncidentAngle(
                            incidAngleCmd -> GetNewDoubleValue(val));
}
