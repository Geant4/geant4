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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLXXInterfaceMessenger.cc
 * \brief Messenger class for the Geant4 INCL++ interface.
 *
 * \date 26th April 2012
 * \author Davide Mancusi
 */

#include "G4INCLXXInterfaceMessenger.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

const G4String G4INCLXXInterfaceMessenger::theUIDirectory = "/process/had/inclxx/";

G4INCLXXInterfaceMessenger::G4INCLXXInterfaceMessenger(G4INCLXXInterfaceStore *anInterfaceStore) :
  theINCLXXInterfaceStore(anInterfaceStore)
{
  // Create a directory for the INCL++ commands
  theINCLXXDirectory = new G4UIdirectory(theUIDirectory);
  theINCLXXDirectory->SetGuidance("Parameters for the INCL++ model");

  // This command controls whether nucleus-nucleus reactions should accurately
  // describe the projectile or the target nucleus (default: projectile)
  accurateNucleusCmd = new G4UIcmdWithAString((theUIDirectory + "accurateNucleus").data(),this);
  accurateNucleusCmd->SetGuidance("Set which nucleus will be accurately described in nucleus-nucleus reactions.");
  accurateNucleusCmd->SetGuidance(" projectile: accurate description of projectile-related quantities");
  accurateNucleusCmd->SetGuidance(" target: accurate description of target-related quantities");
  accurateNucleusCmd->SetGuidance(" Default: projectile");
  accurateNucleusCmd->SetParameterName("AccurateNucleus",true);
  accurateNucleusCmd->SetDefaultValue("projectile");
  accurateNucleusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // This command selects the maximum mass number of clusters to be produced in
  // the INCL++ cascade
  maxClusterMassCmd = new G4UIcmdWithAnInteger((theUIDirectory + "maxClusterMass").data(),this);
  maxClusterMassCmd->SetGuidance("Set the maximum cluster mass.");
  maxClusterMassCmd->SetGuidance(" The INCL++ cascade stage will produce clusters with mass up to the value of this parameter (included)");
  maxClusterMassCmd->SetGuidance(" Allowed range: [2,12]");
  maxClusterMassCmd->SetParameterName("MaxClusterMass",true);
  maxClusterMassCmd->SetDefaultValue(8);
  maxClusterMassCmd->SetRange("MaxClusterMass>=2 && MaxClusterMass<=12");
  maxClusterMassCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // This command sets the energy below which PreCoumpound will be used
  cascadeMinEnergyPerNucleonCmd = new G4UIcmdWithADoubleAndUnit((theUIDirectory + "cascadeMinEnergyPerNucleon").data(),this);
  cascadeMinEnergyPerNucleonCmd->SetGuidance("Set the minimum energy per nucleon at which cascade will be used.");
  cascadeMinEnergyPerNucleonCmd->SetGuidance(" INCL++ will rely on PreCompound for reactions induced by projectiles slower than the given energy (per nucleon, where applicable)");
  cascadeMinEnergyPerNucleonCmd->SetParameterName("cascadeMinEnergyPerNucleon",true);
  cascadeMinEnergyPerNucleonCmd->SetDefaultValue(1.*MeV);
  cascadeMinEnergyPerNucleonCmd->SetRange("cascadeMinEnergyPerNucleon>=0");
  cascadeMinEnergyPerNucleonCmd->SetUnitCategory("Energy");
  cascadeMinEnergyPerNucleonCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // This command allows the user to set several INCL++ internals in one shot
  inclPhysicsCmd = new G4UIcmdWithAString((theUIDirectory + "setPhysics").data(),this);
  inclPhysicsCmd->SetGuidance("Set a global configuration for INCL++.");
  inclPhysicsCmd->SetGuidance(" default: default configuration, most recent options");
  inclPhysicsCmd->SetGuidance(" incl42: try to mimick the behaviour from INCL4.2");
  inclPhysicsCmd->SetParameterName("type", false);
  inclPhysicsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // This command allows the user to change the de-excitation model to be used
  // with INCL++
  useAblaCmd = new G4UIcommand((theUIDirectory + "useAbla").data(),this);
  useAblaCmd->SetGuidance("Use ABLA++ as de-excitation model after INCL++.");
  useAblaCmd->AvailableForStates(G4State_Idle);
}

G4INCLXXInterfaceMessenger::~G4INCLXXInterfaceMessenger() {
  delete theINCLXXDirectory;
  delete accurateNucleusCmd;
  delete maxClusterMassCmd;
  delete cascadeMinEnergyPerNucleonCmd;
  delete inclPhysicsCmd;
  delete useAblaCmd;
}

void G4INCLXXInterfaceMessenger::SetNewValue(G4UIcommand *command, G4String newValues) {
  if(command==accurateNucleusCmd) {
    newValues.toLower();
    if(newValues=="projectile") {
      theINCLXXInterfaceStore->SetAccurateProjectile(true);
    } else if(newValues=="target") {
      theINCLXXInterfaceStore->SetAccurateProjectile(false);
    }
  } else if(command==maxClusterMassCmd) {
    const G4int parameter = maxClusterMassCmd->GetNewIntValue(newValues);
    theINCLXXInterfaceStore->SetMaxClusterMass(parameter);
  } else if(command==cascadeMinEnergyPerNucleonCmd) {
    const G4double parameter = cascadeMinEnergyPerNucleonCmd->GetNewDoubleValue(newValues);
    theINCLXXInterfaceStore->SetCascadeMinEnergyPerNucleon(parameter);
  } else if(command==inclPhysicsCmd) {
    theINCLXXInterfaceStore->SetINCLPhysics(newValues);
  } else if(command==useAblaCmd) {
    theINCLXXInterfaceStore->UseAblaDeExcitation();
  }
}
