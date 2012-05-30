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
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.1_rc11
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLXXInterfaceMessenger.cc
 * \brief Messenger class for the Geant4 INCL++ interface.
 *
 * Created on: 26th April 2012
 *     Author: Davide Mancusi
 */

#include "G4INCLXXInterfaceMessenger.hh"
#include <sstream>

G4INCLXXInterfaceMessenger::G4INCLXXInterfaceMessenger(G4INCLXXInterfaceConfig *anInterfaceConfig) :
  theINCLXXInterfaceConfig(anInterfaceConfig)
{
  // Create a directory for the INCL++ commands
  theINCLXXDirectory = new G4UIdirectory("/inclxx/");
  theINCLXXDirectory->SetGuidance("Controls for the INCL++ interface.");

  // This command controls whether nucleus-nucleus reactions should be
  // simulated in inverse kinematics (default: true)
  inverseKinematicsCmd = new G4UIcmdWithABool ("/inclxx/useInverseKinematics",this);
  inverseKinematicsCmd->SetGuidance("Set how the model simulates nucleus-nuclues reactions.");
  inverseKinematicsCmd->SetGuidance(" True: simulate nucleus-nucleus reactions in inverse kinematics (i.e. treat the projectile as the target, and vice-versa)");
  inverseKinematicsCmd->SetGuidance(" False: don't perform any switch");
  inverseKinematicsCmd->SetGuidance(" Default: true.");
  inverseKinematicsCmd->SetParameterName("UseInverseKinematics",true);
  inverseKinematicsCmd->SetDefaultValue(true);

  // This command selects the de-excitation model to be used
/*  deexcitationModelCmd = new G4UIcmdWithAString ("/inclxx/deexcitationModel",this);
  deexcitationModelCmd->SetGuidance("Select the de-excitation model to be used for the afterburning stage.");
  deexcitationModelCmd->SetGuidance(" g4 (G4ExcitationHandler; default)");
//  deexcitationModelCmd->SetGuidance(" abla (ABLA model)");
  deexcitationModelCmd->SetParameterName("DeexcitationModel",true);
  deexcitationModelCmd->SetDefaultValue("g4");*/

  // This command selects the maximum mass of clusters to be produced in the
  // INCL++ cascade
  maxClusterMassCmd = new G4UIcmdWithAnInteger ("/inclxx/maxClusterMass",this);
  maxClusterMassCmd->SetGuidance("Set the maximum cluster mass.");
  maxClusterMassCmd->SetGuidance(" The INCL++ cascade stage will produce clusters with mass up to the value of this parameter (included)");
  maxClusterMassCmd->SetGuidance(" Allowed range: [2,12]");
  maxClusterMassCmd->SetParameterName("MaxClusterMass",true);
  maxClusterMassCmd->SetDefaultValue(8);
  maxClusterMassCmd->SetRange("MaxClusterMass>=2 && MaxClusterMass<=12");

  // This command selects the maximum mass of clusters to be produced in the
  // INCL++ cascade
  maxProjMassCmd = new G4UIcmdWithAnInteger ("/inclxx/maxProjMass",this);
  maxProjMassCmd->SetGuidance("Set the maximum projectile mass.");
  maxProjMassCmd->SetGuidance(" Use INCL++ for nucleus-induced reactions up to this projectile mass (included).");
  maxProjMassCmd->SetParameterName("MaxProjMass",true);
  maxProjMassCmd->SetDefaultValue(18);
  maxProjMassCmd->SetRange("MaxProjMass>=2 && MaxProjMass<=18");

}

G4INCLXXInterfaceMessenger::~G4INCLXXInterfaceMessenger() {
  delete theINCLXXDirectory;
  delete inverseKinematicsCmd;
  delete maxClusterMassCmd;
  delete maxProjMassCmd;
}

void G4INCLXXInterfaceMessenger::SetNewValue(G4UIcommand *command, G4String newValues) {
  if(command==inverseKinematicsCmd) {
    const G4bool parameter = inverseKinematicsCmd->GetNewBoolValue(newValues);
    if(parameter)
      G4cout << "INCL++ interface: will use inverse kinematics in nucleus-nucleus reactions." << G4endl;
    else
      G4cout << "INCL++ interface: will not use inverse kinematics in nucleus-nucleus reactions." << G4endl;
    theINCLXXInterfaceConfig->SetUseInverseKinematics(parameter);
/*  } else if(command==deexcitationModelCmd) {
    if(newValues=="g4") {
      G4cout << "INCL++ interface: setting de-excitation model to G4ExcitationHandler." << G4endl;
      theINCLXXInterface->SetDeexcitationModel(G4INCLXXInterface::G4ExcitationHandlerModel);
    } else if(newValues=="abla") {
      G4cout << "INCL++ interface: setting de-excitation model to ABLA." << G4endl;
      theINCLXXInterface->SetDeexcitationModel(G4INCLXXInterface::ABLAModel);
    } else {
      G4cout << "INCL++ interface: unrecognised de-excitation model in /inclxx/deexcitationModel command. Ignoring it." << G4endl;
    }*/
  } else if(command==maxClusterMassCmd) {
    const G4int parameter = maxClusterMassCmd->GetNewIntValue(newValues);
    G4cout << "INCL++ interface: setting maximum cluster mass to " << parameter << "." << G4endl;
    theINCLXXInterfaceConfig->SetMaxClusterMass(parameter);
  } else if(command==maxProjMassCmd) {
    const G4int parameter = maxProjMassCmd->GetNewIntValue(newValues);
    G4cout << "INCL++ interface: setting maximum projectile mass to " << parameter << "." << G4endl;
    theINCLXXInterfaceConfig->SetMaxProjMass(parameter);
  }

  theINCLXXInterfaceConfig->DeleteModels();
}
