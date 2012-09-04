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
// INCL++ revision: v5.1.3
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

  // This command controls whether nucleus-nucleus reactions should accurately
  // describe the projectile or the target nucleus (default: projectile)
  accurateNucleusCmd = new G4UIcmdWithAString ("/inclxx/accurateNucleus",this);
  accurateNucleusCmd->SetGuidance("Set which nucleus will be accurately described in nucleus-nucleus reactions.");
  accurateNucleusCmd->SetGuidance(" projectile: accurate description of projectile-related quantities");
  accurateNucleusCmd->SetGuidance(" target: accurate description of target-related quantities");
  accurateNucleusCmd->SetGuidance(" Default: projectile");
  accurateNucleusCmd->SetParameterName("AccurateNucleus",true);
  accurateNucleusCmd->SetDefaultValue("projectile");

  // This command selects the maximum mass number of clusters to be produced in
  // the INCL++ cascade
  maxClusterMassCmd = new G4UIcmdWithAnInteger ("/inclxx/maxClusterMass",this);
  maxClusterMassCmd->SetGuidance("Set the maximum cluster mass.");
  maxClusterMassCmd->SetGuidance(" The INCL++ cascade stage will produce clusters with mass up to the value of this parameter (included)");
  maxClusterMassCmd->SetGuidance(" Allowed range: [2,12]");
  maxClusterMassCmd->SetParameterName("MaxClusterMass",true);
  maxClusterMassCmd->SetDefaultValue(8);
  maxClusterMassCmd->SetRange("MaxClusterMass>=2 && MaxClusterMass<=12");

}

G4INCLXXInterfaceMessenger::~G4INCLXXInterfaceMessenger() {
  delete theINCLXXDirectory;
  delete accurateNucleusCmd;
  delete maxClusterMassCmd;
}

void G4INCLXXInterfaceMessenger::SetNewValue(G4UIcommand *command, G4String newValues) {
  if(command==accurateNucleusCmd) {
    newValues.toLower();
    if(newValues=="projectile") {
      G4cout << "INCL++ interface: accurate description of the projectile in nucleus-nucleus reactions." << G4endl;
      theINCLXXInterfaceConfig->SetAccurateProjectile(true);
    } else if(newValues=="target") {
      G4cout << "INCL++ interface: accurate description of the target in nucleus-nucleus reactions." << G4endl;
      theINCLXXInterfaceConfig->SetAccurateProjectile(false);
    }
  } else if(command==maxClusterMassCmd) {
    const G4int parameter = maxClusterMassCmd->GetNewIntValue(newValues);
    G4cout << "INCL++ interface: setting maximum cluster mass to " << parameter << "." << G4endl;
    theINCLXXInterfaceConfig->SetMaxClusterMass(parameter);
  }

  theINCLXXInterfaceConfig->DeleteModels();
}
