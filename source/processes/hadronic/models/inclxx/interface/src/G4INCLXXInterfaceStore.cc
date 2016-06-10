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

/** \file G4INCLXXInterfaceStore.cc
 * \brief The G4INCLXXInterfaceStore class implementation
 *
 * \date 24 May 2012
 * \author Davide Mancusi
 */

#include "G4INCLXXInterfaceStore.hh"
#include "G4INCLXXInterfaceMessenger.hh"
#include "G4INCLConfigEnums.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4INCLXXInterface.hh"
#include "G4INCLConfig.hh"
#include "G4AblaInterface.hh"
#include <vector>

G4ThreadLocal G4INCLXXInterfaceStore *G4INCLXXInterfaceStore::theInstance = NULL;

G4INCLXXInterfaceStore::G4INCLXXInterfaceStore() :
  accurateProjectile(true),
  theMaxProjMassINCL(18),
  cascadeMinEnergyPerNucleon(1.*MeV),
  conservationTolerance(5*MeV),
  theINCLModel(NULL),
  theTally(NULL),
  nWarnings(0),
  maxWarnings(50)
{
  constructINCLXXVersionName();
  theINCLXXInterfaceMessenger = new G4INCLXXInterfaceMessenger(this);
}

G4INCLXXInterfaceStore::~G4INCLXXInterfaceStore() {
  delete theINCLXXInterfaceMessenger;
  delete theINCLModel;
}

void G4INCLXXInterfaceStore::DeleteModel() {
  delete theINCLModel; theINCLModel=NULL;
}

G4INCLXXInterfaceStore *G4INCLXXInterfaceStore::GetInstance() {
  if(!theInstance)
    theInstance = new G4INCLXXInterfaceStore;
  return theInstance;
}

void G4INCLXXInterfaceStore::DeleteInstance() {
  delete theInstance;
  theInstance = NULL;
}

G4INCL::INCL *G4INCLXXInterfaceStore::GetINCLModel() {
  if(!theINCLModel) {
    G4INCL::Config *aConfig = new G4INCL::Config(theConfig);
    theINCLModel = new G4INCL::INCL(aConfig);
    // ownership of the aConfig object is taken over by the INCL model engine
  }
  return theINCLModel;
}

void G4INCLXXInterfaceStore::constructINCLXXVersionName() {
  const std::string versionID = G4INCL_VERSION_ID;
  const size_t lastDash = versionID.find_last_of("-");
  versionName = "INCL++ " + versionID.substr(0,lastDash);
}

const std::string &G4INCLXXInterfaceStore::getINCLXXVersionName() {
  return versionName;
}



void G4INCLXXInterfaceStore::SetAccurateProjectile(const G4bool b) {
  if(accurateProjectile!=b) {
    // Parameter is changed, emit a big warning message
    std::stringstream ss;
    ss << "Switching from "
      << (accurateProjectile ? "\"accurate projectile\" mode to \"accurate target\"" : "\"accurate target\" mode to \"accurate projectile\"")
      << " mode."
      << G4endl
      << "Do this ONLY if you fully understand what it does!";
    EmitBigWarning(ss.str());
  }

  // No need to delete the model for this parameter

  accurateProjectile=b;
}

void G4INCLXXInterfaceStore::SetMaxClusterMass(const G4int aMass) {
  const G4int theMaxClusterMass = theConfig.getClusterMaxMass();
  if(theMaxClusterMass!=aMass) {
    // Parameter is changed, emit a big warning message
    std::stringstream ss;
    ss << "Changing maximum cluster mass from "
      << theMaxClusterMass
      << " to "
      << aMass
      << "."
      << G4endl
      << "Do this ONLY if you fully understand what this setting does!";
    EmitBigWarning(ss.str());

    // We must delete the model object to make sure that we use the new
    // parameter
    DeleteModel();

    theConfig.setClusterMaxMass(aMass);
  }
}




G4bool G4INCLXXInterfaceStore::GetAccurateProjectile() const { return accurateProjectile; }

G4double G4INCLXXInterfaceStore::GetCascadeMinEnergyPerNucleon() const { return cascadeMinEnergyPerNucleon; }

G4INCL::Config &G4INCLXXInterfaceStore::GetINCLConfig() {
  DeleteModel(); // in case the Config is modified
  return theConfig;
}

G4double G4INCLXXInterfaceStore::GetConservationTolerance() const { return conservationTolerance; }




G4int G4INCLXXInterfaceStore::GetMaxProjMassINCL() const { return theMaxProjMassINCL; }

void G4INCLXXInterfaceStore::EmitWarning(const G4String &message) {
  if(++nWarnings<=maxWarnings) {
    G4cout << "[INCL++] Warning: " << message << G4endl;
    if(nWarnings==maxWarnings) {
      G4cout << "[INCL++] INCL++ has already emitted " << maxWarnings << " warnings and will emit no more." << G4endl;
    }
  }
}

void G4INCLXXInterfaceStore::EmitBigWarning(const G4String &message) const {
  G4cout
    << G4endl
    << "================================================================================"
    << G4endl
    << "                                 INCL++ WARNING                                 "
    << G4endl
    << message
    << G4endl
    << "================================================================================"
    << G4endl
    << G4endl;
}

void G4INCLXXInterfaceStore::SetCascadeMinEnergyPerNucleon(const G4double anEnergy) {
  if(cascadeMinEnergyPerNucleon!=anEnergy) {
    // Parameter is changed, emit a big warning message
    std::stringstream ss;
    ss << "Changing minimim cascade energy from "
      << cascadeMinEnergyPerNucleon / MeV
      << " to "
      << anEnergy / MeV
      << " MeV."
      << G4endl
      << "Do this ONLY if you fully understand what this setting does!";
    EmitBigWarning(ss.str());
  }

  // No need to delete the model object

  cascadeMinEnergyPerNucleon=anEnergy;
}

void G4INCLXXInterfaceStore::SetConservationTolerance(const G4double aTolerance) {
  conservationTolerance = aTolerance;
}

G4INCLXXVInterfaceTally *G4INCLXXInterfaceStore::GetTally() const { return theTally; }

void G4INCLXXInterfaceStore::SetTally(G4INCLXXVInterfaceTally * const aTally) { theTally = aTally; }

void G4INCLXXInterfaceStore::SetINCLPhysics(const G4String &option) {
  if(option == "default") {
    theConfig.init();
  } else if(option == "incl42") {
    const G4String message = "Changing INCL++ physics to mimick INCL4.2. Do this ONLY if you fully understand the implications!";
    EmitBigWarning(message);

    theConfig.setPotentialType(G4INCL::ConstantPotential);
    theConfig.setPionPotential(false);
    theConfig.setLocalEnergyBBType(G4INCL::NeverLocalEnergy);
    theConfig.setLocalEnergyPiType(G4INCL::NeverLocalEnergy);
    theConfig.setBackToSpectator(false);
    theConfig.setClusterAlgorithm(G4INCL::NoClusterAlgorithm);
    theConfig.setCoulombType(G4INCL::NoCoulomb);
    // UseRealMasses intentionally left out because it creates problems with
    // energy conservation
    // theConfig.setUseRealMasses(false);
    theConfig.setCrossSectionsType(G4INCL::INCL46CrossSections);
  } else {
    G4Exception("G4INCLXXInterfaceStore::SetINCLPhysics", "INCLXX0001", FatalErrorInArgument,
                "SetINCLPhysics argument must be one of: default, incl42"
                );
  }
}

void G4INCLXXInterfaceStore::UseAblaDeExcitation() {
  // Get hold of pointers to the INCL++ model interfaces
  std::vector<G4HadronicInteraction *> const &interactions = G4HadronicInteractionRegistry::Instance()
    ->FindAllModels(G4INCLXXInterfaceStore::GetInstance()->getINCLXXVersionName());
  for(std::vector<G4HadronicInteraction *>::const_iterator iInter=interactions.begin(), e=interactions.end();
      iInter!=e; ++iInter) {
    G4INCLXXInterface *theINCLInterface = dynamic_cast<G4INCLXXInterface*>(*iInter);
    if(theINCLInterface) {
      // Instantiate the ABLA model
      G4HadronicInteraction *interaction = G4HadronicInteractionRegistry::Instance()->FindModel("ABLA");
      G4AblaInterface *theAblaInterface = dynamic_cast<G4AblaInterface*>(interaction);
      if(!theAblaInterface)
        theAblaInterface = new G4AblaInterface;
      // Couple INCL++ to ABLA
      G4cout << "Coupling INCLXX to ABLA" << G4endl;
      theINCLInterface->SetDeExcitation(theAblaInterface);
    }
  }
}
