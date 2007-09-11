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
// $Id: singleProcessRelayTest.cc,v 1.4 2007-09-11 03:01:44 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, May 2007. Functor demonstration.
//
#include "G4GPRProcessWrappers.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"

#include "G4GPRSeedT.hh"
#include "G4GPRElementSuperStore.hh"
#include "G4GPRSimpleGenerator.hh"

#include "G4GPRSingleProcessRelayT.hh"

#include "G4GPRPhysicsListManagerSuperStore.hh"
#include "G4GPRNode.hh"
#include "G4GPRManager.hh"
#include "G4Gamma.hh"

using namespace G4GPRProcessWrappers;

// Process function
G4VParticleChange* SeedFunction(const G4Track&, const G4Step&) {
  G4cout<<"Execute SeedFunction"<<G4endl;
  return 0;
}

G4VParticleChange* RelayFunction(G4GPRDoItWrapper& wrappedSeedFunction, const G4Track& track, const G4Step& step) {
  G4cout<<"Execute RelayFunction"<<G4endl;
  G4cout<<"jane "<<wrappedSeedFunction.GetIdentifier();
  wrappedSeedFunction(track, step);
  return 0;
}

int main(int argc, char** argv) {

  typedef G4GPRSeedT<G4GPRProcessLists::DiscreteDoIt> Seed;
  typedef G4GPRSingleProcessRelayT<G4GPRProcessLists::DiscreteDoIt> Relay;

  Seed* seed = new Seed("Seed", &SeedFunction, G4GPRPlacement::First);
  Relay* relay = new Relay("Relay", &RelayFunction, G4GPRPlacement::First);
  
  G4ParticleDefinition* def = G4Gamma::Definition();

  G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];
  G4GPRPhysicsList* physicsList = physicsListManager->GetDefaultList();

  G4GPRElementStore* elementStore = &(*G4GPRElementSuperStore::Instance())[def][physicsList];


  elementStore->G4GPRManagerT<Seed>::Register(seed);
  elementStore->G4GPRManagerT<Relay>::Register(relay);

  // Generate process list
  typedef std::vector< G4GPRProcessWrappers::Wrappers<G4GPRProcessLists::DiscreteDoIt>::SeedWrapper > ProcessList;

  G4Track* dummyTrk = new G4Track; 
  G4Step* dummyStep = new G4Step;

  G4GPRManager gprManager(def);

  ProcessList* result(0);
  gprManager.GetList<G4GPRProcessLists::DiscreteDoIt>(result);
  
  // Iterate over process list
  G4cout<<"jane proc list length "<<result->size()<<G4endl;
  for (ProcessList::iterator iter = result->begin(); iter != result->end(); iter++) {
    G4cout<<"Executing functor :" <<iter->GetIdentifier();
    (*iter)(*dummyTrk, *dummyStep);
    G4cout<<G4endl;
  }

  return 0;
}
