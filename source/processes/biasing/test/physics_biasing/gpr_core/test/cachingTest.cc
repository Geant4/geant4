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
// $Id: cachingTest.cc,v 1.1 2007-08-02 18:12:06 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#include "G4GPRProcessWrappers.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
#include "G4GPRSeedT.hh"
#include "G4GPRElementSuperStore.hh"
#include "G4GPRGenerator.hh"
#include "G4GPRTriggerSuperStore.hh"
#include "G4Gamma.hh"
#include "G4DynamicParticle.hh"
#include "G4GPRKeyManagerT.hh"
#include "G4GPRKeyNode.hh"
#include "G4GPRKeySuperStore.hh"
#include "TestSetup.hh"

using namespace TestSetup;

int main(int argc, char** argv) {

  typedef G4GPRSeedT<G4GPRProcessLists::DiscreteDoIt> Element;

  G4VProcess* vProcess = new VProcess();
  OtherProcess* otherProcess = new OtherProcess();

  Element* element0 = new Element("Element0", otherProcess, G4GPRPlacement::First);  
  Element* element1 = new Element("Element1", otherProcess, &OtherProcess::Method, 1);
  Element* element2 = new Element("Element2", vProcess, &G4VProcess::PostStepDoIt, G4GPRPlacement::Append);
  Element* element3 = new Element("Element3", &MyFunction, G4GPRPlacement::Last);

  // Register element with super store which takes ownership
  G4GPRElementSuperStore* elementSuperStore = G4GPRElementSuperStore::Instance();

  // Random register to test placement
  elementSuperStore->G4GPRManagerT<Element>::Register(element1);
  elementSuperStore->G4GPRManagerT<Element>::Register(element3);  
  elementSuperStore->G4GPRManagerT<Element>::Register(element0);
  elementSuperStore->G4GPRManagerT<Element>::Register(element2);

  G4GPRTriggerSuperStore* triggerSuperStore = G4GPRTriggerSuperStore::Instance();

  // Elements 1 and 2 are on the same trigger.
  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Tracking::StartTracking>::Register(&PrimaryTrackTrigger, element0, &Element::ChangeState);
  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Tracking::StartTracking>::Register(&PrimaryTrackTrigger, element1, &Element::ChangeState);

  // Elements 3 and 4 are on the same trigger.
  MaxEnergyTrigger* maxEnergyTrigger = new MaxEnergyTrigger;
  maxEnergyTrigger->SetMaxEnergy(50*MeV);
  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Stepping::StartStep>::Register(maxEnergyTrigger, element2, &Element::ChangeState);
  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Stepping::StartStep>::Register(maxEnergyTrigger, element3, &Element::ChangeState);
 
  // Create and register key nodes with trigger manager so that know when an element has been activated or deactivated
  G4GPRKeyNode* node1 = new G4GPRKeyNode;
  G4GPRKeyNode* node2 = new G4GPRKeyNode;

  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Tracking::StartTracking>::Register(&PrimaryTrackTrigger, node1, &G4GPRKeyNode::ChangeState);

  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Stepping::StartStep>::Register(maxEnergyTrigger, node2, &G4GPRKeyNode::ChangeState);

  G4GPRKeySuperStore* keySuperStore = G4GPRKeySuperStore::Instance();
  keySuperStore->G4GPRKeyManagerT<Element::List>::AddNode(node1);
  keySuperStore->G4GPRKeyManagerT<Element::List>::AddNode(node2);

  // Dummy track and step
  G4DynamicParticle* dynamic= new G4DynamicParticle(G4Gamma::Gamma(), 0, G4ThreeVector(0, 0, 0));
  G4Track* dummyTrk = new G4Track(dynamic, 0, G4ThreeVector(0,0,0)); 
  G4Step* dummyStep = new G4Step;

  G4GPRGenerator generator;
  ProcessList* listA(0);
  ProcessList* listB(0);
  ProcessList* listC(0);

  // Generate first lists
  GenerateList(generator, &ConditionsA, listA, dummyTrk, dummyStep);
  GenerateList(generator, &ConditionsB, listB, dummyTrk, dummyStep);
  GenerateList(generator, &ConditionsC, listC, dummyTrk, dummyStep);

  // Generated all lists, now test caching
  ProcessList* cachedListA(0);
  ProcessList* cachedListB(0);
  ProcessList* cachedListC(0);

  GenerateList(generator, &ConditionsC, cachedListC, dummyTrk, dummyStep);
  GenerateList(generator, &ConditionsA, cachedListA, dummyTrk, dummyStep);
  GenerateList(generator, &ConditionsB, cachedListB, dummyTrk, dummyStep);

  G4cout<<"jane list ptrs A: "<<listA<<" "<<cachedListA<<G4endl;
  G4cout<<"jane list ptrs B: "<<listB<<" "<<cachedListB<<G4endl;
  G4cout<<"jane list ptrs C: "<<listC<<" "<<cachedListC<<G4endl;

  // Caching works if pointers are the same
  assert (listA == cachedListA);
  assert (listB == cachedListB);
  assert (listC == cachedListC);

  // Print lists for good measure
  G4cout<<"jane printing       listA: "; PrintList(listA); G4cout<<G4endl;
  G4cout<<"jane printing cachedListA: "; PrintList(cachedListA); G4cout<<G4endl;

  G4cout<<"jane printing       listB: "; PrintList(listB); G4cout<<G4endl;
  G4cout<<"jane printing cachedListB: "; PrintList(cachedListB); G4cout<<G4endl;

  G4cout<<"jane printing       listC: "; PrintList(listC); G4cout<<G4endl;
  G4cout<<"jane printing cachedListC: "; PrintList(cachedListC); G4cout<<G4endl;

  return 0;
}
