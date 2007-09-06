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
// $Id: triggerTest.cc,v 1.4 2007-09-06 22:10:10 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#include "G4GPRProcessWrappers.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"

#include "G4GPRSeedT.hh"
#include "G4GPRElementSuperStore.hh"
#include "G4GPRSimpleGenerator.hh"
#include "G4GPRTriggerSuperStore.hh"
#include "G4Gamma.hh"
#include "G4DynamicParticle.hh"
#include "G4GPRPhysicsListManagerSuperStore.hh"
#include "G4GPRNode.hh"
#include "G4GPRManager.hh"

// Regular G4VProcess
struct VProcess : public G4VDiscreteProcess
{
  VProcess():G4VDiscreteProcess("test"){}
  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*){return 0;}

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) {
    G4cout<<"Execute VProcess::PostStepDoIt"<<G4endl;
    return 0;
  }
};

// Alternative process
struct OtherProcess
{
  G4VParticleChange* Method(const G4Track&, const G4Step&)
  {
    G4cout<<"Execute OtherProcess::Method"<<G4endl;
  }

  G4VParticleChange* operator()(const G4Track&, const G4Step&)
  {
    G4cout<<"Execute OtherProcess::operator"<<G4endl;
    return 0;
  }
};

G4bool PrimaryTrackTrigger(G4Track* track) 
{
  G4cout<<"jane executing MyTrigger "<<track->GetTrackID()<<G4endl;
  return (track->GetTrackID() == 0 ? true : false);
}

class MaxEnergyTrigger {

public:
  MaxEnergyTrigger()
    :fMaxEnergy(0)
  {}

  void SetMaxEnergy(G4double maxEnergy) {fMaxEnergy = maxEnergy;}
  G4double GetMaxEnergy() {return fMaxEnergy;}

  G4bool operator()(const G4Track& track, const G4Step& step) {
    G4cout<<"jane max energy trigger "<<track.GetKineticEnergy()<<G4endl;
    return (track.GetKineticEnergy() < fMaxEnergy ? true : false);
  }

  G4bool operator()(const G4Track& track) {
      G4cout<<"jane track energy "<<track.GetKineticEnergy()<<G4endl;
      //return (track.GetKineticEnergy() < fMaxEnergy ? true : false);
  }

private:

  G4double fMaxEnergy;
};

// Process function
G4VParticleChange* MyFunction(const G4Track&, const G4Step&) {
  G4cout<<"Execute MyFunction"<<G4endl;
  return 0;
}

int main(int argc, char** argv) {

  G4VProcess* vProcess = new VProcess();
  OtherProcess* otherProcess = new OtherProcess();

  typedef G4GPRSeedT<G4GPRProcessLists::DiscreteDoIt> Element;

  Element* element0 = new Element("Element0", otherProcess, G4GPRPlacement::First);  
  Element* element1 = new Element("Element1", otherProcess, &OtherProcess::Method, 1);
  Element* element2 = new Element("Element2", vProcess, &G4VProcess::PostStepDoIt, G4GPRPlacement::Append);
  Element* element3 = new Element("Element3", &MyFunction, G4GPRPlacement::Last);

  G4ParticleDefinition* def = G4Gamma::Definition();

  G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];
  G4GPRPhysicsList* physicsList = physicsListManager->GetDefaultList();

  G4GPRElementStore* elementStore = &(*G4GPRElementSuperStore::Instance())[def][physicsList];

  // Random register to test placement
  elementStore->G4GPRManagerT<Element>::Register(element1);
  elementStore->G4GPRManagerT<Element>::Register(element3);  
  elementStore->G4GPRManagerT<Element>::Register(element0);
  elementStore->G4GPRManagerT<Element>::Register(element2);

  G4GPRTriggerStore* triggerStore = &(*G4GPRTriggerSuperStore::Instance())[def][physicsList];

  // Elements 1 and 2 are on the same trigger. They're only active 
  // when tracking a primary (track id = 0). Should be evaluated
  // when starting to track a track.
  triggerStore->G4GPRTriggerManagerT<G4GPRTriggering::Tracking::StartTracking>::Register(&PrimaryTrackTrigger, element0, &Element::ChangeState);
  triggerStore->G4GPRTriggerManagerT<G4GPRTriggering::Tracking::StartTracking>::Register(&PrimaryTrackTrigger, element1, &Element::ChangeState);

  // Elements 3 and 4 are on the same trigger. They only become active when track energy falls below 50 MeV. 
  // Should be evaluated at the start of each step.
  // Restriction - want to be able to have multiple elements associated with
  // one trigger so that don't have to evaluate the same trigger multiple
  // times. However, since can't meaningfully compare virtual member function pointers, 
  // trigger classes must provide appropriate "operator()" method rather than any old method for the
  // trigger.
  MaxEnergyTrigger* maxEnergyTrigger = new MaxEnergyTrigger;
  maxEnergyTrigger->SetMaxEnergy(50*MeV);
  triggerStore->G4GPRTriggerManagerT<G4GPRTriggering::Stepping::StartStep>::Register(maxEnergyTrigger, element2, &Element::ChangeState);
  triggerStore->G4GPRTriggerManagerT<G4GPRTriggering::Stepping::StartStep>::Register(maxEnergyTrigger, element3, &Element::ChangeState);

  G4GPRNode* node1 = new G4GPRNode;
  G4GPRNode* node2 = new G4GPRNode;

  triggerStore->G4GPRTriggerManagerT<G4GPRTriggering::Tracking::StartTracking>::Register(&PrimaryTrackTrigger, node1, &G4GPRNode::FlipState);
  triggerStore->G4GPRTriggerManagerT<G4GPRTriggering::Stepping::StartStep>::Register(maxEnergyTrigger, node2, &G4GPRNode::FlipState);

  G4GPRKeyStore* keyStore = &(*G4GPRKeySuperStore::Instance())[def][physicsList];

  keyStore->G4GPRKeyManagerT<G4GPRProcessLists::DiscreteDoIt>::AddNode(node1);
  keyStore->G4GPRKeyManagerT<G4GPRProcessLists::DiscreteDoIt>::AddNode(node2);

  // Pretend to process two tracks - a primary and a secondary.
  // Test with simple generator : *final lists aren't cached*
  typedef std::vector< G4GPRProcessWrappers::Wrappers<G4GPRProcessLists::DiscreteDoIt>::SeedWrapper > ProcessList;
  ProcessList* result(0);

  G4DynamicParticle* dynamic= new G4DynamicParticle(def, 0, G4ThreeVector(0, 0, 0));
  G4Track* trk = new G4Track(dynamic, 0, G4ThreeVector(0,0,0)); 
  trk->SetKineticEnergy(99*MeV);

  G4Step* step = new G4Step;

  G4GPRManager gprManager(def);

  // Start tracking primary
  gprManager.Fire<G4GPRTriggering::Tracking::StartTracking>(trk);  

  // Start Step
  gprManager.Fire<G4GPRTriggering::Stepping::StartStep>(*trk, *step);  
  
  // Generate process list
  gprManager.GetList<G4GPRProcessLists::DiscreteDoIt>(result);
  
  G4cout<<"Jane, list should be : "<<element0->GetName()<<":"<<element1->GetName()<<G4endl;

  // Iterate over process list
  for (ProcessList::iterator iter = result->begin(); iter != result->end(); iter++) {
    G4cout<<"jane "<<iter->GetIdentifier()<<G4endl;
    (*iter)(*trk, *step);
  }
  
  // Track energy decreased by 50MeV in this step
  trk->SetKineticEnergy(49*MeV);

  // Start new step
  gprManager.Fire<G4GPRTriggering::Stepping::StartStep>(*trk, *step);  

  // Generate new list
  gprManager.GetList<G4GPRProcessLists::DiscreteDoIt>(result);

  G4cout<<"jane, list should be : "<<element0->GetName()<<":"<<element1->GetName()<<":"
	<<element2->GetName()<<":"<<element3->GetName()<<G4endl;

  // Iterate over process list
  for (ProcessList::iterator iter = result->begin(); iter != result->end(); iter++) {
    G4cout<<"jane "<<iter->GetIdentifier()<<G4endl;
    (*iter)(*trk, *step);
  }

  // Pretend done processing primary track, and move onto processing secondary with energy 30MeV
  trk->SetTrackID(1);
  trk->SetKineticEnergy(30*MeV);

  // Start tracking secondary 
  gprManager.Fire<G4GPRTriggering::Tracking::StartTracking>(trk);  

  // Start Step
  gprManager.Fire<G4GPRTriggering::Stepping::StartStep>(*trk, *step);  
  
  // Generate list
  gprManager.GetList<G4GPRProcessLists::DiscreteDoIt>(result);

  G4cout<<"Jane, list should be : "<<element2->GetName()<<":"<<element3->GetName()<<G4endl;

  // Iterate over process list
  for (ProcessList::iterator iter = result->begin(); iter != result->end(); iter++) {
    G4cout<<"jane "<<iter->GetIdentifier()<<G4endl;
    (*iter)(*trk, *step);
  }

  // Cleanup
  //  delete vProcess;
  //  delete otherProcess;

  return 0;
}
