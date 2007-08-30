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
// $Id: seedTest.cc,v 1.3 2007-08-30 19:37:45 tinslay Exp $
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
#include "G4Gamma.hh"
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

// Process function
G4VParticleChange* MyFunction(const G4Track&, const G4Step&) {
  G4cout<<"Execute MyFunction"<<G4endl;
  return 0;
}

int main(int argc, char** argv) {

  G4VProcess* vProcess = new VProcess();
  OtherProcess* otherProcess = new OtherProcess();

  typedef G4GPRSeedT<G4GPRProcessLists::DiscreteDoIt> Element;

  Element* element0 = new Element("Element0", otherProcess, 0);  
  Element* element1 = new Element("Element1", otherProcess, &OtherProcess::Method, 1);
  Element* element2 = new Element("Element2", vProcess, &G4VProcess::PostStepDoIt, 2);
  Element* element3 = new Element("Element3", &MyFunction, 3);
  
  G4ParticleDefinition* def = G4Gamma::Definition();
  G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];
  G4GPRPhysicsList* physicsList = physicsListManager->GetDefaultList();

  G4GPRElementStore* elementStore = &(*G4GPRElementSuperStore::Instance())[def][physicsList];


  elementStore->G4GPRManagerT<Element>::Register(element3);
  elementStore->G4GPRManagerT<Element>::Register(element1);
  elementStore->G4GPRManagerT<Element>::Register(element0);
  elementStore->G4GPRManagerT<Element>::Register(element2);


  // Generate process list
  typedef std::vector< G4GPRProcessWrappers::Wrappers<G4GPRProcessLists::DiscreteDoIt>::SeedWrapper > ProcessList;

  G4GPRManager gprManager(def);
   
  ProcessList* result(0);
  gprManager.GetList<G4GPRProcessLists::DiscreteDoIt>(result);
  
  G4Track* dummyTrk = new G4Track; 
  G4Step* dummyStep = new G4Step;
  G4cout<<"jane result "<<result<<G4endl;
  // Iterate over process list
  for (ProcessList::iterator iter = result->begin(); iter != result->end(); iter++) {
    (*iter)(*dummyTrk, *dummyStep);
  }

  // Cleanup
  delete vProcess;
  delete otherProcess;

  return 0;
}
