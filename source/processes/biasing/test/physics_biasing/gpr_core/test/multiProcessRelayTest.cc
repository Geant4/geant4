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
// $Id: multiProcessRelayTest.cc,v 1.1 2007-08-10 22:23:04 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Functor demonstration.
//
#include "G4GPRProcessWrappers.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"

#include "G4GPRSeedT.hh"
#include "G4GPRElementSuperStore.hh"
#include "G4GPRSimpleGenerator.hh"

#include "G4GPRMultiProcessRelayT.hh"

// Process function
G4double DiscreteGPIL1(const G4Track& track,
		       G4double previousStepSize,
		       G4ForceCondition* condition) 
{
  
  return 1.0;
}

G4double DiscreteGPIL2(const G4Track& track,
		       G4double previousStepSize,
		       G4ForceCondition* condition) 
{
  
  return 2.0;
}

G4double MyMultiRelay(std::vector<G4DiscreteGPILWrapper>& input, const G4Track& track,
		      G4double previousStepSize,
		      G4ForceCondition* condition)
{
  typedef std::vector<G4DiscreteGPILWrapper> Vect;

  G4double sum = 0;
  G4cout<<"jane input length "<<input.size()<<G4endl;

  for (Vect::iterator iter = input.begin(); iter != input.end(); ++iter) {
    sum += (*iter)(track, previousStepSize, condition);
  }
  
  return sum;
}

int main(int argc, char** argv) {

  typedef G4GPRSeedT<G4GPRProcessLists::DiscreteGPIL> Seed;
  typedef G4GPRMultiProcessRelayT<G4GPRProcessLists::DiscreteGPIL> MultiRelay;

  Seed* seed1 = new Seed("Seed1", &DiscreteGPIL1, G4GPRPlacement::First);
  Seed* seed2 = new Seed("Seed1", &DiscreteGPIL2, G4GPRPlacement::Second);

  std::vector<unsigned> processes;
  processes.push_back(G4GPRPlacement::First);
  processes.push_back(G4GPRPlacement::Second);

  MultiRelay* relay = new MultiRelay("MultiRelay", &MyMultiRelay, processes, G4GPRPlacement::First);
  
  G4GPRElementSuperStore* superStore = G4GPRElementSuperStore::Instance();

  superStore->G4GPRManagerT<Seed>::Register(seed1);
  superStore->G4GPRManagerT<Seed>::Register(seed2);

  superStore->G4GPRManagerT<MultiRelay>::Register(relay);
  
  // Generate process list
  typedef std::vector< G4DiscreteGPILWrapper > ProcessList;

  G4Track* dummyTrk = new G4Track; 
  G4Step* dummyStep = new G4Step;

  G4GPRSimpleGenerator generator;  

  ProcessList* result(0);
  generator.Generate<G4GPRProcessLists::DiscreteGPIL>(result);
  G4cout<<"jane generated size "<<result->size()<<G4endl;
  // Iterate over process list
  for (ProcessList::iterator iter = result->begin(); iter != result->end(); ++iter) {
    G4cout<<"Executing functor :" <<iter->GetIdentifier()<<" : ";
    G4cout<< (*iter)(*dummyTrk, 0, 0);
    G4cout<<G4endl;
  }

  return 0;
}
