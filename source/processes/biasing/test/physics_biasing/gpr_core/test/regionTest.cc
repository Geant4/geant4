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
// $Id: regionTest.cc,v 1.1 2007-08-07 22:43:17 tinslay Exp $
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
#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4GRSVolume.hh"
#include "G4GPRTriggerSuperStore.hh"
#include "TestSetup.hh"

using namespace TestSetup;

namespace RegionA {

  G4bool Trigger(const G4Track& track) 
  {
    return (track.GetVolume()->GetLogicalVolume()->GetRegion()->GetName() == "RegionA");
  }

  // Process functions
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute AtRestDoIt_RegionA"<<G4endl;
    return 0;
  }
  
  G4VParticleChange* ContinuousDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute ContinuousDoIt_RegionA"<<G4endl;
    return 0;
  }
  
  G4VParticleChange* DiscreteDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute DiscreteDoIt_RegionA"<<G4endl;
    return 0;
  }
  
  G4double AtRestGPIL(const G4Track& track, G4ForceCondition* condition) 
  {
    G4cout<<"Execute AtRestGPIL_RegionA"<<G4endl;
    return 0;
  }
  
  G4double ContinuousGPIL(const G4Track& track,
			  G4double  previousStepSize,
			  G4double  currentMinimumStep,
			  G4double& proposedSafety,
			  G4GPILSelection* selection)
  {
    G4cout<<"Execute ContinuousGPIL_RegionA"<<G4endl;
    return 0;
  }
  
  G4double DiscreteGPIL(const G4Track& track,
			G4double   previousStepSize,
			G4ForceCondition* condition) {
    G4cout<<"Execute DiscreteGPIL_RegionA"<<G4endl;
    return 0;
  }
}

namespace RegionB {

  G4bool Trigger(const G4Track& track) 
  {
    return (track.GetVolume()->GetLogicalVolume()->GetRegion()->GetName() == "RegionB");
  }
  
  // Process functions
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute AtRestDoIt_RegionB"<<G4endl;
    return 0;
  }
  
  G4VParticleChange* ContinuousDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute ContinuousDoIt_RegionB"<<G4endl;
    return 0;
  }
  
  G4VParticleChange* DiscreteDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute DiscreteDoIt_RegionB"<<G4endl;
    return 0;
  }
  
  G4double AtRestGPIL(const G4Track& track, G4ForceCondition* condition) 
  {
  G4cout<<"Execute AtRestGPIL_RegionB"<<G4endl;
  return 0;
  }
  
  G4double ContinuousGPIL(const G4Track& track,
			  G4double  previousStepSize,
			  G4double  currentMinimumStep,
			     G4double& proposedSafety,
			  G4GPILSelection* selection)
  {
    G4cout<<"Execute ContinuousGPIL_RegionB"<<G4endl;
    return 0;
  }
  
  G4double DiscreteGPIL(const G4Track& track,
			G4double   previousStepSize,
			G4ForceCondition* condition) 
  {
    G4cout<<"Execute DiscreteGPIL_RegionB"<<G4endl;
    return 0;
  }
}
///////////////////////////////

int main(int argc, char** argv) {

  // Construct world
  G4LogicalVolume* logicalWorld = new G4LogicalVolume(new G4Box("World", 2.0*m, 2.0*m, 2.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "World");
  
  G4VPhysicalVolume* world = new G4PVPlacement(0,
                                               G4ThreeVector(),
                                               logicalWorld,
                                               "World",
                                               0,   
                                               false,
                                               0);

  G4LogicalVolume* volA_log = new G4LogicalVolume(new G4Box("World", 1.0*m, 1.0*m, 1.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "volA_log");
  
  G4VPhysicalVolume* volA_phys = new G4PVPlacement(0,
                                               G4ThreeVector(0, 0, -1.0*m),
                                               volA_log,
                                               "volA_phys",
                                               0,   
                                               false,
                                               0);

  G4LogicalVolume* volB_log = new G4LogicalVolume(new G4Box("World", 1.0*m, 1.0*m, 1.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "volB_log");
  
  G4VPhysicalVolume* volB_phys = new G4PVPlacement(0,
                                               G4ThreeVector(0, 0, 1.0*m),
                                               volB_log,
                                               "volB_phys",
                                               0,   
                                               false,
                                               0);

  G4Region* regionA = new G4Region("RegionA");
  volA_log->SetRegion(regionA);

  G4Region* regionB = new G4Region("RegionB");
  volB_log->SetRegion(regionB);

  //  G4Navigator navigator;
  //  navigator.SetWorldVolume(world);

  G4GRSVolume* touchable_A = new G4GRSVolume(volA_phys,NULL,G4ThreeVector(0,0,0));      
  G4GRSVolume* touchable_B = new G4GRSVolume(volB_phys,NULL,G4ThreeVector(0,0,0));      

  G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>* atRestDoIt_A = 
    new G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>("AtRestDoIt_A", &RegionA::AtRestDoIt, G4GPRPlacement::First);

  G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>* atRestDoIt_B = 
    new G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>("AtRestDoIt_B", &RegionB::AtRestDoIt, G4GPRPlacement::First);
  
  G4GPRTriggerSuperStore* triggerSuperStore = G4GPRTriggerSuperStore::Instance();
  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Geometry::NewRegion>::Register(&RegionA::Trigger, atRestDoIt_A, 
										      &G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>::ChangeState);
  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Geometry::NewRegion>::Register(&RegionB::Trigger, atRestDoIt_B, 
										      &G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>::ChangeState);

  // Register element with super store which takes ownership
  G4GPRElementSuperStore* superStore = G4GPRElementSuperStore::Instance();

  superStore->G4GPRManagerT< G4GPRSeedT<G4GPRProcessLists::AtRestDoIt> >::Register(atRestDoIt_A);
  superStore->G4GPRManagerT< G4GPRSeedT<G4GPRProcessLists::AtRestDoIt> >::Register(atRestDoIt_B);
  
  // Process list

  // Dummy track and step
  G4DynamicParticle* dynamic= new G4DynamicParticle(G4Gamma::Gamma(), 0, G4ThreeVector(0, 0, 0));
  G4Track* track = new G4Track(dynamic, 0, G4ThreeVector(0,0,0)); 
  G4Step* step = new G4Step;

  track->SetTouchableHandle(touchable_A);
  
  // In regular processing could keep a pointer to the previous region so know when region has changed.
  // Anyway, manually trigger change in region here.
  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Geometry::NewRegion>::Fire(*track);  

  // Generate list in region A
  typedef std::vector< G4GPRProcessWrappers::Wrappers<G4GPRProcessLists::AtRestDoIt>::SeedWrapper > ProcessList;

  // Use simple generator for the moment - no list caching
  G4GPRSimpleGenerator generator;
  
  ProcessList* result(0);
  generator.Generate<G4GPRProcessLists::AtRestDoIt>(result);
  G4cout<<"jane process list for region A"<<G4endl;
  PrintList(result);

  // Swith to region B
  track->SetTouchableHandle(touchable_B);
  triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Geometry::NewRegion>::Fire(*track); 

  generator.Generate<G4GPRProcessLists::AtRestDoIt>(result);
  G4cout<<"jane process list for region B"<<G4endl;
  PrintList(result);
 
  /*
  generator.Generate<G4GPRProcessLists::DiscreteDoIt>(result);
  
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
  */
}
