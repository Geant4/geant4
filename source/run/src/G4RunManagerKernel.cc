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
//
// $Id: G4RunManagerKernel.cc,v 1.34 2006/06/29 21:13:52 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//

#include "G4RunManagerKernel.hh"

#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "G4ExceptionHandler.hh"
#include "G4PrimaryTransformer.hh"
#include "G4GeometryManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4Version.hh"
#include "G4ios.hh"
#include <vector>

G4RunManagerKernel* G4RunManagerKernel::fRunManagerKernel = 0;

G4RunManagerKernel* G4RunManagerKernel::GetRunManagerKernel()
{ return fRunManagerKernel; }

G4RunManagerKernel::G4RunManagerKernel()
:physicsList(0),currentWorld(0),
 geometryInitialized(false),physicsInitialized(false),
 geometryNeedsToBeClosed(true),geometryToBeOptimized(true),
 physicsNeedsToBeReBuilt(true),verboseLevel(0)
{
  defaultExceptionHandler = new G4ExceptionHandler();
  if(fRunManagerKernel)
  {
    G4Exception("G4RunManagerKernel::G4RunManagerKernel()","MoreThanOneRunManager",
                FatalException,"More than one G4RunManagerKernel is constructed.");
  }
  fRunManagerKernel = this;

  // construction of Geant4 kernel classes
  eventManager = new G4EventManager();
  defaultRegion = new G4Region("DefaultRegionForTheWorld");
  defaultRegion->SetProductionCuts(
    G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());

  // Following line is tentatively moved from SetPhysics method
  // Commented out for introduction of non-static particle definition // G4ParticleTable::GetParticleTable()->SetReadiness();
  // set the initial application state
  G4StateManager::GetStateManager()->SetNewState(G4State_PreInit);

  // version banner
  G4String vs = G4Version;
  vs = vs.substr(1,vs.size()-2);
  versionString = " Geant4 version ";
  versionString += vs;
  versionString += "   ";
  versionString += G4Date;
  G4cout << G4endl
    << "*************************************************************" << G4endl
    << versionString << G4endl
    << "                      Copyright : Geant4 Collaboration" << G4endl
    << "                      Reference : NIM A 506 (2003), 250-303" << G4endl
    << "                            WWW : http://cern.ch/geant4" << G4endl
    << "*************************************************************" << G4endl
    << G4endl;
}

G4RunManagerKernel::~G4RunManagerKernel()
{
  // set the application state to the quite state
  if(verboseLevel>0) G4cout << "G4 kernel has come to Quit state." << G4endl;
  G4StateManager* pStateManager = G4StateManager::GetStateManager();
  pStateManager->SetNewState(G4State_Quit);

  // open geometry for deletion
  G4GeometryManager::GetInstance()->OpenGeometry();

  // deletion of Geant4 kernel classes
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
  {
    delete fSDM;
    if(verboseLevel>1) G4cout << "G4SDManager deleted." << G4endl;
  }
  delete eventManager;
  if(verboseLevel>1) G4cout << "EventManager deleted." << G4endl;
  G4UImanager* pUImanager = G4UImanager::GetUIpointer();
  {
    if(pUImanager) delete pUImanager;
    if(verboseLevel>1) G4cout << "UImanager deleted." << G4endl;
  }
  delete pStateManager; 
  if(verboseLevel>1) G4cout << "StateManager deleted." << G4endl;
  delete defaultExceptionHandler;
  if(verboseLevel>1) G4cout << "RunManagerKernel is deleted." << G4endl;
  fRunManagerKernel = 0;
}

void G4RunManagerKernel::DefineWorldVolume(G4VPhysicalVolume* worldVol,
                                     G4bool topologyIsChanged)
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(!(currentState==G4State_Idle||currentState==G4State_PreInit))
  { 
    G4Exception("G4RunManagerKernel::DefineWorldVolume",
                "DefineWorldVolumeAtIncorrectState",
                JustWarning,
                "Geant4 kernel is not PreInit or Idle state : Method ignored.");
    if(verboseLevel>1) G4cerr << "Current application state is " 
      << stateManager->GetStateString(currentState) << G4endl;
    return;
  }

  // The world volume MUST NOT have a region defined by the user
  if(worldVol->GetLogicalVolume()->GetRegion())
  {
    if(worldVol->GetLogicalVolume()->GetRegion()!=defaultRegion)
    {
      G4cerr << "The world volume has a user-defined region <"
           << worldVol->GetLogicalVolume()->GetRegion()->GetName()
           << ">." << G4endl;
      G4Exception("G4RunManager::DefineWorldVolume",
                "RUN:WorldHasUserDefinedRegion",
                FatalException,
                "World would have a default region assigned by RunManagerKernel.");
    }
  }

  // Remove old world logical volume from the default region, if exist
  if(defaultRegion->GetNumberOfRootVolumes())
  {
    if(defaultRegion->GetNumberOfRootVolumes()>size_t(1))
    {
      G4Exception("G4RunManager::DefineWorldVolume",
                "DefaultRegionHasMoreThanOneVolume",
                FatalException,
                "Default world region should have a unique logical volume.");
    }
    std::vector<G4LogicalVolume*>::iterator lvItr
     = defaultRegion->GetRootLogicalVolumeIterator();
    defaultRegion->RemoveRootLogicalVolume(*lvItr);
    if(verboseLevel>1) G4cout << (*lvItr)->GetName()
     << " is removed from the default region." << G4endl;
  }

  // Accept the world volume
  currentWorld = worldVol; 

  // Set the default region to the world
  G4LogicalVolume* worldLog = currentWorld->GetLogicalVolume();
  worldLog->SetRegion(defaultRegion);
  defaultRegion->AddRootLogicalVolume(worldLog);
  if(verboseLevel>1) G4cout << worldLog->GetName()
   << " is registered to the default region." << G4endl;

  // Set the world volume, notify the Navigator and reset its state
  G4TransportationManager::GetTransportationManager()
      ->GetNavigatorForTracking()
      ->SetWorldVolume(currentWorld); 
  if(topologyIsChanged) geometryNeedsToBeClosed = true;
  
  // Notify the VisManager as well
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) pVVisManager->GeometryHasChanged();

  geometryInitialized = true;
  if(physicsInitialized && currentState!=G4State_Idle)
  { stateManager->SetNewState(G4State_Idle); }
} 
  
void G4RunManagerKernel::SetPhysics(G4VUserPhysicsList* uPhys)
{
  physicsList = uPhys;
  G4ParticleTable::GetParticleTable()->SetReadiness();
  physicsList->ConstructParticle();
  if(verboseLevel>2) G4ParticleTable::GetParticleTable()->DumpTable();
  if(verboseLevel>1)
  {
    G4cout << "List of instantiated particles ============================================" << G4endl;
    G4int nPtcl = G4ParticleTable::GetParticleTable()->entries();
    for(G4int i=0;i<nPtcl;i++)
    {
      G4ParticleDefinition* pd = G4ParticleTable::GetParticleTable()->GetParticle(i);
      G4cout << pd->GetParticleName() << " ";
      if(i%10==9) G4cout << G4endl;
    }
    G4cout << G4endl;
  }
}
  
void G4RunManagerKernel::InitializePhysics()
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(!(currentState==G4State_Idle||currentState==G4State_PreInit))
  { 
    G4Exception("G4RunManagerKernel::InitializePhysics",
                "InitializePhysicsAtIncorrectState",
                JustWarning,
                "Geant4 kernel is not PreInit or Idle state : Method ignored.");
    return;
  }

  if(!physicsList)
  {
    G4Exception("G4RunManagerKernel::InitializePhysics",
                "PhysicsListIsNotDefined",
                FatalException,
                "G4VUserPhysicsList is not defined");
  }

  //if(verboseLevel>1) G4cout << "physicsList->ConstructParticle() start." << G4endl;
  //physicsList->ConstructParticle();

  if(verboseLevel>1) G4cout << "physicsList->Construct() start." << G4endl;
  physicsList->Construct();

  if(verboseLevel>1) G4cout << "physicsList->setCut() start." << G4endl;
  physicsList->SetCuts();
  CheckRegions();
  physicsInitialized = true;
  if(geometryInitialized && currentState!=G4State_Idle)
  { stateManager->SetNewState(G4State_Idle); }
}

G4bool G4RunManagerKernel::RunInitialization()
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();

  if(!geometryInitialized) 
  { 
    G4Exception("G4RunManagerKernel::RunInitialization",
                "GeometryHasNotYetInitialized",
                JustWarning,
                "Geometry has not yet initialized : method ignored.");
    return false;
  }
  
  if(!physicsInitialized) 
  { 
    G4Exception("G4RunManagerKernel::RunInitialization",
                "PhysicsHasNotYetInitialized",
                JustWarning,
                "Physics has not yet initialized : method ignored.");
    return false;
  }

  if( currentState != G4State_Idle )
  { 
    G4Exception("G4RunManagerKernel::RunInitialization",
                "RunInitializationAtIncorrectState",
                JustWarning,
                "Geant4 kernel not in Idle state : method ignored.");
    return false;
  }

  CheckRegions();
  UpdateRegion();
  BuildPhysicsTables();

  if(geometryNeedsToBeClosed) ResetNavigator();
 
  GetPrimaryTransformer()->CheckUnknown();

  stateManager->SetNewState(G4State_GeomClosed);
  return true;
}

void G4RunManagerKernel::RunTermination()
{ G4StateManager::GetStateManager()->SetNewState(G4State_Idle); }

void G4RunManagerKernel::ResetNavigator()
{
  // We have to tweak the navigator's state in case a geometry has been
  // modified between runs. By the following calls we ensure that navigator's
  // state is reset properly. It is required the geometry to be closed
  // and previous optimisations to be cleared.
  
  G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
  if(verboseLevel>1) G4cout << "Start closing geometry." << G4endl;
  geomManager->OpenGeometry();
  geomManager->CloseGeometry(geometryToBeOptimized, verboseLevel>1);
 
  // Reseting Navigator has been moved to G4Eventmanager, so that resetting
  // is now done for every event.  
  // G4ThreeVector center(0,0,0);
  // G4Navigator* navigator =
  //     G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  // navigator->LocateGlobalPointAndSetup(center,0,false);

  geometryNeedsToBeClosed = false;
}

void G4RunManagerKernel::UpdateRegion()
{
  G4RegionStore::GetInstance()->UpdateMaterialList(currentWorld);
  G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable(currentWorld);
}

void G4RunManagerKernel::BuildPhysicsTables()
{ 
  if(G4ProductionCutsTable::GetProductionCutsTable()->IsModified()
  || physicsNeedsToBeReBuilt)
  {
    physicsList->BuildPhysicsTable();
    G4ProductionCutsTable::GetProductionCutsTable()->PhysicsTableUpdated();
    physicsNeedsToBeReBuilt = false;
  }

  if(verboseLevel>1) DumpRegion();
  if(verboseLevel>0) physicsList->DumpCutValuesTable();
  physicsList->DumpCutValuesTableIfRequested();
}

void G4RunManagerKernel::CheckRegions()
{
  G4RegionStore::GetInstance()->SetWorldVolume();

  for(size_t i=0;i<G4RegionStore::GetInstance()->size();i++)
  { 
    G4Region* region = (*(G4RegionStore::GetInstance()))[i];
    if(region->GetWorldPhysical()!=currentWorld) continue;
    G4ProductionCuts* cuts = region->GetProductionCuts();
    if(!cuts)
    {
      G4cerr << "Warning : Region <" << region->GetName()
             << "> does not have specific production cuts," << G4endl
             << "even though it appears in the current tracking world." << G4endl;
      G4cerr << "Default cuts are used for this region." << G4endl;
      region->SetProductionCuts(
          G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
    }
  }
}

void G4RunManagerKernel::DumpRegion(G4String rname) const
{
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(rname);
  if(region) DumpRegion(region);
}

void G4RunManagerKernel::DumpRegion(G4Region* region) const
{
  if(!region)
  {
    for(size_t i=0;i<G4RegionStore::GetInstance()->size();i++)
    { DumpRegion((*(G4RegionStore::GetInstance()))[i]); }
  }
  else
  {
    G4cout << G4endl;
    G4cout << "Region <" << region->GetName() << "> -- appears in <" 
           << region->GetWorldPhysical()->GetName() << "> world volume" << G4endl;
    if(region->GetWorldPhysical()!=currentWorld)
    {
      G4cout << G4endl;
      return;
    }
    G4cout << " Materials : ";
    std::vector<G4Material*>::const_iterator mItr = region->GetMaterialIterator();
    size_t nMaterial = region->GetNumberOfMaterials();
    for(size_t iMate=0;iMate<nMaterial;iMate++)
    {
      G4cout << (*mItr)->GetName() << " ";
      mItr++;
    }
    G4cout << G4endl;
    G4ProductionCuts* cuts = region->GetProductionCuts();
    if(!cuts)
    {
      G4cerr << "Warning : Region <" << region->GetName()
             << "> does not have specific production cuts." << G4endl;
      G4cerr << "Default cuts are used for this region." << G4endl;
      region->SetProductionCuts(
          G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
    }
    G4cout << " Production cuts : "
           << " gamma " << G4BestUnit(cuts->GetProductionCut("gamma"),"Length")
           << "    e- " << G4BestUnit(cuts->GetProductionCut("e-"),"Length")
           << "    e+ " << G4BestUnit(cuts->GetProductionCut("e+"),"Length")
           << G4endl;
  }
}
