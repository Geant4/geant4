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
// $Id: G4GlobalFastSimulationManager_deprecating.cc,v 1.1 2006-11-03 17:33:30 mverderi Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
//  G4GlobalFastSimulationManager_deprecating.cc
//
//  Description:
//    Gathers code to be dropped at the next major release
//
//  History:
//    November 2006
//
//---------------------------------------------------------------

#include "G4GlobalFastSimulationManager_deprecating.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4FastSimulationMessenger.hh"

#include "G4GlobalFastSimulationManager.hh"

G4GlobalFastSimulationManager_deprecating::G4GlobalFastSimulationManager_deprecating(G4GlobalFastSimulationManager* theGlobal)
{
  // It starts closed.
  fClosed=true;
  _global = theGlobal;
}

void G4GlobalFastSimulationManager_deprecating::FastSimulationNeedsToBeClosed()
{
  fClosed=false;
}

void G4GlobalFastSimulationManager_deprecating::CloseFastSimulation()
{
  if(fClosed) return;
  
  G4ParticleTable* theParticleTable;
  
  // Reset the NeededFlavoredWorlds List
  NeededFlavoredWorlds.clearAndDestroy();
  
  // We'll look for models for all particle types
  theParticleTable=G4ParticleTable::GetParticleTable();
  
  // Creates the first world volume clone.
  G4VPhysicalVolume* aClone = GiveMeAWorldVolumeClone();
  
  G4cout << "Closing Fast-Simulation ..." << G4endl;
  
  for (G4int iParticle=0;
       iParticle<theParticleTable->entries(); iParticle++) {
    G4bool Needed = false;
    for (size_t ifsm=0; ifsm< _global->ManagedManagers.size(); ifsm++)
      Needed = Needed || _global->ManagedManagers[ifsm]->
        InsertGhostHereIfNecessary(aClone,
                                   *(theParticleTable->
                                     GetParticle(iParticle)));
    // if some FSM inserted a ghost, keep this clone.
    if(Needed)
      {
	G4Region* aWorldGhostRegion = new G4Region("Ghost world region for " + theParticleTable->GetParticle(iParticle)->GetParticleName());
	aWorldGhostRegion->AddRootLogicalVolume(aClone->GetLogicalVolume());
	NeededFlavoredWorlds.push_back(new G4FlavoredParallelWorld(theParticleTable->GetParticle(iParticle), aClone));
	// and prepare a new one.
	aClone=GiveMeAWorldVolumeClone();
      }
  }
  
  fClosed=true;
}

G4VFlavoredParallelWorld* 
G4GlobalFastSimulationManager_deprecating::
GetFlavoredWorldForThis(G4ParticleDefinition* particle)
{
  for (size_t ipw=0; ipw<NeededFlavoredWorlds.size(); ipw++)
    if(NeededFlavoredWorlds[ipw]->GetTheParticleType()==particle)
      return NeededFlavoredWorlds[ipw];
  return 0;
}


G4VPhysicalVolume* 
G4GlobalFastSimulationManager_deprecating::GiveMeAWorldVolumeClone()
{
  G4VPhysicalVolume* theParallelWorldPhysical;
  
  G4TransportationManager *transportationManager= 
    G4TransportationManager::GetTransportationManager();
  
  G4VSolid *parallelWorldSolid=transportationManager->
    GetNavigatorForTracking()->GetWorldVolume()->
    GetLogicalVolume()->GetSolid();
  
  G4Material *parallelWorldMaterial=transportationManager->
    GetNavigatorForTracking()->GetWorldVolume()->
    GetLogicalVolume()->GetMaterial();
  
  G4LogicalVolume *parallelWorldLog = 
    new G4LogicalVolume(parallelWorldSolid,
                        parallelWorldMaterial,
                        "ParallelWorldLogical",
                        0, 0, 0);
  
  theParallelWorldPhysical= new G4PVPlacement(0,G4ThreeVector(),
                                              "ParallelWorldPhysical",
                                              parallelWorldLog,
                                              0,false,0);
  
  return theParallelWorldPhysical;
}


G4bool 
G4GlobalFastSimulationManager_deprecating::Notify(G4ApplicationState requestedState)
{ 
  G4StateManager * stateManager = G4StateManager::GetStateManager();
  if((stateManager->GetPreviousState()==G4State_Idle) &&
     (requestedState==G4State_GeomClosed) && (!fClosed) )
  {
    G4cerr
      << "ERROR - G4GlobalFastSimulationManager_deprecating::Notify()" << G4endl
      << "        1) you are using ghost volumes;" << G4endl
      << "        2) In this case the G4GlobalFastSimulationManager_deprecating" << G4endl
      << "           MUST be closed BEFORE closing the geometry;" << G4endl
      << "        3) To do this put in your code the call to" << G4endl
      << "           GetGlobalFastSimulationManager()->CloseFastSimulation();"
      << G4endl
      << "           just before closing the geometry." << G4endl;
    G4Exception("G4GlobalFastSimulationManager_deprecating::Notify()", "InvalidState",
                FatalException, "Fast simulation must be closed.");
  }
  return true;
}
