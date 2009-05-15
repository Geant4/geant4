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
// $Id: G4tgbParallelGeomMgr.cc,v 1.1 2009-05-15 16:25:47 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4tgbParallelGeomMgr

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbParallelGeomMgr.hh"
#include "G4tgbParallelWorld.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgbVolume.hh"

#include "G4VUserParallelWorld.hh"
#include "G4tgbDetectorBuilder.hh"
#include "G4UIcommand.hh"
#include "G4tgbVParallelWorldCreator.hh"
#include "G4tgbParallelWorldCreator.hh"

G4tgbParallelGeomMgr* G4tgbParallelGeomMgr::theInstance = 0;

//-------------------------------------------------------------
G4tgbParallelGeomMgr::G4tgbParallelGeomMgr()
{
  theWorldCreator = new G4tgbParallelWorldCreator;
}


//-------------------------------------------------------------
G4tgbParallelGeomMgr::~G4tgbParallelGeomMgr()
{ 
  delete theInstance;
}


//-------------------------------------------------------------
G4tgbParallelGeomMgr* G4tgbParallelGeomMgr::GetInstance()
{
  if( !theInstance )
  {
    theInstance = new G4tgbParallelGeomMgr;
  }
  return theInstance;
}


//-------------------------------------------------------------------
void G4tgbParallelGeomMgr::AddParallelWorldIndex( G4int index )
{
  theIndices.insert(index);
}


//-------------------------------------------------------------------
std::vector<G4VUserParallelWorld*> G4tgbParallelGeomMgr::CreateParalleWorlds()
{
  //  G4cout << "  G4tgbParallelGeomMgr::CreateParalleWorlds() nindex " << theIndices.size() << G4endl;
  std::set<G4int>::iterator ite;
  G4tgbVolumeMgr* tgbVolmgr = G4tgbVolumeMgr::GetInstance();
  //  tgbVolmgr->GetDetectorBuilder();

  //Get mas world name
  G4String massWorldName = G4tgrVolumeMgr::GetInstance()->GetTopVolume(-1)->GetName();
  G4tgrVolumeMgr* tgrVolmgr = G4tgrVolumeMgr::GetInstance();
  for( ite = theIndices.begin(); ite != theIndices.end(); ite++ ){
    G4int index = *ite;
    const G4tgrVolume* tgrMassWorld = tgrVolmgr->GetTopVolume(index);  
    if( tgrMassWorld->GetName() != massWorldName ) {
      G4Exception("G4tgbParallelGeomMgr::CreateParalleWorlds()",
		  "Volumes to be placed in parallel world should be placed in the mass world (parallel world will be built automatically",
		  FatalException,
		  ("Top parallel world is "+tgrMassWorld->GetName()
		  +" , mass world is "+ massWorldName).c_str());
    }
    G4String parallelWorldName = massWorldName+"_parallel_"+G4UIcommand::ConvertToString(index);

    G4VUserParallelWorld* parallelWorld = theWorldCreator->CreateParallelWorld( parallelWorldName, index );
    theParallelWorlds.push_back( parallelWorld );
    //     G4cout << " new parallelWorld " << parallelWorld << " " << theParallelWorlds.size() << G4endl;
    //check that there are not two 

    //--- Change the placement in the world to be placement in this parallel world: the parallel world is constructed by G4TransportationManager after calling G4VUserParallalWorld::GetWorld() (there is not other way to do it)
    std::vector<G4tgrPlace*> placeList = tgrVolmgr->GetDetPlaceList();
    std::vector<G4tgrPlace*>::iterator itep;
    for( itep = placeList.begin(); itep != placeList.end(); itep++ ){
      G4tgrPlace* place = *itep;
      if( place->GetParallelID() == index && 
	  place->GetParentName() == massWorldName ) {
	place->SetParentName( parallelWorldName );
	tgrVolmgr->RegisterMe( place );
	tgrVolmgr->RegisterParentChild( parallelWorldName, place  );
        }
    } 

    //--- Build G4tgbVolume for parallel world
    std::vector<G4String> wl;
    G4tgrVolume* tgrVolPW = new G4tgrVolume( *tgrMassWorld );
    tgrVolPW->SetName( parallelWorldName );    
    tgrVolmgr->RegisterMe( tgrVolPW );
    G4tgbVolume* tgbVolPW = new G4tgbVolume( tgrVolPW );
    tgbVolmgr->RegisterMe( tgbVolPW );
    //-    BuildPhysicsProcess( tgrVoltop->GetName(), index );
  }

  return theParallelWorlds;
}

#include "G4ParallelWorldScoringProcess.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"

//--------------------------------------------------------------------
void G4tgbParallelGeomMgr::BuildPhysicsProcess( const G4String& volName, const G4int index )
{
  //  G4cout << " BuildPhysicsProcess "  << index << G4endl;
  G4ParallelWorldScoringProcess* theParallelWorldScoringProcess
    = new G4ParallelWorldScoringProcess("parallelWorldProcess_"+G4UIcommand::ConvertToString(index));
  theParallelWorldScoringProcess->SetParallelWorld(volName);
  
  
  G4ParticleTable::G4PTblDicIterator* theParticleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    pmanager->AddProcess(theParallelWorldScoringProcess);
    pmanager->SetProcessOrderingToLast(theParallelWorldScoringProcess, idxAtRest);
    pmanager->SetProcessOrdering(theParallelWorldScoringProcess, idxAlongStep, 1);
    pmanager->SetProcessOrderingToLast(theParallelWorldScoringProcess, idxPostStep);
  }  
}

