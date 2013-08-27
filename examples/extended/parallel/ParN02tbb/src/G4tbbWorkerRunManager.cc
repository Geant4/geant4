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
#include "G4tbbWorkerRunManager.hh"

#include "G4tbbRunManager.hh"
#include "G4RunManagerKernel.hh"

#include "G4VUserDetectorConstruction.hh"

#include "G4AutoLock.hh"
// #include "G4Random.hh" 

G4Mutex countWorkerMtx = G4MUTEX_INITIALIZER;

G4tbbWorkerRunManager::G4tbbWorkerRunManagerInstancesType 
 G4tbbWorkerRunManager::instancesList=
   G4tbbWorkerRunManager::G4tbbWorkerRunManagerInstancesType();

unsigned int G4tbbWorkerRunManager::fWorkerCounter=0;  // Do not use id=0 
                                                    // -- as though master used it
G4bool G4tbbWorkerRunManager::fUseCounterId= true;    // By default use a count for Id

G4tbbWorkerRunManager::G4tbbWorkerRunManager() 
   : G4RunManager(true),  // Ensure use of Run Manager for workers
     fWorkerContext(0), 
     fMasterRM(0)
{
  // Register this RM - needed to eventually free/destroy resources
  instancesList.push(this);

  fMasterRM= G4tbbRunManager::GetMasterTbbRunManager(); 
  fWorkerContext->BuildGeometryAndPhysicsVector(); 

  //  InitSlave();
  fWorkerContext= new G4WorkerThread(); 

  // Lock only to increment count
  { 
     G4AutoLock al( &countWorkerMtx ); 
     ++fWorkerCounter;  
  }
  if( fUseCounterId ) {

     fWorkerContext->SetThreadId(fWorkerCounter);  
  }
  // Could define identify from order of creation - or 
  //  potentially from number of last event executed.
}

G4tbbWorkerRunManager::~G4tbbWorkerRunManager() 
{
  fWorkerContext->DestroyGeometryAndPhysicsVector(); 

  //  These would do something like:  
  //    tbbSlaveDestroyGeometryAndPhysicsVector(); 

  delete fWorkerContext; 
  fWorkerContext= 0;
}

void G4tbbWorkerRunManager::ProcessOneEvent( G4int i_event)
{
   // In MT the seeds are set first

   // Do the rest
   currentEvent = GenerateEvent(i_event);
   eventManager->ProcessOneEvent(currentEvent);
   AnalyzeEvent(currentEvent);
   UpdateScoring();

   // if(i_event<n_select) G4UImanager::GetUIpointer()->ApplyCommand(msgText);
}

G4bool G4tbbWorkerRunManager::DoOneEvent( G4int i_event, G4String msg)
{
   // This method is called by the tbb::task. 
   // It 
   //   i) performs the simulation of a single event
   //  ii) obtains the status
   // iii) stacks the event (!?)

   // Set the UI command text 
   msgText= msg;

   // Do the basic work
   ProcessOneEvent( i_event ); 

   StackPreviousEvent(currentEvent);
   currentEvent = 0;
   //G4cout<<"runAborted is:"<<runAborted<<G4endl;
   return runAborted;
}

void G4tbbWorkerRunManager::DoEventLoop(G4int ,      // n_event, 
                                        const char*, // macroFile , 
                                        G4int )      // n_select)
{
   G4cerr << " ERROR> The event Loop is not handled by the (worker) "
          << " Run Manager - it is handled by TBB" << G4endl;
}



G4tbbWorkerRunManager::G4tbbWorkerRunManagerInstancesType& 
  G4tbbWorkerRunManager::GetInstancesList()
{
  return instancesList; 
}

unsigned int G4tbbWorkerRunManager::NumberOfWorkers()
{
   return instancesList.unsafe_size();
}

void G4tbbWorkerRunManager::DestroyWorkersAndCleanup()
{
   // instancesList.ClearAndDestroy(); 
}

void G4tbbWorkerRunManager::InitializeGeometry() 
{
   // Copied from G4WorkerRunManager::InitializeGeometry() - June 27, 2013 JA

    if(!userDetector)
    {
        G4Exception("G4tbbWorkerRunManager::InitializeGeometry", "Run0034",
                    FatalException, "G4VUserDetectorConstruction is not defined!");
        return;
    }
    //Step1: Call user's ConstructSDandField()
    userDetector->ConstructSDandField();
    userDetector->ConstructParallelSD();

    // TODO: extend to scoring .. JA

    //Step2: Get pointer to the physiWorld (note: needs to get 
    //   the "super pointer, i.e. the one shared by all threads"

    // G4VPhysicalVolume* worldVol = (fMasterRM->kernel)->GetCurrentWorld();
    //Step3:, Call a new "WorkerDefineWorldVolume( pointer from 2-, false); 
    // kernel->WorkerDefineWorldVolume(worldVol,false);
    
    // kernel->SetNumberOfParallelWorld(fMasterRM->Kernel->GetNumberOfParallelWorld());
    geometryInitialized = true;
}

void G4tbbWorkerRunManager::RunTermination()
{
    //Merge partial results into global run
    // G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
    // if(ScM) fMasterRM->MergeScores(ScM);

    // Merge scoring from different workers
    // fMasterRM->MergeRun(currentRun);

    G4RunManager::RunTermination();
    //Signal this thread has finished envent-loop.
    //Note this will return only whan all threads reach this point
}

void G4tbbWorkerRunManager::ConstructScoringWorlds()
{
    G4Exception("G4tbbWorkerRunManager::ConstructScoringWorlds", "Run0035",
               FatalException, "ConstructScoringWorlds is not defined!");
   /****
    // Do the correct stuff ...
    G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
    if(!ScM) return;
    G4int nPar = ScM->GetNumberOfMesh();
    if(nPar<1) return;
    
    G4ScoringManager* masterScM = G4MTRunManager::GetMasterScoringManager();
    assert( masterScM != NULL );
    
    G4ParticleTable::G4PTblDicIterator* particleIterator 
       = G4ParticleTable::GetParticleTable()->GetIterator();
    
    for(G4int iw=0;iw<nPar;iw++)
    {
        G4VScoringMesh* mesh = ScM->GetMesh(iw);
        G4VScoringMesh* masterMesh = masterScM->GetMesh(iw);
        mesh->SetMeshElementLogical(masterMesh->GetMeshElementLogical());
        
        G4VPhysicalVolume* pWorld
          = G4TransportationManager::GetTransportationManager()
            ->GetParallelWorld(ScM->GetWorldName(iw));
        G4ParallelWorldProcess* theParallelWorldProcess
          = new G4ParallelWorldProcess(ScM->GetWorldName(iw));
        theParallelWorldProcess->SetParallelWorld(ScM->GetWorldName(iw));

        particleIterator->reset();
        while( (*particleIterator)() ){
          G4ParticleDefinition* particle = particleIterator->value();
          G4ProcessManager* pmanager = particle->GetProcessManager();
          if(pmanager)
          {
            pmanager->AddProcess(theParallelWorldProcess);
            if(theParallelWorldProcess->IsAtRestRequired(particle))
            { pmanager->SetProcessOrdering(theParallelWorldProcess, idxAtRest, 9999); }
            pmanager->SetProcessOrderingToSecond(theParallelWorldProcess, idxAlongStep);
            pmanager->SetProcessOrdering(theParallelWorldProcess, idxPostStep, 9999);
          } //if(pmanager)
        }//while
        mesh->WorkerConstruct(pWorld);
    }
   ***/ 
}

void G4tbbWorkerRunManager::StoreRNGStatus(const G4String& fn )
{
    std::ostringstream os;
    os << randomNumberStatusDir << "G4Worker"<< fWorkerContext->GetThreadId()<<"_"<<fn <<".rndm";
    G4Random::saveEngineStatus(os.str().c_str());    
}
