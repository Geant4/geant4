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

#include "G4WorkerRunManager.hh"
#include "G4WorkerRunManagerKernel.hh"
#include "G4UImanager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MTRunManager.hh"
#include "G4ScoringManager.hh"
#include "G4TransportationManager.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4WorkerThread.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4UserRunAction.hh"

G4WorkerRunManager::G4WorkerRunManager() : G4RunManager(true) {
    //This constructor should never be called in non-multithreaded mode
#ifndef G4MULTITHREADED
    G4ExceptionDescription msg;
    msg<<"Geant4 code is compiled without multi-threading support (-DG4MULTITHREADED is set to off).";
    msg<<" This type of RunManager can only be used in mult-threaded applications.";
    G4Exception("G4WorkerRunManager::G4WorkerRunManager()","Run0035",FatalException,msg);
#endif
    G4ParticleTable::GetParticleTable()->SlaveG4ParticleTable();
    //Add this worker to the list of workers
    G4MTRunManager::GetMasterRunManager()->AddWorkerRunManager(this);
    G4ScoringManager* masterScM = G4MTRunManager::GetMasterScoringManager();
    if(masterScM) G4ScoringManager::GetScoringManager(); //TLS instance for a worker
}

G4WorkerRunManager::~G4WorkerRunManager() {
    G4cout<<"Destroying WorkerRunManager"<<G4endl;
}

void G4WorkerRunManager::DeleteUserDetector()
{
    //Nothing to do: userDetecotr pointer is owned by master thread
    userDetector = 0;
}


void G4WorkerRunManager::InitializeGeometry() {
    if(!userDetector)
    {
        G4Exception("G4RunManager::InitializeGeometry", "Run0033",
                    FatalException, "G4VUserDetectorConstruction is not defined!");
        return;
    }
    //TODO: add new revised way of doing this after Boston, that is:
    //Step1: Call user's ConstructSDandField()
    userDetector->ConstructSDandField();
    userDetector->ConstructParallelSD();
    //Step2: Get pointer to the physiWorld (note: needs to get the "super pointer, i.e. the one shared by all threads"
    //G4RunManager* masterRM = G4MTRunManager::GetMasterRunManager();
    G4RunManagerKernel* masterKernel = G4MTRunManager::GetMasterRunManagerKernel();
    G4VPhysicalVolume* worldVol = masterKernel->GetCurrentWorld();
    //Step3:, Call a new "SlaveDefineWorldVolume( pointer from 2-, false); //TODO: Change name
    kernel->SlaveDefineWorldVolume(worldVol,false);
    
    nParallelWorlds = userDetector->ConstructParallelGeometries();
    kernel->SetNumberOfParallelWorld(nParallelWorlds);
    geometryInitialized = true;
}

void G4WorkerRunManager::DoEventLoop(G4int n_event, const char* macroFile , G4int n_select)
{
    //This is the same as in the sequential case, just the for-loop indexes are
    //different
    InitializeEventLoop(n_event,macroFile,n_select);
    //Signal this thread can start event loop.
    //Note this will return only when all threads reach this point
    G4MTRunManager::GetMasterRunManager()->ThisWorkerReady();
    // Event loop
    for ( G4int i_event = workerContext->GetThreadId(); i_event<n_event; i_event += workerContext->GetNumberThreads() )
    {
        ProcessOneEvent(i_event);
        TerminateOneEvent();
        if(runAborted) break;
    }
    //Signal this thread has finished envent-loop.
    //Note this will return only whan all threads reach this point
    G4MTRunManager::GetMasterRunManager()->ThisWorkerEndEventLoop();
    
    TerminateEventLoop();

    G4MTRunManager* mtRM = G4MTRunManager::GetMasterRunManager();
    G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
    if(ScM) mtRM->MergeScores(ScM);
    mtRM->MergeRun(currentRun);
}

void G4WorkerRunManager::ProcessOneEvent(G4int i_event)
{
    //Need to reseed random number generator
    long seeds[3] = { G4MTRunManager::GetSeed(i_event*2) , G4MTRunManager::GetSeed(i_event*2+1), 0 };
    G4Random::setTheSeeds(seeds);
    currentEvent = GenerateEvent(i_event);
    eventManager->ProcessOneEvent(currentEvent);
    AnalyzeEvent(currentEvent);
    UpdateScoring();
    if(i_event<n_select_msg) G4UImanager::GetUIpointer()->ApplyCommand(msgText);
}

void G4WorkerRunManager::ConstructScoringWorlds()
{
    // Do the correct stuff ...
    G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
    if(!ScM) return;
    G4int nPar = ScM->GetNumberOfMesh();
    if(nPar<1) return;
    
    G4ScoringManager* masterScM = G4MTRunManager::GetMasterScoringManager();
    assert( masterScM != NULL );
    
    G4ParticleTable::G4PTblDicIterator* particleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
    
    G4int counter = 0;
    for(G4int iw=0;iw<nPar;iw++)
    {
        G4VScoringMesh* mesh = ScM->GetMesh(iw);
        G4VScoringMesh* masterMesh = masterScM->GetMesh(iw);
/************************************************************************
        //Try copy the original mesh to the worker mesh. Only few mesh types are supported, only the correct one will work down here
        meshcopy<G4ScoringBox>(masterMesh,mesh);
        meshcopy<G4ScoringCylinder>(masterMesh,mesh );
*************************************************************************/
        mesh->SetMeshElementLogical(masterMesh->GetMeshElementLogical());
        
        G4VPhysicalVolume* pWorld = G4TransportationManager::GetTransportationManager()->IsWorldExisting(ScM->GetWorldName(iw));
        if(pWorld)
        {
            //Retrieve from master thread parallel world pointer
            G4MTRunManager::masterWorlds_t::const_iterator it = G4MTRunManager::GetMasterWorlds().find(counter++);
            assert( it != G4MTRunManager::GetMasterWorlds().end() );
            
            G4TransportationManager::GetTransportationManager()->RegisterWorld( (*it).second );
            
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
        }//if(!pWorld)
        mesh->WorkerConstruct(pWorld);
    }
}

void G4WorkerRunManager::SetUserInitialization(G4UserWorkerInitialization*)
{
    G4Exception("G4RunManager::SetUserInitialization(G4UserWorkerInitialization*)", "Run3021",
                FatalException, "This method should be used only with an instance of G4MTRunManager");
}

void G4WorkerRunManager::SetUserInitialization(G4VUserActionInitialization*)
{
    G4Exception("G4RunManager::SetUserInitialization(G4VUserActionInitialization*)", "Run3021",
                FatalException, "This method should be used only with an instance of G4MTRunManager");
}

void G4WorkerRunManager::SetUserInitialization(G4VUserDetectorConstruction*)
{
    G4Exception("G4RunManager::SetUserInitialization(G4VUserDetectorConstruction*)", "Run3021",
                FatalException, "This method should be used only with an instance of G4MTRunManager");
}

void G4WorkerRunManager::SetUserInitialization(G4VUserPhysicsList* pl)
{
    pl->InitializeWorker();
    G4RunManager::SetUserInitialization(pl);
}

void G4WorkerRunManager::SetUserAction(G4UserRunAction* userAction)
{
    userRunAction = userAction;
    userRunAction->SetMaster(false);
}

void G4WorkerRunManager::SetupDefaultRNGEngine()
{
    const CLHEP::HepRandomEngine* mrnge = G4MTRunManager::GetMasterRunManager()->getMasterRandomEngine();
    assert(mrnge);//Master has created RNG
    const G4UserWorkerInitialization* uwi =G4MTRunManager::GetMasterRunManager()->GetUserWorkerInitialization();
    uwi->SetupRNGEngine(mrnge);
}


//Forward calls (avoid GCC compilation warnings)
void G4WorkerRunManager::SetUserAction(G4UserEventAction* ua)
{
    G4RunManager::SetUserAction(ua);
}

void G4WorkerRunManager::SetUserAction(G4VUserPrimaryGeneratorAction* ua)
{
    G4RunManager::SetUserAction(ua);
}

void G4WorkerRunManager::SetUserAction(G4UserStackingAction* ua)
{
    G4RunManager::SetUserAction(ua);
}

void G4WorkerRunManager::SetUserAction(G4UserTrackingAction* ua)
{
    G4RunManager::SetUserAction(ua);
}

void G4WorkerRunManager::SetUserAction(G4UserSteppingAction* ua)
{
    G4RunManager::SetUserAction(ua);
}
