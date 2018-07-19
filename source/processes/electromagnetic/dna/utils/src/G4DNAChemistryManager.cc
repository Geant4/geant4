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
// $Id: G4DNAChemistryManager.cc 103042 2017-03-10 11:50:07Z gcosmo $
//
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4DNAChemistryManager.hh"

#include "G4Scheduler.hh"
#include "G4SystemOfUnits.hh"
#include "G4Molecule.hh"
#include "G4VITTrackHolder.hh"
#include "G4H2O.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4DNAWaterIonisationStructure.hh"
#include "G4Electron_aq.hh"
#include "G4MolecularConfiguration.hh"
#include "G4VMoleculeCounter.hh"
#include "G4VUserChemistryList.hh"
#include "G4AutoLock.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4MoleculeFinder.hh"
#include "G4MoleculeTable.hh"
#include "G4PhysChemIO.hh"

using namespace std;

G4DNAChemistryManager* G4DNAChemistryManager::fgInstance;

G4ThreadLocal G4DNAChemistryManager::ThreadLocalData*
  G4DNAChemistryManager::fpThreadData = 0;

G4Mutex chemManExistence;

//------------------------------------------------------------------------------

G4DNAChemistryManager::ThreadLocalData::ThreadLocalData()
{
  fpPhysChemIO = nullptr;
  fThreadInitialized_tl = false;
}

//------------------------------------------------------------------------------

G4DNAChemistryManager::ThreadLocalData::~ThreadLocalData()
{
  if(fpPhysChemIO) delete fpPhysChemIO;
  fThreadInitialized_tl = false;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::SetPhysChemIO(G4VPhysChemIO* physChemIO)
{
  if(fpThreadData->fpPhysChemIO)
    delete fpThreadData->fpPhysChemIO;
  fpThreadData->fpPhysChemIO = physChemIO;
}

//------------------------------------------------------------------------------

G4DNAChemistryManager::G4DNAChemistryManager() :
    G4UImessenger(), G4VStateDependent()
{

//------------------------------------------------------------------------------
/* 
 * The chemistry manager is shared between threads
 * It is initialized both on the master thread and on the worker threads
 */
//------------------------------------------------------------------------------
  
  fpExcitationLevel = 0;
  fpIonisationLevel = 0;
  fpUserChemistryList = 0;
  fMasterInitialized = false;
  fpChemDNADirectory = new G4UIdirectory("/chem/");
  fpActivateChem = new G4UIcmdWithABool("/chem/activate", this);
  fpRunChem = new G4UIcmdWithoutParameter("/chem/run", this);
  //fpGridSize = new G4UIcmdWithADoubleAndUnit("/chem/gridRes", this);
  fpScaleForNewTemperature = new G4UIcmdWithADoubleAndUnit("/chem/temperature",
                                                           this);
  fpSkipReactionsFromChemList =
      new G4UIcmdWithoutParameter("/chem/skipReactionsFromChemList", this);
  fpInitChem = new G4UIcmdWithoutParameter("/chem/init", this);
  //fDefaultGridResolution = -1;
  fBuildPhysicsTable = false;
  fGeometryClosed = false;
  fPhysicsTableBuilt = false;
  fForceThreadReinitialization = false;
  fVerbose = 0;
  fActiveChemistry = false;
  fSkipReactions = false;
  fResetCounterWhenRunEnds = true;
}

//------------------------------------------------------------------------------

G4DNAChemistryManager*
G4DNAChemistryManager::Instance()
{
  if(fgInstance == 0)
  {
    G4AutoLock lock(&chemManExistence);
    if (fgInstance == 0) // MT : double check at initialisation
    {
      fgInstance = new G4DNAChemistryManager();
    }
    lock.unlock();
  }
  
  if(fpThreadData==0) fpThreadData = new ThreadLocalData();
   // make sure thread local data is initialized for all threads
  
  return fgInstance;
}

//------------------------------------------------------------------------------

G4DNAChemistryManager*
G4DNAChemistryManager::GetInstanceIfExists()
{
  return fgInstance;
}

//------------------------------------------------------------------------------

G4DNAChemistryManager::~G4DNAChemistryManager()
{
//	G4cout << "Deleting G4DNAChemistryManager" << G4endl;
  Clear();
  fgInstance = 0;
  /*
   * DEBUG : check that the chemistry manager has well been deregistered
   *  assert(G4StateManager::GetStateManager()->
   *  DeregisterDependent(this) == true);
   */
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::Clear()
{
  if (fpIonisationLevel)
  {
    delete fpIonisationLevel;
    fpIonisationLevel = 0;

  }
  if (fpExcitationLevel)
  {
    delete fpExcitationLevel;
    fpExcitationLevel = 0;
  }
  if (fpUserChemistryList)
  {
    if(fpUserChemistryList->IsPhysicsConstructor() == false)
    {
      delete fpUserChemistryList;
    }
//    else
//    {
//   G4cout << "G4DNAChemistryManager will not delete the chemistry list "
//       "since it inherits from G4VPhysicsConstructor and it is then "
//       "expected to be the responsability to the G4VModularPhysics to handle"
//       " the chemistry list." << G4endl;
//    }
    fpUserChemistryList = 0;
  }
  if (fpChemDNADirectory)
  {
    delete fpChemDNADirectory;
    fpChemDNADirectory = 0;
  }
  if (fpActivateChem)
  {
    delete fpActivateChem;
    fpActivateChem = 0;
  }
  if(fpRunChem)
  {
    delete fpRunChem;
    fpRunChem = 0;
  }
  if(fpSkipReactionsFromChemList)
  {
    delete fpSkipReactionsFromChemList;
    fpSkipReactionsFromChemList = 0;
  }
  if(fpInitChem)
  {
    delete fpInitChem;
    fpInitChem = 0;
  }

  G4DNAMolecularReactionTable::DeleteInstance();
  //G4MoleculeHandleManager::DeleteInstance();
  G4MolecularConfiguration::DeleteManager();
  G4VMoleculeCounter::DeleteInstance();
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::DeleteInstance()
{
  //G4cout << "G4DNAChemistryManager::DeleteInstance" << G4endl;

  G4AutoLock lock(&chemManExistence);

  if(fgInstance)
  {
    G4DNAChemistryManager* deleteMe = fgInstance;
    fgInstance = 0;
    lock.unlock();
    delete deleteMe;
  }
  else
  {
    G4cout << "G4DNAChemistryManager already deleted" << G4endl;
  }
  lock.unlock();
}

//------------------------------------------------------------------------------

G4bool G4DNAChemistryManager::Notify(G4ApplicationState requestedState)
{
  if (requestedState == G4State_Quit)
  {
    if(fVerbose)
    G4cout << "G4DNAChemistryManager::Notify ---> received G4State_Quit"
           << G4endl;
    //DeleteInstance();
    Clear();
  }
//  else if(requestedState == G4State_EventProc)
//    Note: From here we can know that a new event is started
//          But the run and event IDs remain unknown
//  {
//    if(fpThreadData->fpPhysChemIO)
//    {
//      
//    }
//  }
  else if(requestedState == G4State_GeomClosed)
  {
    fGeometryClosed = true;
  }

  else if (requestedState == G4State_Idle)
  {
    G4MoleculeTable::Instance()->PrepareMolecularConfiguration();
  }

  return true;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::SetNewValue(G4UIcommand* command, G4String value)
{
  if (command == fpActivateChem)
  {
    Activated(G4UIcmdWithABool::GetNewBoolValue(value));
  }
  else if (command == fpRunChem)
  {
    Run();
  }
  /*
  else if(command == fpGridSize)
  {
    fDefaultGridResolution = fpGridSize->ConvertToDimensionedDouble(value);
  }*/
  else if (command == fpSkipReactionsFromChemList)
  {
    fSkipReactions = true;
  }
  else if(command == fpScaleForNewTemperature)
  {
    SetGlobalTemperature(fpScaleForNewTemperature->
                         ConvertToDimensionedDouble(value));
  }
  else if(command == fpInitChem)
  {
    Initialize();
    InitializeThread();
  }
}

//------------------------------------------------------------------------------

G4String G4DNAChemistryManager::GetCurrentValue(G4UIcommand* command)
{
  if (command == fpActivateChem)
  {
    return G4UIcmdWithABool::ConvertToString(fActiveChemistry);
  }

  return "";
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::Run()
{
  if (fActiveChemistry)
  {
    InitializeThread();

    if (fMasterInitialized == false)
    {
      G4ExceptionDescription description;
      description << "Global components were not initialized.";
      G4Exception("G4DNAChemistryManager::Run", "MASTER_INIT", FatalException,
                  description);
    }

    if (fpThreadData->fThreadInitialized_tl == 0)
    {
      G4ExceptionDescription description;
      description << "Thread local components were not initialized.";
      G4Exception("G4DNAChemistryManager::Run", "THREAD_INIT", FatalException,
                  description);
    }
    
    G4MoleculeTable::Instance()->Finalize();
    G4Scheduler::Instance()->Process();
    if(fResetCounterWhenRunEnds)
    {
      G4VMoleculeCounter::Instance()->ResetCounter();
    }
    CloseFile();
  }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::Gun(G4ITGun* gun, bool physicsTableToBuild)
{
  fBuildPhysicsTable = physicsTableToBuild;
  G4Scheduler::Instance()->SetGun(gun);
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::Initialize()
{
  //===========================================================================
  // MT MODE
  //===========================================================================
  if(G4Threading::IsMultithreadedApplication())
  {
    //==========================================================================
    // ON WORKER THREAD
    //==========================================================================
    if(G4Threading::IsWorkerThread())
    {
      InitializeThread(); // Will create and initialize G4Scheduler
      return;
    }
    //==========================================================================
    // ON MASTER THREAD
    //==========================================================================
    else
    {
      InitializeMaster();
      return;
    }
  }
  //===========================================================================
  // IS NOT IN MT MODE
  //===========================================================================
  else
  {
    InitializeMaster();
    // In this case: InitializeThread is called when Run() is called
    return;
  }

}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::InitializeMaster()
{
  if (fMasterInitialized == false)
  {
    if(fVerbose)
    {
     G4cout << "G4DNAChemistryManager::InitializeMaster() is called" << G4endl;
    }

    G4Scheduler::Instance(); 
    // creates a concrete object of the scheduler 

    if (fpUserChemistryList)
    {
      fpUserChemistryList->ConstructDissociationChannels();
      if(fSkipReactions == false)
      {
        fpUserChemistryList->ConstructReactionTable(
            G4DNAMolecularReactionTable::GetReactionTable());
      }
      else
      {
        G4DNAMolecularReactionTable::GetReactionTable(); // init pointer
      }
      fMasterInitialized = true;
    }
    else
    {
      if (fActiveChemistry)
      {
        G4ExceptionDescription description;
        description << "No user chemistry list has been provided.";
        G4Exception("G4DNAChemistryManager::InitializeMaster", "NO_CHEM_LIST",
                    FatalException, description);
      }
    }
  }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::InitializeThread()
{
  if (fpThreadData->fThreadInitialized_tl == false
      || fForceThreadReinitialization == true)
  {
    if (fpUserChemistryList)
    {
      if(fVerbose)
      {
        G4cout << "G4DNAChemistryManager::InitializeThread() is called"
               << G4endl;
      }

      if (fBuildPhysicsTable && fPhysicsTableBuilt == false)
      {
        if(fVerbose)
        {
          G4cout << "G4DNAChemistryManager: Build the physics tables for "
              "molecules."
                 << G4endl;
        }

        fpUserChemistryList->BuildPhysicsTable();
        if (fGeometryClosed == false)
        {
          if(fVerbose)
          {
            G4cout << "G4DNAChemistryManager: Close geometry"
                   << G4endl;
          }

          G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
          // G4cout << "Start closing geometry." << G4endl;
          geomManager->OpenGeometry();
          geomManager->CloseGeometry(true, true);
          fGeometryClosed = true;
        }

        fPhysicsTableBuilt = true;
      }
      fpUserChemistryList->ConstructTimeStepModel(
          G4DNAMolecularReactionTable::GetReactionTable());
      G4Scheduler::Instance()->Initialize();

      fpThreadData->fThreadInitialized_tl = true;
    }
    else
    {
      G4ExceptionDescription description;
      description << "No user chemistry list has been provided.";
      G4Exception("G4DNAChemistryManager::InitializeThread", "NO_CHEM_LIST",
                  FatalException, description);
    }

    G4VMoleculeCounter::InitializeInstance();
  }

  InitializeFile();
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::InitializeFile()
{
  if(fVerbose)
  {
    G4cout << "G4DNAChemistryManager::InitializeFile() is called"
           << G4endl;
  }
  
  if(fpThreadData->fpPhysChemIO){
    fpThreadData->fpPhysChemIO->InitializeFile();
  }
}

//------------------------------------------------------------------------------

G4bool G4DNAChemistryManager::IsActivated()
{
  return Instance()->fActiveChemistry;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::Activated(G4bool flag)
{
  Instance()->fActiveChemistry = flag;
}

//------------------------------------------------------------------------------

G4bool G4DNAChemistryManager::IsChemistryActivated()
{
  return fActiveChemistry;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::SetChemistryActivation(G4bool flag)
{
  fActiveChemistry = flag;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::WriteInto(const G4String& output,
                                      ios_base::openmode mode)
{
  if (fVerbose)
  {
    G4cout << "G4DNAChemistryManager: Write chemical stage into "
           << output.data() << G4endl;
  }

  if(fpThreadData->fpPhysChemIO){
    fpThreadData->fpPhysChemIO->WriteInto(output, mode);
  }
  else{
    fpThreadData->fpPhysChemIO = new G4PhysChemIO::FormattedText();
    fpThreadData->fpPhysChemIO->WriteInto(output, mode);
  }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::AddEmptyLineInOuputFile()
{
  if(fpThreadData->fpPhysChemIO){
    fpThreadData->fpPhysChemIO->AddEmptyLineInOuputFile();
  }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::CloseFile()
{
  if(fpThreadData->fpPhysChemIO){
    fpThreadData->fpPhysChemIO->CloseFile();
  }
}

//------------------------------------------------------------------------------

G4DNAWaterExcitationStructure*
G4DNAChemistryManager::GetExcitationLevel()
{
  if (!fpExcitationLevel)
  {
    fpExcitationLevel = new G4DNAWaterExcitationStructure;
  }
  return fpExcitationLevel;
}

//------------------------------------------------------------------------------

G4DNAWaterIonisationStructure*
G4DNAChemistryManager::GetIonisationLevel()
{
  if (!fpIonisationLevel)
  {
    fpIonisationLevel = new G4DNAWaterIonisationStructure;
  }
  return fpIonisationLevel;
}

//------------------------------------------------------------------------------

void
G4DNAChemistryManager::CreateWaterMolecule(ElectronicModification modification,
                                           G4int electronicLevel,
                                           const G4Track* theIncomingTrack)
{
  if(fpThreadData->fpPhysChemIO){
    G4double energy = -1.;

    switch (modification)
    {
      case eDissociativeAttachment:
        energy = 0;
        break;
      case eExcitedMolecule:
        energy = GetExcitationLevel()->ExcitationEnergy(electronicLevel);
        break;
      case eIonizedMolecule:
        energy = GetIonisationLevel()->IonisationEnergy(electronicLevel);
        break;
    }
  
    fpThreadData->fpPhysChemIO->CreateWaterMolecule(modification,
                                      4-electronicLevel,
                                      energy,
                                      theIncomingTrack);
  }

  if(fActiveChemistry)
  {
    G4Molecule * H2O = new G4Molecule (G4H2O::Definition());

    switch (modification)
    {
      case eDissociativeAttachment:
      H2O -> AddElectron(5,1);
      break;
      case eExcitedMolecule :
      H2O -> ExciteMolecule(4-electronicLevel);
      break;
      case eIonizedMolecule :
      H2O -> IonizeMolecule(4-electronicLevel);
      break;
    }

    G4Track * H2OTrack = H2O->BuildTrack(1*picosecond,
        theIncomingTrack->GetPosition());

    H2OTrack -> SetParentID(theIncomingTrack->GetTrackID());
    H2OTrack -> SetTrackStatus(fStopButAlive);
    H2OTrack -> SetKineticEnergy(0.);
    G4VITTrackHolder::Instance()->Push(H2OTrack);
  }
}

//------------------------------------------------------------------------------

void
G4DNAChemistryManager::CreateSolvatedElectron(const G4Track* theIncomingTrack,
                                              G4ThreeVector* finalPosition)
// finalPosition is a pointer because this argument is optional
{
  if(fpThreadData->fpPhysChemIO){
    fpThreadData->fpPhysChemIO->CreateSolvatedElectron(theIncomingTrack,
                                         finalPosition);
  }

  if(fActiveChemistry)
  {
    G4Molecule* e_aq = new G4Molecule(G4Electron_aq::Definition());
    G4Track * e_aqTrack(0);
    if(finalPosition)
    {
      e_aqTrack = e_aq->BuildTrack(picosecond,*finalPosition);
    }
    else
    {
      e_aqTrack = e_aq->BuildTrack(picosecond,theIncomingTrack->GetPosition());
    }
    e_aqTrack -> SetTrackStatus(fAlive);
    e_aqTrack -> SetParentID(theIncomingTrack->GetTrackID());
    G4VITTrackHolder::Instance()->Push(e_aqTrack);
  }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::PushMolecule(G4Molecule*& molecule,
                                         double time,
                                         const G4ThreeVector& position,
                                         int parentID)
{
  // TODO: PhysChemIO - method unused in the released code
  
  if(fActiveChemistry)
  {
    G4Track* track = molecule->BuildTrack(time,position);
    track -> SetTrackStatus(fAlive);
    track -> SetParentID(parentID);
    G4VITTrackHolder::Instance()->Push(track);
  }
  else
  {
    delete molecule;
    molecule = 0;
  }
}

//------------------------------------------------------------------------------


void G4DNAChemistryManager::
PushMoleculeAtParentTimeAndPlace(G4Molecule*& molecule,
                                 const G4Track* theIncomingTrack)
{
  // TODO: PhysChemIO - method unused in the released code

  if(fActiveChemistry)
  {
    G4Track* track = molecule->BuildTrack(theIncomingTrack->GetGlobalTime(),
                                          theIncomingTrack->GetPosition());
    track -> SetTrackStatus(fAlive);
    track -> SetParentID(theIncomingTrack->GetTrackID());
    G4VITTrackHolder::Instance()->Push(track);
  }
  else
  {
    delete molecule;
    molecule = 0;
  }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::SetGlobalTemperature(G4double temp_K)
{
  G4MolecularConfiguration::SetGlobalTemperature(temp_K);
  G4DNAMolecularReactionTable::Instance()->
    ScaleReactionRateForNewTemperature(temp_K);
}

