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
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disappear in the next releases.
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
#include "G4UIcmdWithAnInteger.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4MoleculeFinder.hh"
#include "G4MoleculeTable.hh"
#include "G4PhysChemIO.hh"


G4DNAChemistryManager* G4DNAChemistryManager::fgInstance = nullptr;

G4ThreadLocal G4DNAChemistryManager::ThreadLocalData*
    G4DNAChemistryManager::fpThreadData = nullptr;

G4Mutex chemManExistence;

//------------------------------------------------------------------------------

G4DNAChemistryManager::ThreadLocalData::ThreadLocalData()
{
    fpPhysChemIO = nullptr;
    fThreadInitialized = false;
}

//------------------------------------------------------------------------------

G4DNAChemistryManager::ThreadLocalData::~ThreadLocalData()
{
    fpThreadData = nullptr;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::SetPhysChemIO(std::unique_ptr<G4VPhysChemIO> pPhysChemIO)
{
    fpThreadData->fpPhysChemIO = std::move(pPhysChemIO);
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/*
 * The chemistry manager is shared between threads
 * It is initialized both on the master thread and on the worker threads
 */
//------------------------------------------------------------------------------
G4DNAChemistryManager::G4DNAChemistryManager()
        : G4UImessenger()
        , G4VStateDependent()
        , fpChemDNADirectory(new G4UIdirectory("/chem/"))
        , fpActivateChem(new G4UIcmdWithABool("/chem/activate", this))
        , fpRunChem(new G4UIcmdWithAnInteger("/chem/run", this))
        , fpSkipReactionsFromChemList(new G4UIcmdWithoutParameter("/chem/skipReactionsFromChemList", this))
        , fpScaleForNewTemperature(new G4UIcmdWithADoubleAndUnit("/chem/temperature", this))
        , fpInitChem(new G4UIcmdWithoutParameter("/chem/init", this))
        , fActiveChemistry(false)
        , fMasterInitialized(false)
        , fForceThreadReinitialization(false)
        , fpExcitationLevel(nullptr)
        , fpIonisationLevel(nullptr)
        , fpUserChemistryList(nullptr)
        , fOwnChemistryList(false)
        , fUseInStandalone(false)
        , fPhysicsTableBuilt(false)
        , fSkipReactions(false)
        , fGeometryClosed(false)
        , fVerbose(0)
        , fResetCounterWhenRunEnds(true)
{
    fpRunChem->SetParameterName("Number of runs to execute for the chemistry module"
                                "(this works when used in standalone", true, true);
    fpRunChem->SetDefaultValue(1);
    fpScaleForNewTemperature->SetUnitCategory("Temperature");
}

//------------------------------------------------------------------------------

G4DNAChemistryManager* G4DNAChemistryManager::Instance()
{
    if (fgInstance == nullptr)
    {
        G4AutoLock lock(&chemManExistence);
        if (fgInstance == nullptr) // MT : double check at initialisation
        {
            fgInstance = new G4DNAChemistryManager();
        }
        lock.unlock();
    }

    // make sure thread local data is initialized for all threads
    if (fpThreadData == nullptr)
    {
        fpThreadData = new ThreadLocalData();
    }

    assert(fpThreadData != nullptr);

    return fgInstance;
}

//------------------------------------------------------------------------------

G4DNAChemistryManager* G4DNAChemistryManager::GetInstanceIfExists()
{
    return fgInstance;
}

//------------------------------------------------------------------------------

G4DNAChemistryManager::~G4DNAChemistryManager()
{
    Clear();
    fgInstance = nullptr;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::Clear()
{
    fpIonisationLevel.reset();
    fpExcitationLevel.reset();

    if (fpUserChemistryList)
    {
        Deregister(*fpUserChemistryList);
    }

    fpChemDNADirectory.reset();
    fpActivateChem.reset();
    fpRunChem.reset();

    fpSkipReactionsFromChemList.reset();
    fpInitChem.reset();

    if (fpThreadData != nullptr)
    {
        delete fpThreadData;
        fpThreadData = nullptr;
    }

    G4DNAMolecularReactionTable::DeleteInstance();
    G4MolecularConfiguration::DeleteManager();
    G4VMoleculeCounter::DeleteInstance();
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::DeleteInstance()
{
    G4AutoLock lock(&chemManExistence);

    if (fgInstance != nullptr)
    {
        G4DNAChemistryManager* pDeleteMe = fgInstance;
        fgInstance = nullptr;
        lock.unlock();
        delete pDeleteMe;
    }
    else
    {
        G4cerr << "G4DNAChemistryManager already deleted" << G4endl;
    }
    lock.unlock();
}

//------------------------------------------------------------------------------

G4bool G4DNAChemistryManager::Notify(G4ApplicationState requestedState)
{
    if (requestedState == G4State_Quit)
    {
        if (fVerbose)
        {
            G4cout << "G4DNAChemistryManager::Notify ---> received G4State_Quit"
                   << G4endl;
        }
        Clear();
    }
    else if (requestedState == G4State_GeomClosed)
    {
        fGeometryClosed = true;
    }
    else if (requestedState == G4State_Idle)
    {
        InitializeThreadSharedData();
    }

    return true;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::SetNewValue(G4UIcommand* pCommand, G4String value)
{
    if (pCommand == fpActivateChem.get())
    {
        SetChemistryActivation(G4UIcmdWithABool::GetNewBoolValue(value));
    }
    else if (pCommand == fpRunChem.get())
    {
        int nbExec = value.empty() ? 1 : G4UIcommand::ConvertToInt(value);
        for (int i = 0 ; i < nbExec ; ++i)
        {
            Run();
        }
    }
    else if (pCommand == fpSkipReactionsFromChemList.get())
    {
        fSkipReactions = true;
    }
    else if (pCommand == fpScaleForNewTemperature.get())
    {
        SetGlobalTemperature(fpScaleForNewTemperature->ConvertToDimensionedDouble(value));
    }
    else if (pCommand == fpInitChem.get())
    {
        Initialize();
        InitializeThread();
    }
}

//------------------------------------------------------------------------------

G4String G4DNAChemistryManager::GetCurrentValue(G4UIcommand* pCommand)
{
    if (pCommand == fpActivateChem.get())
    {
        return G4UIcmdWithABool::ConvertToString(fActiveChemistry);
    }
    else if (pCommand == fpScaleForNewTemperature.get())
    {
        return fpScaleForNewTemperature->ConvertToStringWithBestUnit(G4MolecularConfiguration::GetGlobalTemperature());
    }
    else if (pCommand == fpSkipReactionsFromChemList.get())
    {
        return G4UIcmdWithABool::ConvertToString(fSkipReactions);
    }

    return "";
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::InitializeThreadSharedData()
{
    if (!G4Threading::IsMasterThread())
    {
        return;
    }

    G4MoleculeTable::Instance()->PrepareMolecularConfiguration();
    G4MoleculeTable::Instance()->Finalize();
}

//------------------------------------------------------------------------------
void G4DNAChemistryManager::Run()
{
    if (!fActiveChemistry)
    {
        return;
    }

    InitializeThread();

    if (!fMasterInitialized)
    {
        G4ExceptionDescription description;
        description << "Global components were not initialized.";
        G4Exception("G4DNAChemistryManager::Run", "MASTER_INIT", FatalException,
                    description);
    }

    if (!fpThreadData->fThreadInitialized)
    {
        G4ExceptionDescription description;
        description << "Thread local components were not initialized.";
        G4Exception("G4DNAChemistryManager::Run", "THREAD_INIT", FatalException,
                    description);
    }
    
    G4MoleculeTable::Instance()->Finalize();
    G4Scheduler::Instance()->Process();
    if (fResetCounterWhenRunEnds)
    {
        G4VMoleculeCounter::Instance()->ResetCounter();
    }
    CloseFile();
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::UseAsStandalone(G4bool flag)
{
    fUseInStandalone = flag;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::SetGun(G4ITGun* pChemGun)
{
    G4Scheduler::Instance()->SetGun(pChemGun);
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::Initialize()
{
    //===========================================================================
    // MT MODE
    //===========================================================================
    if (G4Threading::IsMultithreadedApplication())
    {
        //==========================================================================
        // ON WORKER THREAD
        //==========================================================================
        if (G4Threading::IsWorkerThread())
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
            InitializeThreadSharedData();
            return;
        }
    }
    //===========================================================================
    // IS NOT IN MT MODE
    //===========================================================================
    else
    {
        InitializeMaster();
        InitializeThreadSharedData();
        // In this case: InitializeThread is called when Run() is called
        return;
    }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::InitializeMaster()
{
    if (fMasterInitialized)
    {
        return;
    }

    if (fVerbose)
    {
        G4cout << "G4DNAChemistryManager::InitializeMaster() is called" << G4endl;
    }


    if (fpUserChemistryList == nullptr)
    {
        G4ExceptionDescription description;
        description << "No user chemistry list has been provided.";
        G4Exception("G4DNAChemistryManager::InitializeMaster", "NO_CHEM_LIST",
                    FatalException, description);
    }else
    {
      fpUserChemistryList->ConstructDissociationChannels();
      if (!fSkipReactions)
      {
          fpUserChemistryList->ConstructReactionTable(G4DNAMolecularReactionTable::GetReactionTable());
      }
      else
      {
          G4DNAMolecularReactionTable::GetReactionTable(); // init pointer
      }
    }

    G4Scheduler::Instance();
    // creates a concrete object of the scheduler
    fMasterInitialized = true;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::HandleStandaloneInitialization()
{
    if (!fUseInStandalone || fPhysicsTableBuilt)
    {
        return;
    }

    if (fVerbose)
    {
        G4cout << "G4DNAChemistryManager: Build the physics tables for "
                  "molecule definition only."
               << G4endl;
    }

    fpUserChemistryList->BuildPhysicsTable();

    if (!fGeometryClosed)
    {
        if (fVerbose)
        {
            G4cout << "G4DNAChemistryManager: Close geometry" << G4endl;
        }

        G4GeometryManager* pGeomManager = G4GeometryManager::GetInstance();
        pGeomManager->OpenGeometry();
        pGeomManager->CloseGeometry(true, true);
        fGeometryClosed = true;
    }

    fPhysicsTableBuilt = true;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::InitializeThread()
{
    if (fpThreadData->fThreadInitialized && !fForceThreadReinitialization)
    {
        return;
    }

    if (fpUserChemistryList == nullptr)
    {
        G4ExceptionDescription description;
        description << "No user chemistry list has been provided.";
        G4Exception("G4DNAChemistryManager::InitializeThread", "NO_CHEM_LIST",
                    FatalException, description);
    }else
    {
        HandleStandaloneInitialization();// To make coverty happy
        fpUserChemistryList->ConstructTimeStepModel(G4DNAMolecularReactionTable::GetReactionTable());
    }

    if (fVerbose)
    {
        G4cout << "G4DNAChemistryManager::InitializeThread() is called"
               << G4endl;
    }

    G4Scheduler::Instance()->Initialize();

    fpThreadData->fThreadInitialized = true;

    G4VMoleculeCounter::InitializeInstance();

    InitializeFile();
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::InitializeFile()
{
    if (fVerbose)
    {
        G4cout << "G4DNAChemistryManager::InitializeFile() is called"
               << G4endl;
    }

    if (fpThreadData->fpPhysChemIO)
    {
        fpThreadData->fpPhysChemIO->InitializeFile();
    }
}

//------------------------------------------------------------------------------

G4bool G4DNAChemistryManager::IsActivated()
{
    return fgInstance ? fgInstance->IsChemistryActivated() : false;
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
                                      std::ios_base::openmode mode)
{
    if (fVerbose)
    {
        G4cout << "G4DNAChemistryManager: Write chemical stage into "
               << output.data() << G4endl;
    }

    if (!fpThreadData->fpPhysChemIO)
    {
        fpThreadData->fpPhysChemIO.reset(new G4PhysChemIO::FormattedText());
    }

    fpThreadData->fpPhysChemIO->WriteInto(output, mode);

}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::AddEmptyLineInOutputFile()
{
    if (fpThreadData->fpPhysChemIO)
    {
        fpThreadData->fpPhysChemIO->AddEmptyLineInOutputFile();
    }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::CloseFile()
{
    if (fpThreadData->fpPhysChemIO)
    {
        fpThreadData->fpPhysChemIO->CloseFile();
    }
}

//------------------------------------------------------------------------------

G4DNAWaterExcitationStructure* G4DNAChemistryManager::GetExcitationLevel()
{
    if (!fpExcitationLevel)
    {
        fpExcitationLevel.reset(new G4DNAWaterExcitationStructure);
    }
    return fpExcitationLevel.get();
}

//------------------------------------------------------------------------------

G4DNAWaterIonisationStructure* G4DNAChemistryManager::GetIonisationLevel()
{
    if (!fpIonisationLevel)
    {
        fpIonisationLevel.reset(new G4DNAWaterIonisationStructure);
    }
    return fpIonisationLevel.get();
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::CreateWaterMolecule(ElectronicModification modification,
                                                G4int electronicLevel,
                                                const G4Track* pIncomingTrack)
{
    if (fpThreadData->fpPhysChemIO)
    {
        G4double energy = -1.;

        switch (modification)
        {
        case eDissociativeAttachment:
            energy = 0.;
            break;
        case eExcitedMolecule:
            energy = GetExcitationLevel()->ExcitationEnergy(electronicLevel);
            break;
        case eIonizedMolecule:
            energy = GetIonisationLevel()->IonisationEnergy(electronicLevel);
            break;
        }

        fpThreadData->fpPhysChemIO->CreateWaterMolecule(modification,
                                                        4 - electronicLevel,
                                                        energy,
                                                        pIncomingTrack);
    }

    if (fActiveChemistry)
    {
        G4Molecule* pH2OMolecule = new G4Molecule(G4H2O::Definition());

        switch (modification)
        {
        case eDissociativeAttachment:
            pH2OMolecule->AddElectron(5, 1);
            break;
        case eExcitedMolecule:
            pH2OMolecule->ExciteMolecule(4 - electronicLevel);
            break;
        case eIonizedMolecule:
            pH2OMolecule->IonizeMolecule(4 - electronicLevel);
            break;
        }

        G4Track* pH2OTrack = pH2OMolecule->BuildTrack(picosecond,
                                                      pIncomingTrack->GetPosition());

        pH2OTrack->SetParentID(pIncomingTrack->GetTrackID());
        pH2OTrack->SetTrackStatus(fStopButAlive);
        pH2OTrack->SetKineticEnergy(0.);
        PushTrack(pH2OTrack);
    }
}

//------------------------------------------------------------------------------
// pFinalPosition is optional
void G4DNAChemistryManager::CreateSolvatedElectron(const G4Track* pIncomingTrack,
                                                   G4ThreeVector* pFinalPosition)
{
    if (fpThreadData->fpPhysChemIO)
    {
        fpThreadData->fpPhysChemIO->CreateSolvatedElectron(pIncomingTrack,
                                                           pFinalPosition);
    }

    if (fActiveChemistry)
    {
        PushMolecule(std::unique_ptr<G4Molecule>(new G4Molecule(G4Electron_aq::Definition())),
                     picosecond,
                     pFinalPosition ? *pFinalPosition : pIncomingTrack->GetPosition(),
                     pIncomingTrack->GetTrackID());
    }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::PushMolecule(std::unique_ptr<G4Molecule> pMolecule,
                                         double time,
                                         const G4ThreeVector& position,
                                         int parentID)
{
    assert(fActiveChemistry
           && "To inject chemical species, the chemistry must be activated. "
              "Check chemistry activation before injecting species.");
    G4Track* pTrack = pMolecule->BuildTrack(time, position);
    pTrack->SetTrackStatus(fAlive);
    pTrack->SetParentID(parentID);
    pMolecule.release();
    PushTrack(pTrack);
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::SetGlobalTemperature(G4double temp_K)
{
    G4MolecularConfiguration::SetGlobalTemperature(temp_K);
    G4DNAMolecularReactionTable::Instance()->ScaleReactionRateForNewTemperature(temp_K);
}

//------------------------------------------------------------------------------
// [[deprecated]] : chemistry list should never be nullptr
void G4DNAChemistryManager::SetChemistryList(G4VUserChemistryList* pChemistryList)
{
    fpUserChemistryList.reset(pChemistryList);
    fOwnChemistryList = false;
    SetChemistryActivation(true);
}

void G4DNAChemistryManager::SetChemistryList(G4VUserChemistryList& chemistryList)
{
    fpUserChemistryList.reset(&chemistryList);
    fOwnChemistryList = false;
    SetChemistryActivation(true);
}

void G4DNAChemistryManager::SetChemistryList(std::unique_ptr<G4VUserChemistryList> pChemistryList)
{
    fpUserChemistryList = std::move(pChemistryList);
    fOwnChemistryList = true;
    SetChemistryActivation(true);
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::Deregister(G4VUserChemistryList& chemistryList)
{
    if (fpUserChemistryList.get() == &chemistryList)
    {
        if (!fpUserChemistryList->IsPhysicsConstructor() || fOwnChemistryList)
        {
            fpUserChemistryList.reset();
        }

        fpUserChemistryList.release();
    }
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::PushTrack(G4Track* pTrack)
{
    G4ITTrackHolder::Instance()->Push(pTrack);
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::SetVerbose(G4int verbose)
{
    fVerbose = verbose;
}

//------------------------------------------------------------------------------

G4bool G4DNAChemistryManager::IsCounterResetWhenRunEnds() const
{
    return fResetCounterWhenRunEnds;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::ResetCounterWhenRunEnds(G4bool resetCounterWhenRunEnds)
{
    fResetCounterWhenRunEnds = resetCounterWhenRunEnds;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::ForceRebuildingPhysicsTable()
{
    fPhysicsTableBuilt = false;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::ForceMasterReinitialization()
{
    fMasterInitialized = false;
    InitializeMaster();
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::ForceThreadReinitialization()
{
    fForceThreadReinitialization = true;
}

//------------------------------------------------------------------------------

void G4DNAChemistryManager::TagThreadForReinitialization()
{
    fpThreadData->fThreadInitialized = false;
}
