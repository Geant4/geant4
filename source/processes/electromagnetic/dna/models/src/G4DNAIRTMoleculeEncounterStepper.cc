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
/*
 * G4DNAIRTMoleculeEncounterStepper.cc
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
*/

#include "G4DNAIRTMoleculeEncounterStepper.hh"
#include "G4VDNAReactionModel.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4H2O.hh"
#include "G4memory.hh"
#include "G4UnitsTable.hh"
#include "G4MoleculeFinder.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Scheduler.hh"
#include "G4ITReaction.hh"

using namespace std;
using namespace CLHEP;

//#define DEBUG_MEM

#ifdef DEBUG_MEM
#include "G4MemStat.hh"
using namespace G4MemStat;
#endif

G4DNAIRTMoleculeEncounterStepper::Utils::Utils(const G4Track& tA,
                                            const G4MolecularConfiguration* pMoleculeB)
    : fpTrackA(tA)
    , fpMoleculeB(pMoleculeB)
{
    fpMoleculeA = GetMolecule(tA);
    fDA = fpMoleculeA->GetDiffusionCoefficient();
    fDB = fpMoleculeB->GetDiffusionCoefficient();
    fConstant = 8 * (fDA + fDB + 2 * sqrt(fDA * fDB));
}

G4DNAIRTMoleculeEncounterStepper::G4DNAIRTMoleculeEncounterStepper()
    : G4VITTimeStepComputer()
    , fHasAlreadyReachedNullTime(false)
    , fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
    , fReactionModel(nullptr)
    , fVerbose(0)
{
	fpTrackContainer = G4ITTrackHolder::Instance();
	fReactionSet = G4ITReactionSet::Instance();
}

G4DNAIRTMoleculeEncounterStepper::~G4DNAIRTMoleculeEncounterStepper() = default;

void G4DNAIRTMoleculeEncounterStepper::Prepare()
{
    fSampledMinTimeStep = DBL_MAX;
    if(G4Scheduler::Instance()->GetGlobalTime() == G4Scheduler::Instance()->GetStartTime()){
        G4VITTimeStepComputer::Prepare();
        G4MoleculeFinder::Instance()->UpdatePositionMap();
    }
}

void G4DNAIRTMoleculeEncounterStepper::InitializeForNewTrack()
{
    if (fReactants)
    {
        fReactants.reset();
    }
    fSampledMinTimeStep = DBL_MAX;
    fHasAlreadyReachedNullTime = false;
}

template<typename T>
inline bool IsInf(T value)
{
    return std::numeric_limits<T>::has_infinity
        && value == std::numeric_limits<T>::infinity();
}

G4double
G4DNAIRTMoleculeEncounterStepper::CalculateStep(const G4Track& trackA,
                                             const G4double& userMinTimeStep)
{

    auto pMoleculeA = GetMolecule(trackA);
    InitializeForNewTrack();
    fUserMinTimeStep = userMinTimeStep;

#ifdef G4VERBOSE
    if (fVerbose)
    {
        G4cout
            << "_______________________________________________________________________"
            << G4endl;
        G4cout << "G4DNAMoleculeEncounterStepper::CalculateStep" << G4endl;
        G4cout << "Check done for molecule : " << pMoleculeA->GetName()
            << " (" << trackA.GetTrackID() << ") "
            << G4endl;
    }
#endif

    //__________________________________________________________________
    // Retrieve general informations for making reactions
    auto pMolConfA = pMoleculeA->GetMolecularConfiguration();

    const auto pReactantList = fMolecularReactionTable->CanReactWith(pMolConfA);

    if (!pReactantList)
    {
#ifdef G4VERBOSE
        //    DEBUG
        if (fVerbose > 1)
        {
            G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
            G4cout << "!!! WARNING" << G4endl;
            G4cout << "G4MoleculeEncounterStepper::CalculateStep will return infinity "
                "for the reaction because the molecule "
                << pMoleculeA->GetName()
                << " does not have any reactants given in the reaction table."
                << G4endl;
            G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
        }
#endif
        return DBL_MAX;
    }

    G4int nbReactives = (G4int)pReactantList->size();

    if (nbReactives == 0)
    {
#ifdef G4VERBOSE
        //    DEBUG
        if (fVerbose)
        {
            // TODO replace with the warning mode of G4Exception
            G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
            G4cout << "!!! WARNING" << G4endl;
            G4cout << "G4MoleculeEncounterStepper::CalculateStep will return infinity "
                "for the reaction because the molecule "
                << pMoleculeA->GetName()
                << " does not have any reactants given in the reaction table."
                << "This message can also result from a wrong implementation of the reaction table."
                << G4endl;
            G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
        }
#endif
        return DBL_MAX;
    }

    fReactants.reset(new vector<G4Track*>());
    fReactionModel->Initialise(pMolConfA, trackA);

    //__________________________________________________________________
    // Start looping on possible reactants
    for (G4int i = 0; i < nbReactives; ++i)
    {
        auto pMoleculeB = (*pReactantList)[i];

        //______________________________________________________________
        // Retrieve reaction range
        const G4double R = fReactionModel->GetReactionRadius(i);

        //______________________________________________________________
        // Use KdTree algorithm to find closest reactants
        G4KDTreeResultHandle resultsNearest(
            G4MoleculeFinder::Instance()->FindNearest(pMoleculeA,
                                                      pMoleculeB->GetMoleculeID()));

        if (resultsNearest == 0) continue;

        G4double r2 = resultsNearest->GetDistanceSqr();
        Utils utils(trackA, pMoleculeB);

        if (r2 <= R * R) // ==> Record in range
        {
            // Entering in this condition may due to the fact that molecules are very close
            // to each other
            // Therefore, if we only take the nearby reactant into account, it might have already
            // reacted. Instead, we will take all possible reactants that satisfy the condition r<R

            if (fHasAlreadyReachedNullTime == false)
            {
                fReactants->clear();
                fHasAlreadyReachedNullTime = true;
            }

            fSampledMinTimeStep = 0.;
            G4KDTreeResultHandle resultsInRange(
                G4MoleculeFinder::Instance()->FindNearestInRange(pMoleculeA,
                                                                 pMoleculeB->GetMoleculeID(),
                                                                 R));
            CheckAndRecordResults(utils,
#ifdef G4VERBOSE
                                  R,
#endif
                                  resultsInRange);
        }
        else
        {
            G4double r = sqrt(r2);
            G4double tempMinET = pow(r - R, 2) / utils.fConstant;
            // constant = 16 * (fDA + fDB + 2*sqrt(fDA*fDB))

            if (tempMinET <= fSampledMinTimeStep)
            {
                if (fUserMinTimeStep < DBL_MAX/*IsInf(fUserMinTimeStep) == false*/
                    && tempMinET <= fUserMinTimeStep) // ==> Record in range
                {
                    if (fSampledMinTimeStep > fUserMinTimeStep)
                    {
                        fReactants->clear();
                    }

                    fSampledMinTimeStep = fUserMinTimeStep;

                    G4double range = R + sqrt(fUserMinTimeStep*utils.fConstant);

                    G4KDTreeResultHandle resultsInRange(
                        G4MoleculeFinder::Instance()->
                        FindNearestInRange(pMoleculeA,
                                           pMoleculeB->GetMoleculeID(),
                                           range));

                    CheckAndRecordResults(utils,
#ifdef G4VERBOSE
                                          range,
#endif
                                          resultsInRange);
                }
                else // ==> Record nearest
                {
                    if (tempMinET < fSampledMinTimeStep)
                        // to avoid cases where fSampledMinTimeStep == tempMinET
                    {
                        fSampledMinTimeStep = tempMinET;
                        fReactants->clear();
                    }

                    CheckAndRecordResults(utils,
#ifdef G4VERBOSE
                                          R,
#endif
                                          resultsNearest);
                }
            }
        }
    }

#ifdef G4VERBOSE
    if (fVerbose)
    {
        G4cout << "G4MoleculeEncounterStepper::CalculateStep will finally return :"
            << G4BestUnit(fSampledMinTimeStep, "Time") << G4endl;

        if (fVerbose > 1)
        {
            G4cout << "Selected reactants for trackA: " << pMoleculeA->GetName()
                << " (" << trackA.GetTrackID() << ") are: ";

            vector<G4Track*>::iterator it;
            for (it = fReactants->begin(); it != fReactants->end(); it++)
            {
                G4Track* trackB = *it;
                G4cout << GetMolecule(trackB)->GetName() << " ("
                    << trackB->GetTrackID() << ") \t ";
            }
            G4cout << G4endl;
        }
    }
#endif
    return fSampledMinTimeStep;
}




void G4DNAIRTMoleculeEncounterStepper::CheckAndRecordResults(const Utils& utils,
#ifdef G4VERBOSE
                                                          const G4double R,
#endif
                                                          G4KDTreeResultHandle& results)
{
    if (results == 0)
    {
#ifdef G4VERBOSE
        if (fVerbose > 1)
        {
            G4cout << "No molecule " << utils.fpMoleculeB->GetName()
                << " found to react with " << utils.fpMoleculeA->GetName()
                << G4endl;
        }
#endif
        return;
    }

    for (results->Rewind(); !results->End(); results->Next())
    {
        G4IT* reactiveB = results->GetItem<G4IT>();

        if (reactiveB == 0)
        {
            continue;
        }

        G4Track *trackB = reactiveB->GetTrack();

        if (trackB == 0)
        {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription
                << "The reactant B found using the MoleculeFinder does not have a valid "
                "track attached to it. If this is done on purpose, please do "
                "not record this molecule in the MoleculeFinder."
                << G4endl;
            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults",
                        "MoleculeEncounterStepper001", FatalErrorInArgument,
                        exceptionDescription);
            continue;
        }

        if (trackB->GetTrackStatus() != fAlive)
        {
            continue;
        }

        if (trackB == &utils.fpTrackA)
        {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription
                << "A track is reacting with itself (which is impossible) ie fpTrackA == trackB"
                << G4endl;
            exceptionDescription << "Molecule A (and B) is of type : "
                << utils.fpMoleculeA->GetName() << " with trackID : "
                << utils.fpTrackA.GetTrackID() << G4endl;

            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults",
                        "MoleculeEncounterStepper003", FatalErrorInArgument,
                        exceptionDescription);

        }

        if (fabs(trackB->GetGlobalTime() - utils.fpTrackA.GetGlobalTime())
        > utils.fpTrackA.GetGlobalTime() * (1 - 1 / 100))
        {
            // DEBUG
            G4ExceptionDescription exceptionDescription;
            exceptionDescription
                << "The interacting tracks are not synchronized in time" << G4endl;
            exceptionDescription
                << "trackB->GetGlobalTime() != fpTrackA.GetGlobalTime()" << G4endl;

            exceptionDescription << "fpTrackA : trackID : " << utils.fpTrackA.GetTrackID()
                << "\t Name :" << utils.fpMoleculeA->GetName()
                << "\t fpTrackA->GetGlobalTime() = "
                << G4BestUnit(utils.fpTrackA.GetGlobalTime(), "Time") << G4endl;

            exceptionDescription << "trackB : trackID : " << trackB->GetTrackID()
                << "\t Name :" << utils.fpMoleculeB->GetName()
                << "\t trackB->GetGlobalTime() = "
                << G4BestUnit(trackB->GetGlobalTime(), "Time") << G4endl;

            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults",
                        "MoleculeEncounterStepper004", FatalErrorInArgument,
                        exceptionDescription);
        }

#ifdef G4VERBOSE
        if (fVerbose > 1)
        {

            G4double r2 = results->GetDistanceSqr();
            G4cout << "\t ************************************************** " << G4endl;
            G4cout << "\t Reaction between "
                << utils.fpMoleculeA->GetName() << " (" << utils.fpTrackA.GetTrackID() << ") "
                << " & " << utils.fpMoleculeB->GetName() << " (" << trackB->GetTrackID() << "), "
                << "Interaction Range = "
                << G4BestUnit(R, "Length") << G4endl;
            G4cout << "\t Real distance between reactants  = "
                << G4BestUnit((utils.fpTrackA.GetPosition() - trackB->GetPosition()).mag(), "Length") << G4endl;
            G4cout << "\t Distance between reactants calculated by nearest neighbor algorithm = "
                << G4BestUnit(sqrt(r2), "Length") << G4endl;

        }
#endif

        fReactants->push_back(trackB);
    }
}

void G4DNAIRTMoleculeEncounterStepper::SetReactionModel(G4VDNAReactionModel* pReactionModel)
{
    fReactionModel = pReactionModel;
}

G4VDNAReactionModel* G4DNAIRTMoleculeEncounterStepper::GetReactionModel()
{
    return fReactionModel;
}

void G4DNAIRTMoleculeEncounterStepper::SetVerbose(int flag)
{
    fVerbose = flag;
}

G4double G4DNAIRTMoleculeEncounterStepper::CalculateMinTimeStep(G4double currentGlobalTime, G4double definedMinTimeStep){

    G4bool start = true;
    G4bool active = false;

    fUserMinTimeStep = definedMinTimeStep;

    if(fReactionSet->Empty()){
        if(currentGlobalTime == G4Scheduler::Instance()->GetStartTime()){

            for (auto pTrack : *fpTrackContainer->GetMainList())
            {
                if (pTrack == nullptr)
                {
                    G4ExceptionDescription exceptionDescription;
                    exceptionDescription << "No track found.";
                    G4Exception("G4Scheduler::CalculateMinStep", "ITScheduler006",
                                FatalErrorInArgument, exceptionDescription);
                    continue;
                }

                G4TrackStatus trackStatus = pTrack->GetTrackStatus();
                if (trackStatus == fStopAndKill || trackStatus == fStopButAlive)
                {
                    start = false;
                    continue;
                }
                active = true;
            }

            if(start == true){
                return -1;
            }else if(active == false){
                G4Scheduler::Instance()->Stop();
                return fSampledMinTimeStep;
            }else{
                return fSampledMinTimeStep;
            }

        }else{
            for (auto pTrack : *fpTrackContainer->GetMainList())
            {
                pTrack->SetGlobalTime(G4Scheduler::Instance()->GetEndTime());
            }
            return fSampledMinTimeStep;
        }
    }

    auto fReactionSetInTime = fReactionSet->GetReactionsPerTime();
    fSampledMinTimeStep = fReactionSetInTime.begin()->get()->GetTime() - currentGlobalTime;

    return fSampledMinTimeStep;
}
