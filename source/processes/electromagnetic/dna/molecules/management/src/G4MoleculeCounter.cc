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

#include "G4MoleculeCounter.hh"
#include "G4MoleculeTable.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scheduler.hh" // TODO: remove this dependency
#include <iomanip>

using namespace std;

namespace G4{
namespace MoleculeCounter {

bool TimePrecision::operator()(const double& a, const double& b) const
{
    if (std::fabs(a - b) < fPrecision)
    {
        return false;
    }
    else
    {
        return a < b;
    }
}

G4ThreadLocal double TimePrecision::fPrecision = 0.5 * picosecond;
}
}

//------------------------------------------------------------------------------
G4MoleculeCounter* G4MoleculeCounter::Instance()
{
    if (!fpInstance)
    {
        fpInstance = new G4MoleculeCounter();
    }
    return dynamic_cast<G4MoleculeCounter*>(fpInstance);
}

//------------------------------------------------------------------------------

G4MoleculeCounter::G4MoleculeCounter()
{
    fVerbose = 0;
    fCheckTimeIsConsistentWithScheduler = true;
}

//------------------------------------------------------------------------------

G4MoleculeCounter::~G4MoleculeCounter() = default;

//------------------------------------------------------------------------------

void G4MoleculeCounter::Initialize()
{
    auto mol_iterator = G4MoleculeTable::Instance()->GetConfigurationIterator();
    while ((mol_iterator)())
    {
        if (IsRegistered(mol_iterator.value()->GetDefinition()) == false)
        {
            continue;
        }

        fCounterMap[mol_iterator.value()]; // initialize the second map
    }
}

//------------------------------------------------------------------------------

void G4MoleculeCounter::SetTimeSlice(double timeSlice)
{
    G4::MoleculeCounter::TimePrecision::fPrecision = timeSlice;
}

//------------------------------------------------------------------------------

G4bool G4MoleculeCounter::SearchTimeMap(Reactant* molecule)
{
    if (fpLastSearch == nullptr)
    {
        fpLastSearch.reset(new Search());
    }
    else
    {
        if (fpLastSearch->fLowerBoundSet &&
            fpLastSearch->fLastMoleculeSearched->first == molecule)
        {
            return true;
        }
    }

    auto mol_it = fCounterMap.find(molecule);
    fpLastSearch->fLastMoleculeSearched = mol_it;

    if (mol_it != fCounterMap.end())
    {
        fpLastSearch->fLowerBoundTime = fpLastSearch->fLastMoleculeSearched->second
                .end();
        fpLastSearch->fLowerBoundSet = true;
    }
    else
    {
        fpLastSearch->fLowerBoundSet = false;
    }

    return false;
}

//------------------------------------------------------------------------------

int G4MoleculeCounter::SearchUpperBoundTime(double time,
                                            bool sameTypeOfMolecule)
{
    auto mol_it = fpLastSearch->fLastMoleculeSearched;
    if (mol_it == fCounterMap.end())
    {
        return 0;
    }

    NbMoleculeAgainstTime& timeMap = mol_it->second;
    if (timeMap.empty())
    {
        return 0;
    }

    if (sameTypeOfMolecule == true)
    {
        if (fpLastSearch->fLowerBoundSet && fpLastSearch->fLowerBoundTime != timeMap.end())
        {
            if (fpLastSearch->fLowerBoundTime->first < time)
            {
                auto upperToLast = fpLastSearch->fLowerBoundTime;
                upperToLast++;

                if (upperToLast == timeMap.end())
                {
                    return fpLastSearch->fLowerBoundTime->second;
                }

                if (upperToLast->first > time)
                {
                    return fpLastSearch->fLowerBoundTime->second;
                }
            }
        }
    }

    auto up_time_it = timeMap.upper_bound(time);

    if (up_time_it == timeMap.end())
    {
        NbMoleculeAgainstTime::reverse_iterator last_time = timeMap.rbegin();
        return last_time->second;
    }
    if (up_time_it == timeMap.begin())
    {
        return 0;
    }

    up_time_it--;

    fpLastSearch->fLowerBoundTime = up_time_it;
    fpLastSearch->fLowerBoundSet = true;

    return fpLastSearch->fLowerBoundTime->second;
}

//------------------------------------------------------------------------------

int G4MoleculeCounter::GetNMoleculesAtTime(Reactant* molecule,
                                           double time)
{
    G4bool sameTypeOfMolecule = SearchTimeMap(molecule);
    return SearchUpperBoundTime(time, sameTypeOfMolecule);
}

//------------------------------------------------------------------------------

void G4MoleculeCounter::AddAMoleculeAtTime(Reactant* molecule,
                                           G4double time,
                                           const G4ThreeVector* /*position*/,
                                           int number)
{
    if (fDontRegister[molecule->GetDefinition()])
    {
        return;
    }

    if (fVerbose)
    {
        G4cout << "G4MoleculeCounter::AddAMoleculeAtTime : " << molecule->GetName()
               << " at time : " << G4BestUnit(time, "Time") << G4endl;
    }

    auto counterMap_i = fCounterMap.find(molecule);

    if (counterMap_i == fCounterMap.end())
    {
        fCounterMap[molecule][time] = number;
    }
    else if (counterMap_i->second.empty())
    {
        counterMap_i->second[time] = number;
    }
    else
    {
        NbMoleculeAgainstTime::reverse_iterator end = counterMap_i->second.rbegin();

        if (end->first <= time ||
            fabs(end->first - time) <= G4::MoleculeCounter::TimePrecision::fPrecision)
            // Case 1 = new time comes after last recorded data
            // Case 2 = new time is about the same as the last recorded one
        {
            double newValue = end->second + number;
            counterMap_i->second[time] = newValue;
        }
        else
        {
            //      if(fabs(time - G4Scheduler::Instance()->GetGlobalTime()) >
            //         G4Scheduler::Instance()->GetTimeTolerance())
            {
                G4ExceptionDescription errMsg;
                errMsg << "Time of species "
                       << molecule->GetName() << " is "
                       << G4BestUnit(time, "Time") << " while "
                       << " global time is "
                       << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
                       << G4endl;
                G4Exception("G4MoleculeCounter::AddAMoleculeAtTime",
                            "TIME_DONT_MATCH",
                            FatalException, errMsg);
            }
        }
    }
}

//------------------------------------------------------------------------------

void G4MoleculeCounter::RemoveAMoleculeAtTime(const G4MolecularConfiguration* pMolecule,
                                              G4double time,
                                              const G4ThreeVector* /*position*/,
                                              int number)
{
    if (fDontRegister[pMolecule->GetDefinition()])
    {
        return;
    }

    if (fVerbose)
    {
        G4cout << "G4MoleculeCounter::RemoveAMoleculeAtTime : "
               << pMolecule->GetName() << " at time : " << G4BestUnit(time, "Time")
               << G4endl;
    }

    if (fCheckTimeIsConsistentWithScheduler)
    {
        if (fabs(time - G4Scheduler::Instance()->GetGlobalTime()) >
            G4Scheduler::Instance()->GetTimeTolerance())
        {
            G4ExceptionDescription errMsg;
            errMsg << "Time of species "
                   << pMolecule->GetName() << " is "
                   << G4BestUnit(time, "Time") << " while "
                   << " global time is "
                   << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
                   << G4endl;
            G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime",
                        "TIME_DONT_MATCH",
                        FatalException, errMsg);
        }
    }

    NbMoleculeAgainstTime& nbMolPerTime = fCounterMap[pMolecule];

    if (nbMolPerTime.empty())
    {
        pMolecule->PrintState();
        Dump();
        G4String errMsg =
                "You are trying to remove molecule " + pMolecule->GetName() +
                " from the counter while this kind of molecules has not been registered yet";
        G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime", "",
                    FatalErrorInArgument, errMsg);

        return;
    }
    else
    {
        NbMoleculeAgainstTime::reverse_iterator it = nbMolPerTime.rbegin();

        if (it == nbMolPerTime.rend())
        {
            it--;

            G4String errMsg =
                    "There was no " + pMolecule->GetName() + " recorded at the time or even before the time asked";
            G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime", "",
                        FatalErrorInArgument, errMsg);
        }

        if (time - it->first < -G4::MoleculeCounter::TimePrecision::fPrecision)
        {
            Dump();
            G4ExceptionDescription errMsg;
            errMsg << "Is time going back?? " << pMolecule->GetName()
                   << " is being removed at time " << G4BestUnit(time, "Time")
                   << " while last recorded time was "
                   << G4BestUnit(it->first, "Time") << ".";
            G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime",
                        "RETURN_TO_THE_FUTUR",
                        FatalErrorInArgument,
                        errMsg);
        }

        double finalN = it->second - number;

        if (finalN < 0)
        {
            Dump();
            G4ExceptionDescription errMsg;
            errMsg << "After removal of " << number << " species of "
                   << pMolecule->GetName() << " the final number at time "
                   << G4BestUnit(time, "Time") << " is less than zero and so not valid."
                   << " Global time is "
                   << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
                   << ". Previous selected time is "
                   << G4BestUnit(it->first, "Time")
                   << G4endl;
            G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime",
                        "N_INF_0",
                        FatalException, errMsg);
        }

        nbMolPerTime[time] = finalN;
    }
}

//------------------------------------------------------------------------------

G4MoleculeCounter::RecordedMolecules G4MoleculeCounter::GetRecordedMolecules()
{
    if (fVerbose > 1)
    {
        G4cout << "Entering in G4MoleculeCounter::RecordMolecules" << G4endl;
    }

    RecordedMolecules output(new ReactantList());

    for (const auto & it : fCounterMap)
    {
        output->push_back(it.first);
    }
    return output;
}

//------------------------------------------------------------------------------

RecordedTimes G4MoleculeCounter::GetRecordedTimes()
{
    RecordedTimes output(new std::set<G4double>);

    for(const auto& it : fCounterMap)
    {
        for(const auto& it2 : it.second)
        {
            //time = it2->first;
            output->insert(it2.first);
        }
    }

    return output;
}

//------------------------------------------------------------------------------

// >>DEV<<
//void G4MoleculeCounter::SignalReceiver(G4SpeciesInCM* /*speciesInCM*/,
//                                       size_t moleculeID,
//                                       int /*number*/,
//                                       G4SpeciesInCM::SpeciesChange speciesChange,
//                                       int diff)
//{
//  switch(speciesChange)
//  {
//    case G4SpeciesInCM::eAdd:
//      AddAMoleculeAtTime(G4MoleculeTable::Instance()->GetConfiguration((int)moleculeID),
//                         G4Scheduler::Instance()->GetGlobalTime(),
//                         diff);
//      break;
//    case G4SpeciesInCM::eRemove:
//      RemoveAMoleculeAtTime(G4MoleculeTable::Instance()->GetConfiguration((int)moleculeID),
//                         G4Scheduler::Instance()->GetGlobalTime(),
//                         diff);
//      break;
//  }
//}

//------------------------------------------------------------------------------

void G4MoleculeCounter::Dump()
{
    for (const auto& it : fCounterMap)
    {
        auto pReactant = it.first;

        G4cout << " --- > For " << pReactant->GetName() << G4endl;

        for (const auto& it2 : it.second)
        {
            G4cout << " " << G4BestUnit(it2.first, "Time")
                   << "    " << it2.second << G4endl;
        }
    }
}

//------------------------------------------------------------------------------

void G4MoleculeCounter::ResetCounter()
{
    if (fVerbose)
    {
        G4cout << " ---> G4MoleculeCounter::ResetCounter" << G4endl;
    }
    fCounterMap.clear();
    fpLastSearch.reset(0);
}

//------------------------------------------------------------------------------

const NbMoleculeAgainstTime& G4MoleculeCounter::GetNbMoleculeAgainstTime(Reactant* molecule)
{
    return fCounterMap[molecule];
}

//------------------------------------------------------------------------------

void G4MoleculeCounter::SetVerbose(G4int level)
{
    fVerbose = level;
}

//------------------------------------------------------------------------------

G4int G4MoleculeCounter::GetVerbose()
{
    return fVerbose;
}

//------------------------------------------------------------------------------

void G4MoleculeCounter::DontRegister(const G4MoleculeDefinition* molDef)
{
    fDontRegister[molDef] = true;
}

//------------------------------------------------------------------------------

bool G4MoleculeCounter::IsRegistered(const G4MoleculeDefinition* molDef)
{
    if (fDontRegister.find(molDef) == fDontRegister.end())
    {
        return true;
    }
    return false;
}

//------------------------------------------------------------------------------

void G4MoleculeCounter::RegisterAll()
{
    fDontRegister.clear();
}

G4bool G4MoleculeCounter::IsTimeCheckedForConsistency() const
{
    return fCheckTimeIsConsistentWithScheduler;
}

void G4MoleculeCounter::CheckTimeForConsistency(G4bool flag)
{
    fCheckTimeIsConsistentWithScheduler = flag;
}

