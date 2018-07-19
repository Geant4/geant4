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
// $Id: G4MoleculeCounter.cc 101354 2016-11-15 08:27:51Z gcosmo $
//

#include <iomanip>
#include "G4MoleculeCounter.hh"
#include "G4MoleculeTable.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scheduler.hh" // TODO: remove this dependency

using namespace std;

G4ThreadLocal double compDoubleWithPrecision::fPrecision = 0;

//------------------------------------------------------------------------------
G4MoleculeCounter* G4MoleculeCounter::Instance(){
  if (!fpInstance) fpInstance = new G4MoleculeCounter();
  return dynamic_cast<G4MoleculeCounter*>(fpInstance);
}

//------------------------------------------------------------------------------

G4MoleculeCounter::G4MoleculeCounter()
{
  fVerbose = 0;
  fCheckTimeIsConsistentWithScheduler = true;
  if(compDoubleWithPrecision::fPrecision == 0)
  {
    compDoubleWithPrecision::fPrecision = 0.5*picosecond;
  }
}

//------------------------------------------------------------------------------

G4MoleculeCounter::~G4MoleculeCounter()
{
}
  
//------------------------------------------------------------------------------

void G4MoleculeCounter::Initialize()
{
//  G4cout << "G4MoleculeCounter::Initialize" << G4endl;
  
  G4ConfigurationIterator mol_iterator = G4MoleculeTable::Instance()
      ->GetConfigurationIterator();
  while ((mol_iterator)())
  {
    if(IsRegistered(mol_iterator.value()->GetDefinition()) == false)
    {
      continue;
    }

    //    G4cout << "G4MoleculeCounter::Initialize" << G4endl;
    //    G4cout << mol_iterator->value()->GetName() << G4endl;
    fCounterMap[mol_iterator.value()]; // initialize the second map
  }
}

//------------------------------------------------------------------------------

void G4MoleculeCounter::SetTimeSlice(double timeSlice)
{
  compDoubleWithPrecision::fPrecision = timeSlice;
}

//------------------------------------------------------------------------------

G4bool G4MoleculeCounter::SearchTimeMap(G4MolecularConfiguration* molecule)
{
  if (fpLastSearch.get() == 0)
  {
    fpLastSearch.reset(new Search());
  }
  else
  {
    if (fpLastSearch->fLowerBoundSet &&
        fpLastSearch->fLastMoleculeSearched->first == molecule)
      return true;
  }

  CounterMapType::iterator mol_it = fCounterMap.find(molecule);
  fpLastSearch->fLastMoleculeSearched = mol_it;

  if (mol_it != fCounterMap.end()) // TODO
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
  CounterMapType::iterator mol_it = fpLastSearch->fLastMoleculeSearched;
  if (mol_it == fCounterMap.end()) return 0; // RETURN

  NbMoleculeAgainstTime& timeMap = mol_it->second;
  if (timeMap.empty()) return 0;

  NbMoleculeAgainstTime::iterator end_time = timeMap.end();

  if (sameTypeOfMolecule == true)
  {
    //G4cout << "SAME MOLECULE" << G4endl;
    if (fpLastSearch->fLowerBoundSet && fpLastSearch->fLowerBoundTime
        != end_time)
    {
      if (fpLastSearch->fLowerBoundTime->first < time)
      {
        NbMoleculeAgainstTime::iterator upperToLast = fpLastSearch
            ->fLowerBoundTime;
        upperToLast++;

        if (upperToLast == end_time)
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
  /*
   else
   {
   G4cout << "--> Molecule has changed" << G4endl;
   }
   */
  //G4cout << "Searching" << G4endl;
  // With upper bound
  NbMoleculeAgainstTime::iterator up_time_it = timeMap.upper_bound(time);

  if (up_time_it == end_time)
  {
    NbMoleculeAgainstTime::reverse_iterator last_time = timeMap.rbegin();

//    {
      //G4cout << "RETURN LAST : " << G4BestUnit(time, "Time") << G4endl;
      return last_time->second;
//    }

//    G4cout << "RETURN 0 (1)" << G4endl;
//    return 0; // RETURN
  }
  if (up_time_it == timeMap.begin())
  {
//    G4cout << "RETURN 0 (2)" << G4endl;
    return 0; // RETURN
  }

  //G4cout << "Going back : " << up_time_it->first << "-->";

  up_time_it--;

//  G4cout << up_time_it->first << G4endl;

  fpLastSearch->fLowerBoundTime = up_time_it;
  fpLastSearch->fLowerBoundSet = true;

//  G4cout << "returning : " << fpLastSearch->fLowerBoundTime->second << G4endl;

  return fpLastSearch->fLowerBoundTime->second;
}

//------------------------------------------------------------------------------

int G4MoleculeCounter::GetNMoleculesAtTime(G4MolecularConfiguration* molecule,
                                           double time)
{
  G4bool sameTypeOfMolecule = SearchTimeMap(molecule);
  return SearchUpperBoundTime(time, sameTypeOfMolecule);
}

//------------------------------------------------------------------------------

void G4MoleculeCounter::AddAMoleculeAtTime(G4MolecularConfiguration* molecule,
                                           G4double time,
                                           const G4ThreeVector* /*position*/,
                                           int number)
{
  if (fDontRegister[molecule->GetDefinition()]) return;

  if (fVerbose){
    G4cout << "G4MoleculeCounter::AddAMoleculeAtTime : " << molecule->GetName()
           << " at time : " << G4BestUnit(time, "Time") << G4endl;
  }

  CounterMapType::iterator counterMap_i =
      fCounterMap.find(molecule);

  if (counterMap_i == fCounterMap.end()){
    fCounterMap[molecule][time] = number;
  }
  else if (counterMap_i->second.empty()){
    counterMap_i->second[time] = number;
  }
  else{
    NbMoleculeAgainstTime::reverse_iterator end = counterMap_i->second.rbegin();

    if (end->first <= time ||
        fabs(end->first - time) <= compDoubleWithPrecision::fPrecision)
      // Case 1 = new time comes after last recorded data
      // Case 2 = new time is about the same as the last recorded one
    {
      double newValue =  end->second + number;
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
        G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime",
                    "TIME_DONT_MATCH",
                    FatalException, errMsg);
      }
    }
  }
}

//------------------------------------------------------------------------------

void
G4MoleculeCounter::RemoveAMoleculeAtTime(G4MolecularConfiguration* molecule,
                                         G4double time,
                                         const G4ThreeVector* /*position*/,
                                         int number)
{
  if (fDontRegister[molecule->GetDefinition()]) return;

  if (fVerbose)
  {
    G4cout << "G4MoleculeCounter::RemoveAMoleculeAtTime : "
           << molecule->GetName() << " at time : " << G4BestUnit(time, "Time")
           << G4endl;
  }

  if(fCheckTimeIsConsistentWithScheduler)
  {
    if(fabs(time - G4Scheduler::Instance()->GetGlobalTime()) >
       G4Scheduler::Instance()->GetTimeTolerance())
    {
      G4ExceptionDescription errMsg;
      errMsg << "Time of species "
          << molecule->GetName() << " is "
          << G4BestUnit(time, "Time") << " while "
          << " global time is "
          << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
          << G4endl;
      G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime",
                  "TIME_DONT_MATCH",
                  FatalException, errMsg);
    }
  }

  NbMoleculeAgainstTime& nbMolPerTime = fCounterMap[molecule];

  if (nbMolPerTime.empty())
  {
    molecule->PrintState();
    Dump();
    G4String errMsg =
        "You are trying to remove molecule " + molecule->GetName()
        + " from the counter while this kind of molecules has not been registered yet";
    G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime", "",
                FatalErrorInArgument, errMsg);

    return;
  }
  else
  {
    NbMoleculeAgainstTime::reverse_iterator it = nbMolPerTime.rbegin();
    
    if (it == nbMolPerTime.rend()){
      it--;
      
      G4String errMsg = "There was no " + molecule->GetName()
      + " recorded at the time or even before the time asked";
      G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime", "",
                  FatalErrorInArgument, errMsg);
    }

    if (time - it->first < -compDoubleWithPrecision::fPrecision){
      Dump();
      G4ExceptionDescription errMsg;
      errMsg << "Is time going back?? " << molecule->GetName()
             << " is being removed at time " << G4BestUnit(time, "Time")
             << " while last recorded time was "
             << G4BestUnit(it->first, "Time") << ".";
      G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime",
                  "RETURN_TO_THE_FUTUR",
                  FatalErrorInArgument,
                  errMsg);
    }

    double finalN = it->second - number;

    if(finalN < 0){
      Dump();
      G4ExceptionDescription errMsg;
      errMsg << "After removal of " << number << " species of "
          << molecule->GetName() << " the final number at time "
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

  CounterMapType::iterator it;
  RecordedMolecules output (new vector<G4MolecularConfiguration*>);

  for(it = fCounterMap.begin(); it != fCounterMap.end(); it++)
  {
    output->push_back(it->first);
  }
  return output;
}

//------------------------------------------------------------------------------

RecordedTimes G4MoleculeCounter::GetRecordedTimes()
{
  RecordedTimes output(new std::set<G4double>);

  //G4double time;

  CounterMapType::iterator it;
  CounterMapType::const_iterator ite;

  NbMoleculeAgainstTime::iterator it2;
  NbMoleculeAgainstTime::const_iterator ite2;

  // iterate on each molecule
  for (it = fCounterMap.begin(), ite = fCounterMap.end(); it != ite; ++it)
  {
    // iterate on each time
    for (it2 = (it->second).begin(), ite2 = (it->second).end(); it2 != ite2;
        ++it2)
    {
      //time = it2->first;
      output->insert(it2->first);
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
  CounterMapType::iterator it = fCounterMap.begin();
  CounterMapType::iterator end = fCounterMap.end();

  for(;it!=end;++it)
  {
    G4MolecularConfiguration* molConf = it->first;

    G4cout << " --- > For " << molConf->GetName() << G4endl;
    NbMoleculeAgainstTime::iterator it2 = it->second.begin();
    NbMoleculeAgainstTime::iterator end2 = it->second.end();

    for(;it2!=end2;++it2)
    {
      G4cout << " " << G4BestUnit(it2->first, "Time")
             << "    " << it2->second << G4endl;
    }
  }
}
