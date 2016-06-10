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
// $Id: G4MoleculeCounter.cc 84858 2014-10-21 16:08:22Z gcosmo $
//
#include "G4MoleculeCounter.hh"
#include "G4MoleculeTable.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

using namespace std;

G4ThreadLocal double compDoubleWithPrecision::fPrecision = 5e-13;
G4bool G4MoleculeCounter::fUse = FALSE;

G4MoleculeCounter::G4MoleculeCounter()
{
  fVerbose = 0;
}

G4ThreadLocal G4MoleculeCounter* G4MoleculeCounter::fpInstance = 0;

G4MoleculeCounter* G4MoleculeCounter::GetMoleculeCounter()
{
  if (!fpInstance) fpInstance = new G4MoleculeCounter();

  return fpInstance;
}

G4MoleculeCounter* G4MoleculeCounter::Instance()
{
  if (!fpInstance) fpInstance = new G4MoleculeCounter();

  return fpInstance;
}

void G4MoleculeCounter::DeleteInstance()
{
  if (fpInstance)
  {
    delete fpInstance;
    fpInstance = 0;
  }
}

void G4MoleculeCounter::InitializeInstance()
{
  if(fpInstance) fpInstance->Initialize();
}

void G4MoleculeCounter::Initialize()
{
  G4MoleculeModelIterator mol_iterator = G4MoleculeTable::Instance()
      ->GetModelIterator();
  while ((mol_iterator)())
  {
    //    G4cout << "G4MoleculeCounter::Initialize" << G4endl;
    //    G4cout << mol_iterator->value()->GetName() << G4endl;
    fCounterMap[*(mol_iterator.value())]; // initialize the second map
  }
}

void G4MoleculeCounter::SetTimeSlice(double timeSlice)
{
  compDoubleWithPrecision::fPrecision = timeSlice;
}

G4bool G4MoleculeCounter::SearchTimeMap(const G4Molecule &molecule)
{
  if (fpLastSearch.get() == 0)
  {
    fpLastSearch.reset(new Search());
  }
  else
  {
    if (fpLastSearch->fLowerBoundSet && fpLastSearch->fLastMoleculeSearched
        ->first
                                        == molecule) return true;
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
      //G4cout << fpLastSearch->fLowerBoundTime->first << G4endl;
//      G4cout << "fpLastSearch->fLowerBoundTime != timeMap.end() " << time << G4endl;
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
   G4cout << "\n" << G4endl;
   G4cout << "Molecule has changed" << G4endl;
   G4cout << "\n" << G4endl;
   }
   */
  //G4cout << "Searching" << G4endl;
  // With upper bound
  NbMoleculeAgainstTime::iterator up_time_it = timeMap.upper_bound(time);

  if (up_time_it == end_time)
  {
    NbMoleculeAgainstTime::reverse_iterator last_time = timeMap.rbegin();

    if (last_time->first <= time)
    {
      //G4cout << "RETURN LAST : " << G4BestUnit(time, "Time") << G4endl;
      return last_time->second;
    }

//    G4cout << "RETURN 0 (1)" << G4endl;
    return 0; // RETURN
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

int G4MoleculeCounter::GetNMoleculesAtTime(const G4Molecule &molecule,
                                           double time)
{
  G4bool sameTypeOfMolecule = SearchTimeMap(molecule);
  return SearchUpperBoundTime(time, sameTypeOfMolecule);
  // NbMoleculeAgainstTime::iterator low_time_it = timeMap.lower_bound(time);
}

void G4MoleculeCounter::AddAMoleculeAtTime(const G4Molecule& molecule,
                                           G4double time)
{
  if (fDontRegister[molecule.GetDefinition()]) return;

  if (fVerbose)
  {
    G4cout << "G4MoleculeCounter::AddAMoleculeAtTime : " << molecule.GetName()
           << " at time : " << G4BestUnit(time, "Time") << G4endl;
  }

  CounterMapType::iterator counterMap_i = fCounterMap.find(molecule);

  if (counterMap_i == fCounterMap.end())
  {
    // DEBUG
    // if(fVerbose)  G4cout << " !! ***** Map is empty " << G4endl;
    fCounterMap[molecule][time] = 1;
  }
  else if (counterMap_i->second.empty())
  {
    // DEBUG
    // if(fVerbose)  G4cout << " !! ***** Map is empty " << G4endl;
    counterMap_i->second[time] = 1;
  }
  else
  {
    NbMoleculeAgainstTime::iterator end = counterMap_i->second.end();
    end--;

    // DEBUG
    // if(fVerbose)
    // G4cout<<"!! End Time = "<< G4BestUnit(end->first, "Time") <<G4endl;

    if (end->first <= time)
    {
      counterMap_i->second[time] = end->second + 1;
    }
    else
    {
      NbMoleculeAgainstTime::iterator it = counterMap_i->second.lower_bound(
          time);

      while (it->first > time && it != counterMap_i->second.begin())
      {
        // DEBUG
        // if(fVerbose)
        // G4cout<<"!!  ********** Is going back!!!!"<<G4endl;
        it--;
      }

      if (it == counterMap_i->second.begin() && it->first > time)
      {
        // DEBUG
        // if(fVerbose)
        // G4cout<<"!!  ********** Illegal !!!!"<<G4endl;
        return;
      }

      // DEBUG
      // if(fVerbose)
      // {
      //   G4cout<<"!! PREVIOUS NB = "<< it->second <<G4endl;
      //   G4cout<<"!! PREVIOUS TIME = "<< G4BestUnit(it->first,"Time") <<G4endl;
      // }
      counterMap_i->second[time] = it->second + 1;
    }
  }

  // DEBUG
  // if(fVerbose)
  // G4cout<<"!! NB = "<< fCounterMap[molecule][time]<<G4endl;
}

void G4MoleculeCounter::RemoveAMoleculeAtTime(const G4Molecule& molecule,
                                              G4double time)
{
  if (fDontRegister[molecule.GetDefinition()]) return;

  if (fVerbose)
  {
    G4cout << "G4MoleculeCounter::RemoveAMoleculeAtTime : "
           << molecule.GetName() << " at time : " << G4BestUnit(time, "Time")
           << G4endl;
  }

  NbMoleculeAgainstTime& nbMolPerTime = fCounterMap[molecule];

  if (nbMolPerTime.empty())
  {
    molecule.PrintState();
    G4String errMsg =
        "You are trying to remove molecule " + molecule.GetName()
        + " from the counter while this kind of molecules has not been registered yet";
    G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime", "",
                FatalErrorInArgument, errMsg);

    return;
  }
  else
  {
    NbMoleculeAgainstTime::iterator it;

    if (nbMolPerTime.size() == 1)
    {
      it = nbMolPerTime.begin();
      // DEBUG
      // if(fVerbose)
      // G4cout << "!! fCounterMap[molecule].size() == 1" << G4endl;
    }
    else
    {
      it = nbMolPerTime.lower_bound(time);
    }

    if (it == nbMolPerTime.end())
    {
      // DEBUG
      // if(fVerbose)
      // G4cout << " ********** NO ITERATOR !!!!!!!!! " << G4endl;
      it--;

      if (time < it->first)
      {
        G4String errMsg = "There was no " + molecule.GetName()
                          + " record at the time or even before the time asked";
        G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime", "",
                    FatalErrorInArgument, errMsg);
      }
    }

    // DEBUG
    // if(fVerbose)
    // {
    //// G4cout << "G4MoleculeCounter::RemoveAMoleculeAtTime " << G4endl;
    //   G4cout<<"!! Molecule = " << molecule.GetName() << G4endl;
    //   G4cout<<"!! At Time = "<< G4BestUnit(time,"Time") <<G4endl;
    //   G4cout<<"!! PREVIOUS TIME = "<< G4BestUnit(it->first,"Time")<<G4endl;
    //   G4cout<<"!! PREVIOUS Nb = "<< it->second <<G4endl;
    // }

    // If valgrind problem on the line below, it means that the pointer "it"
    // points nowhere
    if (nbMolPerTime.value_comp()(*it, *nbMolPerTime.begin()))
    {
      // DEBUG
      // if(fVerbose)
      //  G4cout<<"!! ***** In value_comp ... " << G4endl;
      it++;
      if (time < it->first)
      {
        G4String errMsg = "There was no " + molecule.GetName()
                          + " record at the time or even before the time asked";
        G4Exception("G4MoleculeCounter::RemoveAMoleculeAtTime", "",
                    FatalErrorInArgument, errMsg);
      }
    }

    while (it->first - time > compDoubleWithPrecision::fPrecision
        && it != nbMolPerTime.begin())
    {
      // DEBUG
      // if(fVerbose)
      // {
      //   G4cout<<"!! ***** Is going back!!!!"<<G4endl;
      //   G4cout<<"!! PREVIOUS TIME = "<< G4BestUnit(it-> first,"Time") <<G4endl;
      // }
      it--;
    }

    if (it == nbMolPerTime.begin() && it->first > time)
    {
      // DEBUG
      //            if(fVerbose)
      //                G4cout<<"!!  ********** Illegal !!!!"<<G4endl;
      return;
    }

    // DEBUG
    // if(fVerbose)
    // {
    //   G4cout<<"!! PREVIOUS NB = "<< (*it).second <<G4endl;
    //   G4cout<<"!! PREVIOUS TIME = "<< G4BestUnit(it->first,"Time")<<G4endl;
    // }
    nbMolPerTime[time] = it->second - 1;
  }

  // DEBUG
  // if(fVerbose)
  // {
  //   G4cout<<"!! NB = "<< nbMolPerTime[time]<<G4endl;
  // }
}

G4MoleculeCounter::RecordedMolecules G4MoleculeCounter::GetRecordedMolecules()
{
  if (fVerbose > 1)
  {
    G4cout << "Entering in G4MoleculeCounter::RecordMolecules" << G4endl;
  }

  CounterMapType::iterator it;
  RecordedMolecules output (new vector<G4Molecule>);

  for(it = fCounterMap.begin(); it != fCounterMap.end(); it++)
  {
    output->push_back(it->first);
  }
  return output;
}

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
