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
// Author: Christian Velten (2025)

#ifndef G4VUSERMOLECULEREACTIONCOUNTER_HH
#define G4VUSERMOLECULEREACTIONCOUNTER_HH 1

#include "G4DNAChemistryManager.hh"
#include "G4MoleculeCounterTemplates.hh"
#include "G4Scheduler.hh"
#include "G4UnitsTable.hh"
#include "G4VMoleculeReactionCounter.hh"

//------------------------------------------------------------------------------

template<class TIndex>
class G4VUserMoleculeReactionCounter : public G4VMoleculeReactionCounter
{
    static_assert(std::is_base_of<G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex, TIndex>::value,
      "TIndex must be derived from G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex! "
      "No forward declaration is allowed.");

  protected:
    struct Search;

  public:
    G4VUserMoleculeReactionCounter();
    G4VUserMoleculeReactionCounter(const G4String&,
                                   MoleculeReactionCounterType = MoleculeReactionCounterType::Basic);
    ~G4VUserMoleculeReactionCounter() override = default;

  public:
    void Initialize() final;
    void InitializeUser() override = 0;
    void ResetCounter() override;
    void Dump() const override;
    void DumpCounterMapIndices() const override;

    void AbsorbCounter(const G4VMoleculeCounterInternalBase*) override;

    std::unique_ptr<G4VMoleculeReactionCounterIndex> BuildSimpleIndex(const G4DNAMolecularReactionData*) const override = 0;

    void RecordReaction(std::unique_ptr<G4VMoleculeReactionCounterIndex>, G4double, G4int = 1) override;

    std::set<const G4DNAMolecularReactionData*> GetRecordedReactions() const override;
    std::set<G4double> GetRecordedTimes() const override;

  protected:
    std::map<TIndex, InnerCounterMapType> fCounterMap{};

  public:
    const std::map<TIndex, InnerCounterMapType>& GetCounterMap() const { return fCounterMap; }
    std::vector<TIndex> GetMapIndices() const;

    virtual G4int GetNbReactionsAtTime(const TIndex&, G4double) const;
    virtual G4int GetNbReactionsAtTime(Search&, const TIndex&, G4double) const;
    virtual std::vector<G4int> GetNbReactionsAtTimes(const TIndex&,
                                                     const std::vector<G4double>&) const;

    //-SEARCH-----------------------------------------------------------------------
  protected:
    struct Search
    {
        Search() : fLowerBoundSet(false) {}
        typename std::map<TIndex, InnerCounterMapType>::const_iterator fLastIndexSearched;
        InnerCounterMapType::const_iterator fLowerBoundTime;
        G4bool fLowerBoundSet;
    };
    G4bool SearchIndexUpdated(Search&, const TIndex&) const;
    G4int SearchUpperBoundTime(Search&, G4double, G4bool) const;
};

//------------------------------------------------------------------------------

// #include "G4VUserMoleculeReactionCounter.icc"

//------------------------------------------------------------------------------

template<typename T>
G4VUserMoleculeReactionCounter<T>::G4VUserMoleculeReactionCounter() : G4VMoleculeReactionCounter()
{}

//------------------------------------------------------------------------------

template<typename T>
G4VUserMoleculeReactionCounter<T>::G4VUserMoleculeReactionCounter(const G4String& name,
                                                                  MoleculeReactionCounterType type)
  : G4VMoleculeReactionCounter(name, type)
{}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeReactionCounter<TIndex>::Initialize()
{
  InitializeUser();
  fIsInitialized = true;
}

//------------------------------------------------------------------------------

template<typename TIndex>
G4int G4VUserMoleculeReactionCounter<TIndex>::GetNbReactionsAtTime(const TIndex& index, G4double time) const
{
  Search search = {};
  return GetNbReactionsAtTime(search, index, time);
}

//------------------------------------------------------------------------------

template<typename TIndex>
G4int G4VUserMoleculeReactionCounter<TIndex>::GetNbReactionsAtTime(Search& search,
                                                                   const TIndex& index,
                                                                   G4double time) const
{
  G4bool sameIndex = !SearchIndexUpdated(search, index);
  return SearchUpperBoundTime(search, time, sameIndex);
}

//------------------------------------------------------------------------------

template<typename TIndex>
std::vector<G4int> G4VUserMoleculeReactionCounter<TIndex>::GetNbReactionsAtTimes(
  const TIndex& index, const std::vector<G4double>& times) const
{
  Search search = {};
  std::vector<G4int> counts = {};
  for (auto time : times)
    counts.push_back(GetNbReactionsAtTime(search, index, time));
  return counts;
}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeReactionCounter<TIndex>::RecordReaction(
  std::unique_ptr<G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex> pIndex,
  G4double time, G4int number)
{
  const TIndex* mapIndex = dynamic_cast<TIndex*>(pIndex.get());

  if (IsTimeAboveUpperBound(time)) {
    if (fVerbose > 3) {
      G4cout << "G4VUserMoleculeReactionCounter<"
             << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">(" << GetName()
             << ")::RecordReaction : " << mapIndex->GetReactionData()->GetReactionID()
             << " at time : " << G4BestUnit(time, "Time") << G4endl;
      G4cout << ":: [IsTimeAboveUpperBound] Skipping since IsTimeAboveUpperBound == true" << G4endl;
    }
    return;
  }
  else if (IsTimeBelowLowerBound(time)) {
    if (fVerbose > 3) {
      G4cout << "G4VUserMoleculeReactionCounter<"
             << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">(" << GetName()
             << ")::RecordReaction : " << mapIndex->GetReactionData()->GetReactionID()
             << " at time : " << G4BestUnit(time, "Time") << G4endl;
      G4cout << ":: [IsTimeBelowLowerBound] Skipping since IsTimeBelowLowerBound == true" << G4endl;
    }
    return;
  }

  if (fVerbose > 2) {
    G4cout << "G4VUserMoleculeReactionCounter<"
           << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">(" << GetName()
           << ")::RecordReaction : " << mapIndex->GetReactionData()->GetReactionID()
           << " at time : " << G4BestUnit(time, "Time") << G4endl;
  }

  auto [it, indexIsNew] = fCounterMap.emplace(*mapIndex, InnerCounterMapType{fTimeComparer});

  if (indexIsNew || it->second.empty()) {
    it->second.emplace(time, number);
    // it->second[time] = number;
  }
  else {
    if (G4MoleculeCounterManager::Instance()->GetResetCountersBeforeEvent())
    // can only do consistency check if the counters are cleared before each event
    {
      auto end = it->second.rbegin();

      if ((end->first <= time || std::fabs(end->first - time) <= fTimeComparer.GetPrecisionAtTime(time)))
      // Case 1 = new time comes after last recorded data
      // Case 2 = new time is about the same as the last recorded one
      {
        // it->second[time] = end->second + number;
        auto [it_time, _] = it->second.emplace(time, end->second);
        it_time->second += number;
      }
      else {
        G4ExceptionDescription errMsg;
        errMsg << "Time of reaction " << mapIndex->GetReactionData()->GetReactionID() << " is "
               << G4BestUnit(time, "Time") << "while the global time is "
               << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
               << "(last counter time: " << G4BestUnit(end->first, "Time") << ")" << G4endl;
        G4Exception(G4String("G4VUserMoleculeReactionCounter<"
                             + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
                             + ">::RecordReaction"),
                    "TIME_DONT_MATCH", FatalException, errMsg);
      }
    }
    else {
      // since counters are not cleared after chemical run (i.e., after event)
      // there will already be numbers in the map, so...
      // (1) find the closest time
      // (2) emplace entry using closest value as init + number
      // (3) add number to all "future" entries as well
      auto it_closest = G4::MoleculeCounter::FindClosestEntryForKey(it->second, time);
      auto [it_new, _] = it->second.emplace(time, it_closest->second);
      do {
        it_new->second += number;
      } while (++it_new != it->second.end());
    }
  }
}

//------------------------------------------------------------------------------

template<typename TIndex>
std::vector<TIndex> G4VUserMoleculeReactionCounter<TIndex>::GetMapIndices() const
{
  if (fVerbose > 2) {
    G4cout << "Entering in G4VUserMoleculeReactionCounter::GetMapIndices" << G4endl;
  }
  return G4::MoleculeCounter::GetMapIndices(fCounterMap);
}

//------------------------------------------------------------------------------

template<typename T>
std::set<const G4DNAMolecularReactionData*> G4VUserMoleculeReactionCounter<T>::GetRecordedReactions() const
{
  if (fVerbose > 2) {
    G4cout << "Entering in G4VUserMoleculeReactionCounter::GetRecordedReactions" << G4endl;
  }
  std::set<const G4DNAMolecularReactionData*> output{};
  for (const auto& it : fCounterMap) {
    output.insert(it.first.GetReactionData());
  }
  return output;
}

//------------------------------------------------------------------------------

template<typename T>
std::set<G4double> G4VUserMoleculeReactionCounter<T>::GetRecordedTimes() const
{
  return G4::MoleculeCounter::GetRecordedTimes<T>(fCounterMap);
}

//------------------------------------------------------------------------------

template<typename T>
void G4VUserMoleculeReactionCounter<T>::Dump() const
{
  DumpCounterMapIndices();
  G4::MoleculeCounter::DumpCounterMapContents<T>(fCounterMap);
}

template<typename T>
void G4VUserMoleculeReactionCounter<T>::DumpCounterMapIndices() const
{
  G4::MoleculeCounter::DumpCounterMapIndices<T>(fCounterMap);
}

//------------------------------------------------------------------------------

template<typename T>
void G4VUserMoleculeReactionCounter<T>::ResetCounter()
{
  if (fVerbose > 1) {
    G4cout << "G4VUserMoleculeReactionCounter<" << G4::MoleculeCounter::GetTemplateTypeName<T>()
           << ">(" << GetName() << ")::ResetCounter" << G4endl;
  }
  fCounterMap.clear();
}

//------------------------------------------------------------------------------

template<typename TIndex>
G4bool G4VUserMoleculeReactionCounter<TIndex>::SearchIndexUpdated(Search& search, const TIndex& index) const
{
  if (search.fLowerBoundSet && !(search.fLastIndexSearched->first < index)
      && !(index < search.fLastIndexSearched->first))
  {
    return true;
  }

  auto mol_it = fCounterMap.find(index);
  search.fLastIndexSearched = mol_it;

  if (mol_it != fCounterMap.end()) {
    search.fLowerBoundTime = search.fLastIndexSearched->second.end();
    search.fLowerBoundSet = true;
  }
  else {
    search.fLowerBoundSet = false;
  }

  return false;
}

//------------------------------------------------------------------------------

template<typename T>
G4int G4VUserMoleculeReactionCounter<T>::SearchUpperBoundTime(Search& search, G4double time, G4bool sameIndex) const
{
  auto mol_it = search.fLastIndexSearched;
  if (mol_it == fCounterMap.end()) {
    return 0;
  }

  InnerCounterMapType const& timeMap = mol_it->second;
  if (timeMap.empty()) {
    return 0;
  }

  if (sameIndex) {
    if (search.fLowerBoundSet && search.fLowerBoundTime != timeMap.end()) {
      if (search.fLowerBoundTime->first < time) {
        auto upperToLast = search.fLowerBoundTime;
        upperToLast++;

        if (upperToLast == timeMap.end()) {
          return search.fLowerBoundTime->second;
        }

        if (upperToLast->first > time) {
          return search.fLowerBoundTime->second;
        }
      }
    }
  }

  auto up_time_it = timeMap.upper_bound(time);

  if (up_time_it == timeMap.end()) {
    auto last_time = timeMap.rbegin();
    return last_time->second;
  }
  if (up_time_it == timeMap.begin()) {
    return 0;
  }

  up_time_it--;

  search.fLowerBoundTime = up_time_it;
  search.fLowerBoundSet = true;

  return search.fLowerBoundTime->second;
}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeReactionCounter<TIndex>::AbsorbCounter(
  const G4VMoleculeCounterInternalBase* pCounterBase)
{
  if (pCounterBase == nullptr) {
    G4ExceptionDescription errMsg;
    errMsg << "Could not cast the pointer to type G4VUserMoleculeReactionCounter<"
           << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">!\n"
           << "Because the pointer is nullptr!" << G4endl;
    G4Exception(G4String("G4VUserMoleculeReactionCounter<"
                         + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AbsorbCounter"),
                "BAD_REFERENCE", FatalException, errMsg);
  }

  auto pCounter = dynamic_cast<G4VUserMoleculeReactionCounter<TIndex> const*>(pCounterBase);

  if (pCounter == nullptr) {
    G4ExceptionDescription errMsg;
    errMsg << "Could not cast the pointer to type G4VUserMoleculeReactionCounter<"
           << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">!\n"
           << "Because the objects aren't of the same type!" << G4endl;
    G4Exception(G4String("G4VUserMoleculeReactionCounter<"
                         + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AbsorbCounter"),
                "BAD_REFERENCE", FatalException, errMsg);
  }

  if (pCounter->GetType() != GetType()) {
    G4ExceptionDescription errMsg;
    errMsg << "You are trying to absorb a counter with different type!" << G4endl;
    G4Exception(G4String("G4VUserMoleculeReactionCounter<"
                         + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AbsorbCounter"),
                "TYPE_DIFF", JustWarning, errMsg);
  }

  for (auto const& worker_it : pCounter->GetCounterMap()) {
    auto [master_it, indexIsNew] =
      fCounterMap.emplace(worker_it.first, InnerCounterMapType{fTimeComparer});

    G4int currentNumber = 0, previousNumber = 0;
    for (auto const& [time, number] : worker_it.second) {
      currentNumber = number - previousNumber;
      previousNumber = number;

      if (master_it->second.empty()) {
        master_it->second.emplace(time, currentNumber);
      }
      else {  // at least one element exists, so we can try to find the closest key
        auto it_closest = G4::MoleculeCounter::FindClosestEntryForKey(master_it->second, time);
        auto [it, _] = master_it->second.emplace(time, it_closest->second);
        do {
          it->second += currentNumber;
        } while (++it != master_it->second.end());
      }
    }
  }
}

//------------------------------------------------------------------------------

#endif
