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
#ifndef G4DNAEventSet_hh
#define G4DNAEventSet_hh 1
#include <list>
#include <map>
#include <G4memory.hh>
#include "G4Track.hh"
#include <set>
#include "G4DNAMesh.hh"
#include <variant>
class G4DNAMolecularReactionTable;
class G4DNAMolecularReactionData;
class Event
{
 public:
  using Index        = G4Voxel::Index;
  using MolType      = const G4MolecularConfiguration*;
  using JumpingData  = std::pair<MolType, Index>;
  using ReactionData = const G4DNAMolecularReactionData;

  // to test C++17
  // using Data = std::variant<std::unique_ptr<JumpingData>, ReactionData*>
  using Data = std::pair<std::unique_ptr<JumpingData>, ReactionData*>;

  Event(G4double time, unsigned int key, ReactionData*);
  Event(G4double time, unsigned int key, std::unique_ptr<JumpingData>&&);

  virtual ~Event();
  G4double GetTime() const { return fTimeStep; }
  unsigned int GetKey() const { return fKey; }
  void PrintEvent() const;

  JumpingData* GetJumpingData() const { return std::get<0>(fData).get(); }

  ReactionData* GetReactionData() const { return std::get<1>(fData); }

 private:
  G4double fTimeStep;
  unsigned int fKey;
  Data fData;
};

struct comparatorEventSet
{
  G4bool operator()(std::unique_ptr<Event> const& rhs,
                    std::unique_ptr<Event> const& lhs) const;
};

class IEventSet
{
 public:
  IEventSet()  = default;
  ~IEventSet() = default;
};

class G4DNAEventSet : public IEventSet
{
 public:
  using Key      = unsigned int;
  using EventSet = std::set<std::unique_ptr<Event>, comparatorEventSet>;
  using EventMap = std::map<Key, EventSet::iterator>;
  G4DNAEventSet();
  virtual ~G4DNAEventSet();

  void CreateEvent(G4double time, Key index,
                   Event::ReactionData* pReactionData);
  void CreateEvent(G4double time, Key index,
                   std::unique_ptr<Event::JumpingData> jum);

  void AddEvent(std::unique_ptr<Event> pEvent);
  void RemoveEventSet()
  {
    fEventSet.clear();
    fEventMap.clear();
  }
  void RemoveEventOfVoxel(const size_t& key);

  EventSet::iterator end() { return fEventSet.end(); }
  EventSet::iterator begin() { return fEventSet.begin(); }

  EventSet::reverse_iterator rend() { return fEventSet.rend(); }
  EventSet::reverse_iterator rbegin() { return fEventSet.rbegin(); }

  EventSet::const_iterator end() const { return fEventSet.end(); }

  EventSet::const_iterator begin() const { return fEventSet.begin(); }

  size_t size() { return fEventSet.size(); }

  G4bool Empty() { return fEventSet.empty(); }

  void RemoveEvent(EventSet::iterator iter);
  [[maybe_unused]] void PrintEventSet();

 private:
  EventSet fEventSet;
  EventMap fEventMap;
};
#endif
