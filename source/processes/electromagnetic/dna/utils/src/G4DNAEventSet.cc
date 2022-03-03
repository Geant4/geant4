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
#include <memory>
#include "G4DNAEventSet.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4DNAMolecularReactionTable.hh"

Event::Event(G4double time, unsigned int key, ReactionData* pReactionData)
  : fTimeStep(time)
  , fKey(key)
  , fData(std::pair<std::unique_ptr<JumpingData>, ReactionData*>(nullptr,
                                                                 pReactionData))
{}

Event::Event(G4double time, unsigned int key,
             std::unique_ptr<JumpingData>&& jumping)
  : fTimeStep(time)
  , fKey(key)
  , fData(std::pair<std::unique_ptr<JumpingData>, ReactionData*>(
      std::move(jumping), nullptr))
{}

Event::~Event() = default;

void Event::PrintEvent() const
{
  G4cout << "****PrintEvent::TimeStep : " << G4BestUnit(fTimeStep, "Time")
         << " key : " << fKey << " action : ";
  if(std::get<0>(fData) == nullptr)
  {
    G4cout << std::get<1>(fData)->GetReactant1()->GetName() << " + "
           << std::get<1>(fData)->GetReactant2()->GetName() << " -> "
           << std::get<1>(fData)->GetProducts()->size() << G4endl;
  }
  else
  {
    G4cout << std::get<0>(fData)->first->GetName() << " jumping to "
           << std::get<0>(fData)->second << G4endl;
  }
}

G4bool comparatorEventSet::operator()(std::unique_ptr<Event> const& rhs,
                                      std::unique_ptr<Event> const& lhs) const
{
  return rhs->GetTime() < lhs->GetTime();
}

G4DNAEventSet::G4DNAEventSet()
  : fEventSet(comparatorEventSet())
{}

void G4DNAEventSet::CreateEvent(G4double time, Key key,
                                Event::ReactionData* pReactionData)
{
  auto pEvent = std::make_unique<Event>(time, key, pReactionData);
  AddEvent(std::move(pEvent));
}

void G4DNAEventSet::CreateEvent(G4double time, Key key,
                                std::unique_ptr<Event::JumpingData> jum)
{
  auto pEvent = std::make_unique<Event>(time, key, std::move(jum));
  AddEvent(std::move(pEvent));
}

void G4DNAEventSet::RemoveEventOfVoxel(const size_t& key)
{
  auto it = fEventMap.find(key);
  if(it != fEventMap.end())
  {
    fEventSet.erase(it->second);
    fEventMap.erase(it);
  }
}

void G4DNAEventSet::RemoveEvent(EventSet::iterator iter)
{
  auto key = (*iter)->GetKey();
  RemoveEventOfVoxel(key);
}

void G4DNAEventSet::AddEvent(std::unique_ptr<Event> pEvent)
{
  // idea is no 2 events in one key (or index)
  auto key = pEvent->GetKey();
  RemoveEventOfVoxel(key);
  auto it        = fEventSet.emplace(std::move(pEvent));
  fEventMap[key] = std::get<0>(it);
}

//_____________________________________________________________________________________

G4DNAEventSet::~G4DNAEventSet() { RemoveEventSet(); }

[[maybe_unused]] void G4DNAEventSet::PrintEventSet()
{
  G4cout << "G4DNAEventSet::PrintEventSet()" << G4endl;
  for(const auto& it : fEventSet)
  {
    (*it).PrintEvent();
  }
  G4cout << "End PrintEventSet()" << G4endl;
  G4cout << G4endl;
}