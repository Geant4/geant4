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
/*
 * G4ITReactionInfo.cc
 *
 *  Created on: 1 févr. 2015
 *      Author: matkara
 */

#include <G4ITReaction.hh>
#include "globals.hh"

//G4ThreadLocal std::set<G4ITReaction*>* G4ITReaction::gAll(0);

G4ThreadLocal G4ITReactionSet* G4ITReactionSet::fpInstance(nullptr);

G4bool compReactionPerTime::operator()(G4ITReactionPtr rhs,
                                     G4ITReactionPtr lhs) const
{
  G4double time1 = rhs->GetTime();
  G4double time2 = lhs->GetTime();
  if (time1 == time2)
  {
    return rhs->GetHash() < lhs->GetHash();
  }
  return rhs->GetTime() < lhs->GetTime();
}


template<class T>
inline void hash_combine(std::size_t& seed,
                         const T& v)
{
//  std::hash<T> hasher;
  seed ^= /*std::hash<T>(v)*/ v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

std::size_t G4ITReaction::GetHash() const
{
  std::size_t hash = 0;
  hash_combine(hash, fReactants.first->GetTrackID());
  hash_combine(hash, fReactants.first->GetTrackID());
  return hash;
}

G4ITReaction::G4ITReaction(G4double time, G4Track* trackA, G4Track* trackB) :
  G4enable_shared_from_this<G4ITReaction>(),
  fTime(time),
  fReactants(trackA,trackB)
{
  //if(gAll == 0) gAll = new std::set<G4ITReaction*>();
  //gAll->insert(this);
  fReactionPerTimeIt = 0;
}

G4ITReaction::~G4ITReaction()
{
  //gAll->erase(this);
  if(fReactionPerTimeIt) delete fReactionPerTimeIt;
}

void G4ITReaction::RemoveMe()
{
  G4ITReactionPtr backMeUp = this->shared_from_this();
  for(auto it = fReactionPerTrack.begin() ;
      it != fReactionPerTrack.end() ; ++it)
  {
    it->first->RemoveThisReaction(it->second);
  }
  fReactionPerTrack.clear();

  if(fReactionPerTimeIt)
  {
   G4ITReactionSet::Instance()->GetReactionsPerTime().erase(*fReactionPerTimeIt);
   delete fReactionPerTimeIt;
   fReactionPerTimeIt = nullptr;
  }
}

G4bool G4ITReactionPerTrack::RemoveThisReaction(G4ITReactionList::iterator it)
{
  fReactions.erase(it);
  if(fReactions.empty())
  {
    G4ITReactionSet::Instance()->RemoveReactionPerTrack(this->shared_from_this());
    return true;
  }
  return false;
}
