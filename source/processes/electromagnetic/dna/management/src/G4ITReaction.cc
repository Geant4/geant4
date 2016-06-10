/*
 * G4ITReactionInfo.cc
 *
 *  Created on: 1 févr. 2015
 *      Author: matkara
 */

#include <G4ITReaction.hh>
#include "globals.hh"

//G4ThreadLocal std::set<G4ITReaction*>* G4ITReaction::gAll(0);


G4ITReaction::G4ITReaction(double time, G4Track* trackA, G4Track* trackB) :
  G4enable_shared_from_this<G4ITReaction>(),
  fTime(time),
  fReactants(trackA,trackB)
{
  //if(gAll == 0) gAll = new std::set<G4ITReaction*>();
  //gAll->insert(this);
}

G4ITReaction::~G4ITReaction()
{
  //gAll->erase(this);
}

void G4ITReaction::RemoveMe(G4ITReactionSet* reactionSet)
{
  G4ITReactionPtr backMeUp = this->shared_from_this();
  for(G4ReactionPerTrackIt::iterator it = fReactionPerTrack.begin() ;
      it != fReactionPerTrack.end() ; ++it)
  {
    // G4cout << it->first.get() << G4endl;
    // assert(it->first.get() != 0);
    it->first->RemoveReaction(it->second, reactionSet);
  }
  fReactionPerTrack.clear();
}

bool G4ITReactionPerTrack::RemoveReaction(G4ITReactionList::iterator it,
                                          G4ITReactionSet* reactionSet)
{
  // G4cout << "G4ITReactionPerTrack::RemoveReaction" << G4endl;
  fReactions.erase(it);
  if(fReactions.empty())
  {
    // G4cout << "G4ITReactionPerTrack is empty" << G4endl;
    reactionSet->RemoveReactionPerTrack(this->shared_from_this());
    return true;
  }
  return false;
}
