//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4CollisionManager.hh"
#include "G4CollisionInitialState.hh"

G4CollisionManager::G4CollisionManager()
{
  theCollisionList = new G4ListOfCollisions;
}


G4CollisionManager::~G4CollisionManager()
{
  ClearAndDestroy();
  delete theCollisionList;
}


void G4CollisionManager::AddCollision(G4double time, G4KineticTrack * proj,
				      G4KineticTrack * target)

{
      if ( time < DBL_MAX )
      {
         G4CollisionInitialState *
	     collision = new G4CollisionInitialState(time, proj, target);
         theCollisionList->push_back(collision);
      } else {
         G4cerr << "G4Scatterer invalid TimeTo Interaction : " << time;
         G4cerr <<"    projectile "<<proj->Get4Momentum()<<" "
                                   <<proj->GetDefinition()->GetParticleName()<<G4endl;
         if (target) G4cerr <<"    target     "
				<<target->Get4Momentum()<<" "
                                <<target->GetDefinition()->GetParticleName()<<G4endl;
         G4cerr <<"G4Scatterer error message end"<< G4endl;
	 G4Exception("G4Scatterer::AddCollision()");
      }
}


void G4CollisionManager::RemoveCollision(G4CollisionInitialState * collision)
{
  theCollisionList->erase(std::find(theCollisionList->begin(),
				      theCollisionList->end(),
				      collision));
  delete collision;
  collision = NULL; // prevent future use of the pointer
}


void G4CollisionManager::RemoveTracksCollisions(G4KineticTrackVector * ktv)
{
  if(ktv == NULL)
    return;
  if(ktv->empty())
    return;

  G4CollisionInitialState * collision;
  std::vector<G4CollisionInitialState *>::iterator collIter, collIter2;
  std::vector<G4KineticTrack *>::iterator trackIter;
  G4ListOfCollisions toRemove;

  for(collIter = theCollisionList->begin();
      collIter != theCollisionList->end(); ++collIter)
  {
    collision = *collIter;
    for(trackIter = ktv->begin(); trackIter != ktv->end(); ++trackIter)
    {
      if((collision->GetTarget() == *trackIter) ||
	 (collision->GetPrimary() == *trackIter))
      {  // cannot remove the collision from the list inside the loop. Save and do it later
	toRemove.push_back(collision);
	break;  // exit from the "trackIter" loop
      }
    }
  }

  // now remove the collisions
  for(collIter = toRemove.begin(); collIter != toRemove.end(); ++collIter)
  {
    collision = *collIter;
    collIter2 = std::find(theCollisionList->begin(),
			    theCollisionList->end(), collision);
    theCollisionList->erase(collIter2);  // remove from list...
    delete collision;                    // ...and delete the collision
  }
}


void G4CollisionManager::ClearAndDestroy()
{
  std::vector<G4CollisionInitialState *>::iterator i;
  for(i = theCollisionList->begin(); i != theCollisionList->end(); ++i)
    delete *i;
  theCollisionList->clear();
}


G4CollisionInitialState * G4CollisionManager::GetNextCollision()
{
  G4CollisionInitialState * theNext=0;
  G4double nextTime = DBL_MAX;
  std::vector<G4CollisionInitialState *>::iterator i;
  for(i = theCollisionList->begin(); i != theCollisionList->end(); ++i)
  {
    if(nextTime > (*i)->GetCollisionTime())
    {
      nextTime = (*i)->GetCollisionTime();
      theNext = *i;
    }
  }
  #ifdef debug_G4CollisionManager
  if(theNext == 0 && theCollisionList->size()!=0)
  {
    G4double debugTime = DBL_MAX;
    G4cerr <<"G4CollisionManager::GetNextCollision - Fatal"<<G4endl;
    G4cerr <<" number of collisions left "<<theCollisionList->size()<<G4endl;
    for(i = theCollisionList->begin(); i != theCollisionList->end(); ++i)
    {
      G4cerr <<" Time to collision "<<(*i)->GetCollisionTime()<<" "<<G4endl;
      G4cerr <<"    projectile "<<(*i)->GetPrimary()->Get4Momentum()<<" "
                                <<(*i)->GetPrimary()->GetDefinition()->GetParticleName()<<G4endl;
      if ((*i)->GetTarget()) G4cerr <<"    target     "<<(*i)->GetTarget()->Get4Momentum()<<" "
                                <<(*i)->GetTarget()->GetDefinition()->GetParticleName()<<G4endl;
    }
    G4cerr <<"G4CollisionManager::GetNextCollision - End of message"<<G4endl;
  }
  #endif

  return theNext;
}


void G4CollisionManager::Print()
{
  std::vector<G4CollisionInitialState *>::iterator i;

  G4cout << "CollisionManager: " << theCollisionList->size()
	 << " entries at " << theCollisionList << G4endl;
  G4CollisionInitialState * collision;
  for(i = theCollisionList->begin(); i != theCollisionList->end(); ++i)
  {
    collision = *i;
    G4cout << "  collision " << collision << " time: "
	   << collision->GetCollisionTime()/second << " proj: "
	   << collision->GetPrimary() << " trgt: "
	   << collision->GetTarget() << G4endl;
  }
}
