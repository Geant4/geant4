
#include "G4UppActionHandler.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4UppActionCollision.hh"


G4UppActionHandler::G4UppActionHandler(const G4UppTrackVector& allTracks)
{
  allTracksPtr = &allTracks;
}


bool G4UppActionHandler::CmpAction::operator() (const G4VUppAction* a, const G4VUppAction* b) const
{ 
  // G4cout << "less: " << a->getActionTime() << " vs. " << b->getActionTime() << G4endl;
  return (a->getActionTime() < b->getActionTime()); 
}


void G4UppActionHandler::cleanUp()
{
  // G4cout << "cleanup: " << q.size() << " elements" << endl;
  G4int count=0;
  for (queuetype::iterator i=q.begin(); i!=q.end(); i++) {
    // (*i)->dump();
    if (!( (*i)->isValid() )) {
      G4cout << "(debug) removing: ";
      (*i)->dump();
      delete *i;
      q.erase(i--);
      count++;
    }
  }
  G4cout << "(debug) cleanup: " << count << " actions removed" << endl;
}

  

void G4UppActionHandler::updateActions()
{
  const G4double cutOff = 8.0*fermi;
  G4int count=0;
  G4cout << "(debug) update..." << G4endl;
  for (int i=0; i<allTracksPtr->size()-1; i++) {
    for (int j=i+1; j<allTracksPtr->size(); j++) {
      G4UppTrack* trk1 = (*allTracksPtr)[i];
      G4UppTrack* trk2 = (*allTracksPtr)[j];
      // cout << i << "  " << j << endl;
      if ( (trk1->hasChanged()) || (trk2->hasChanged()) ) {
	if ( (i!=j) 
	     && ( !lastPartners(*trk1,*trk2) ) 
	     && ( !sameGroup(*trk1,*trk2) ) ) {
	  G4double minDist = minimumDistance(*trk1,*trk2);
	  if ( (minDist > 0.0) && (minDist < cutOff) ) {
	    G4double collTime = collisionTime(*trk1,*trk2);
	    if (collTime>0.0) {
	      G4UppTrackVector collidingParticles;
	      collidingParticles.push_back(trk1);
	      collidingParticles.push_back(trk2);
	      G4VUppAction* actionPtr = new G4UppActionCollision(collTime + allTracksPtr->getGlobalTime(),
								 *allTracksPtr,
								 collidingParticles,
								 aScatterer);
	      G4cout << "(debug)  new ";
	      actionPtr->dump();
	      q.push_back(actionPtr);
	      count++;
	    }
	  }
	}
      }
    }
  }
  sort(q.begin(),q.end(),CmpAction());
  cout << "(debug) " << count << " collision(s) updated" << endl;
}


G4bool G4UppActionHandler::lastPartners(const G4UppTrack& i, const G4UppTrack& j) const
{
  if (i.isLastInteractionPartner(&j)) return true;
  if (j.isLastInteractionPartner(&i)) return true;
  return false;
}


G4bool G4UppActionHandler::sameGroup(const G4UppTrack& i, const G4UppTrack& j) const
{
  if ( (i.getNumberOfCollisions()==0) && (j.getNumberOfCollisions()==0)) {
    if (i.getNonInteractionGroup()==j.getNonInteractionGroup()) {
      return true;
    }
  }
  return false;
}


void G4UppActionHandler::addAction(G4VUppAction* a)
{
  q.push_back(a);
  sort(q.begin(),q.end(),CmpAction());
}


void G4UppActionHandler::dump() const
{
  G4cout << "actionhandler with " << q.size() << " elements" << G4endl;
  for (queuetype::iterator i=q.begin(); i!=q.end(); i++) {
    (*i)->dump();
  }
}


const G4VUppAction* G4UppActionHandler::getFirstAction() const
{
  return q.front();
}


void G4UppActionHandler::deleteFirstAction()
{
  delete q.front();
  q.erase(q.begin());
}


G4double G4UppActionHandler::minimumDistance(const G4UppTrack& i, 
					     const G4UppTrack& j) const
{
  G4LorentzVector pi = i.Get4Momentum();
  G4LorentzVector pj = j.Get4Momentum();
  G4LorentzVector qi = i.Get4Position();
  G4LorentzVector qj = j.Get4Position();
  // cout << "i " << qi << endl;
  // cout << "i " << pi << endl;
  // cout << "k " << qj << endl;
  // cout << "k " << pj << endl;
  G4ThreeVector Beta = - (pi.vect() + pj.vect()) * (1.0 / (pi.e() + pj.e()));
  G4LorentzVector pij = (pi-pj);
  pij.boost(Beta);
  G4LorentzVector qij = (qi-qj);
  qij.boost(Beta);
  // cout << "qij " << qij*qij << " qp " << qij*pij << " pp " << pij*pij << endl;
  G4double dist = ((qij*pij)*(qij*pij))/(pij*pij) - qij*qij;
  // cout << "distance  " << dist << endl;
  return dist;
}


G4double G4UppActionHandler::collisionTime(const G4UppTrack& i, 
					   const G4UppTrack& j) const
{
  G4LorentzVector pi = i.Get4Momentum();
  G4LorentzVector pj = j.Get4Momentum();
  G4LorentzVector qi = i.Get4Position();
  G4LorentzVector qj = j.Get4Position();
  G4ThreeVector beta_i = pi.vect() * (1.0/pi.e());
  G4ThreeVector beta_j = pj.vect() * (1.0/pj.e());
  G4double time = - (qi-qj)*(beta_i-beta_j)/((beta_i-beta_j)*(beta_i-beta_j));
  // cout << "time : " << time << endl;
  return time;
}

