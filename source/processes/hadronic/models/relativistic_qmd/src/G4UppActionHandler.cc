
#include "G4UppActionHandler.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4UppActionCollision.hh"


bool CmpAction::operator() (const G4VUppAction* a, const G4VUppAction* b) const
{ 
  // G4cout << "less: " << a->getActionTime() << " vs. " << b->getActionTime() << G4endl;
  return (a->getActionTime() < b->getActionTime()); 
}


void G4UppActionHandler::CleanUp()
{
  G4cout << "cleanup: " << q.size() << " elements" << endl;
  G4int c=0;
  for (queuetype::iterator i=q.begin(); i!=q.end(); i++) {
    (*i)->dump();
    if (!((*i)->isValid())) {
      cout << "remove: ";
      (*i)->dump();
      delete *i;
      q.erase(i--);
      c++;
    }
  }
  G4cout << "(debug) cleanup: " << c << "   " << q.size() << " left." << endl;
}

  

void G4UppActionHandler::UpdateActions(const G4UppTrackVector& v)
{
  const G4double CutOff = 8.0*fermi;
  G4int c=0;
  cout << "(debug) update..." << endl;
  for (int i=0; i<v.size()-1; i++) {
    for (int j=i+1; j<v.size(); j++) {
      // cout << i << "  " << j << endl;
      if ( (v[i]->hasChanged()) || (v[j]->hasChanged())) {
	c++;
	// cout << c << endl;
	if ( (i!=j) 
	     && (!lastPartners(*v[i],*v[j])) 
	     && (!sameGroup(*v[i],*v[j])) ) {
	  G4double minDist = MinimumDistance(*v[i],*v[j]);
	  if ( (minDist > 0.0) && (minDist < CutOff) ) {
	    G4double collTime = CollisionTime(*v[i],*v[j]);
	    if (collTime>0.0) {
	      G4UppTrackVector inParts;
	      inParts.push_back(v[i]);
	      inParts.push_back(v[j]);
	      q.push_back(new G4UppActionCollision(collTime + v.getGlobalTime(),inParts));
	      cout << "update!!" << i << " hits " << j << endl;
	    }
	  }
	}
      }
    }
  }
  // cout << "vor sort" << endl;
  // cout << "actions: " << q.size() << endl;
  sort(q.begin(),q.end(),CmpAction());
  cout << "(debug) " << c << "pairs checked" << endl;
}


G4bool G4UppActionHandler::lastPartners(const G4UppTrack& i, const G4UppTrack& j) const
{
  if (i.isLastPartner(&j)) return true;
  if (j.isLastPartner(&i)) return true;
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


void G4UppActionHandler::dump(const G4UppTrackVector& t) const
{
  G4cout << "actionhandler: " << q.size() << " elements" << endl;
  for (queuetype::iterator i=q.begin(); i!=q.end(); i++) {
    (*i)->dump(t);
  }
}


const G4VUppAction& G4UppActionHandler::getFirstAction() const
{
  cout << "getAction..." << endl;
  return *q.front();
}


void G4UppActionHandler::deleteFirstAction()
{
  delete q.front();
  q.erase(q.begin());
}


G4double G4UppActionHandler::MinimumDistance(const G4UppTrack& i, 
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


G4double G4UppActionHandler::CollisionTime(const G4UppTrack& i, 
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

