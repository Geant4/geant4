
#include "G4UppTrackVector.hh"
#include "G4UppTrackChange.hh"
#include "G4UppGlobal.hh"


void G4UppTrackVector::add(const G4KineticTrackVector& v)
{
  for (G4int i=0; i<v.length(); i++) {
    this->push_back(new G4UppTrack(*v[i]));
  }
}

    

void G4UppTrackVector::add(const G4RWTPtrOrderedVector<G4Nucleon>& v, const G4int g)
{
  for (G4int i=0; i<v.length(); i++) {
    G4UppTrack* newTrack = new G4UppTrack(*v[i],g);
    // cout << "got nuc: " << newTrack->GetDefinition()->GetParticleName() << ":";
    // cout << newTrack->Get4Momentum().e() << ":";
    // cout << newTrack->GetActualMass()<< endl;
    this->push_back(newTrack);
  }
  if (uppdebug) cout << "(debug) " << v.length() << " Nucleons added." << endl;
  if (uppdebug) cout << "(debug) nonCollisionGroup: " << g << endl;
}


G4int G4UppTrackVector::getIndex(const G4UppTrack* tPtr) const
{
  for (G4int i=0; i<this->size(); i++) {
    if ((*this)[i] == tPtr) return i;
  }
  return -1;
}


void G4UppTrackVector::dump() const
{
  cout << "dumping..." << this->size() << endl;
  for (G4int i=0; i<this->size(); i++) {
    G4cout << i << ": " << endl;
    (*this)[i]->dump();
  }
}


G4UppTrackChange G4UppTrackVector::update(const G4UppTrackChange& aChange)
{
  G4UppTrackChange backup;
  
  // aChange.dump();

  for (G4int i=0; i<aChange.oldParticles.size(); i++) {
    backup.newParticles.push_back(aChange.oldParticles[i]);
    G4int Index = getIndex(aChange.oldParticles[i]);
    // cout << "Erase " << Index << endl;
    this->erase(this->begin()+Index);
  }  

  for (G4int j=0; j<aChange.newParticles.size(); j++) {
    cout << "(info) adding " << aChange.newParticles[j]->GetDefinition()->GetParticleName() << endl;
    this->push_back(aChange.newParticles[j]);
    this->back()->setChanged(true);
    backup.oldParticles.push_back(this->back());
  }
  return backup;
}


G4bool G4UppTrackVector::isPauliBlocked(const G4UppTrackVector& t) const
{
  for (G4int i=0; i<t.size(); i++) {
    if (isPauliBlocked(t[i]))
      return true;
  }
  return false;
}


void G4UppTrackVector::resetChangedFlags()
{
  for (G4int i=0; i<this->size(); i++) {
    (*this)[i]->setChanged(false);
  }
}
    

G4bool G4UppTrackVector::isPauliBlocked(const G4UppTrack* aTrack) const
{ 
  const G4double afit=1.49641; 
  const G4double bfit=0.208736;
  const G4double gw=0.25;

  G4bool isBlocked=false;
  G4double r2 = 0.0;
  G4double p2;
  G4double rho = 0.0;
  G4double rhob = 0.0;
  G4double rhob0 = 0.0;
  G4double f = 0.0;
  G4double pgw = 1.0/hbarc/hbarc/gw;
	 

  if (aTrack->GetDefinition()->GetBaryonNumber() == 1) {
    G4int i = getIndex(aTrack);
    for (G4int j=0; j<size(); j++) {
      //    r2 = (rx(i)-rx(j))**2+(ry(i)-ry(j))**2+(rz(i)-rz(j))**2;
      r2 = (aTrack->GetPosition() - (*this)[j]->GetPosition() ).mag2();
      if ((aTrack->GetDefinition()==(*this)[j]->GetDefinition())) {
	p2 = (G4ThreeVector(aTrack->Get4Momentum()) - 
	      G4ThreeVector((*this)[j]->Get4Momentum()) ).mag2();
	//      p2 = pow(px(i)-px(j),2)+pow(py(i)-py(j),2)
	//           +pow(pz(i)-pz(j),2);
	p2 *= 0.25;
	rho += rho + exp(-(2.0*gw*r2));
	f += f + exp(-(gw*r2)-pgw*p2);
      }
      rhob += rhob + exp(-(2.0*gw*r2));
    }
  }
  
  G4double test = afit + bfit*rho;
  if (test > f) 
    isBlocked=true;

  //  rhob = rhob*(2.0*gw/pi)**1.5/rho0;
  // Baryonendichte in Einheiten von rho0
  G4cout << "paulib: f " << f;
  G4cout << " rho " << rho;
  G4cout << " rhob " << rhob;
  G4cout << " blocking " << isBlocked << G4endl;
  
  return isBlocked;
}

