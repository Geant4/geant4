

#include "G4UppSimpleFieldtransport.hh"
#include "G4LorentzVector.hh"
#include "G4UppGlobal.hh"


void G4UppSimpleFieldtransport::Propagate(G4UppTrackVector& a, const G4double dTime) const
{
  if (uppdebug) cout << "(debug) fieldtransport dt=" << dTime/fermi << endl;
  if (uppdebug) {
    cout << "(debug) from " << a.getGlobalTime()/fermi;
    cout << " to " << (a.getGlobalTime()+dTime)/fermi << endl;
  }
  for (G4int i=0; i<a.size(); i++) {
    G4LorentzVector aMomentum = a[i]->Get4Momentum();
    G4double E = aMomentum.e();
    // sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)    
    G4LorentzVector aPosition = a[i]->Get4Position();
    aPosition.setT(aPosition.t() + dTime);
    aPosition.setVect(aPosition.vect() + aMomentum.vect() * (dTime / E));
    a[i]->Set4Momentum(aMomentum);
    a[i]->Set4Position(aPosition);
  }
  a.setGlobalTime(a.getGlobalTime() + dTime);
} 
