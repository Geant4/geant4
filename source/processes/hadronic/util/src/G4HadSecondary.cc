#include "G4HadSecondary.hh"

G4HadSecondary::G4HadSecondary(G4DynamicParticle * aT, G4double aWeight) :
    theP(aT), theWeight(aWeight), theTime(-1)
{
  if(aT->GetKineticEnergy()<0)
  {
    G4cout << "!!!!! Got it all right !!!!! "<<G4endl;
    G4Exception("");
  }
}
