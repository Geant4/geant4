#include "G4HadSecondary.hh"
#include "G4HadronicException.hh"

G4HadSecondary::G4HadSecondary(G4DynamicParticle * aT, G4double aWeight) :
    theP(aT), theWeight(aWeight), theTime(-1)
{
  if(aT->GetKineticEnergy()<0)
  {
    throw(G4HadronicException(__FILE__, __LINE__, 
    "ATTEMPTING TO CREATE A SECONDARY WITH NGATIVE KINETIC ENERGY.") );
  }
}
