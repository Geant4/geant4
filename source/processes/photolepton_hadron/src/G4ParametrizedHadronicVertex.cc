#include "G4ParametrizedHadronicVertex.hh"

G4VParticleChange * G4ParametrizedHadronicVertex::
ApplyYourself(const G4Nucleus & theTarget, const G4Track &thePhoton)
{   
    G4double theKineticEnergy = thePhoton.GetKineticEnergy();
    if(RandFlat::shootBit())
    {
      if(theKineticEnergy<20*GeV) return theLowEPionMinus.ApplyYourself(thePhoton, theTarget);
      return theHighEPionMinus.ApplyYourself(thePhoton, theTarget);
    }
    else
    {
      if(theKineticEnergy<20*GeV) return theLowEPionPlus.ApplyYourself(thePhoton, theTarget);
      return theHighEPionPlus.ApplyYourself(thePhoton, theTarget);
    }
    return NULL;
}
