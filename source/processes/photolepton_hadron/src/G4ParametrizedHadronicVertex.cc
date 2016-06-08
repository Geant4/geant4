// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParametrizedHadronicVertex.cc,v 1.3 2000/08/03 09:06:33 gcosmo Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// --------------------------------------------------------------
#include "G4ParametrizedHadronicVertex.hh"

G4VParticleChange * G4ParametrizedHadronicVertex::
ApplyYourself(G4Nucleus & theTarget, const G4Track &thePhoton)
{   
    G4double theKineticEnergy = thePhoton.GetKineticEnergy();
    if(RandBit::shootBit())
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
