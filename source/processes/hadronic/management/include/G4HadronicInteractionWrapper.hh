#ifndef G4HadronicInteractionWrapper_h
#define G4HadronicInteractionWrapper_h 1

#include "G4HadronicWhiteBoard.hh"
#include "G4HadFinalState.hh"
#include "G4HadronicInteraction.hh"

class G4HadronicInteractionWrapper
{
  public:
  G4HadFinalState * ApplyInteraction(G4HadProjectile & thePro, 
                                     G4Nucleus & targetNucleus,
                                     G4HadronicInteraction * theInteraction);
};

#endif
