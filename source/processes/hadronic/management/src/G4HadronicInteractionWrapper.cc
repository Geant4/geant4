#include "G4HadronicInteractionWrapper.hh"

  G4HadFinalState * G4HadronicInteractionWrapper::
  ApplyInteraction(G4HadProjectile & thePro, 
                   G4Nucleus & targetNucleus,
                   G4HadronicInteraction * theInteraction)
  {
    static G4HadronicWhiteBoard & theBoard = G4HadronicWhiteBoard::Instance();
    theBoard.SetProjectile(thePro);
    theBoard.SetTargetNucleus(targetNucleus);
    return theInteraction->ApplyYourself( thePro, targetNucleus);
  }
