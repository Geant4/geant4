#ifndef G4AtomicDeexcitation_h
#define G4AtomicDeexcitation_h 1
#include "globals.hh"
#include "g4std/vector"
#include "G4AtomicTransitionManager.hh"
#include "G4DynamicParticle.hh"


class G4AtomicDeexcitation {

public:

  G4AtomicDeexcitation();
  ~G4AtomicDeexcitation();
   G4std::vector<G4DynamicParticle*>* GenerateParticles(G4int Z, G4int shellId);
			
private:

  //static  G4AtomicTransitionManager* transitionManager;
 
  const G4int SelectTypeOfTransition(G4int Z, G4int shellId);
 
  G4DynamicParticle* GenerateFluorescence(G4int Z, G4int shellId,G4int provShellId);
 
  G4DynamicParticle* GenerateAuger(G4int Z, G4int shellId);
 
  G4int newShellId;
 
};

#endif
