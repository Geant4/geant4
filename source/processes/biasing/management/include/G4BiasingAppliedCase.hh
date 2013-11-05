#ifndef G4BiasingAppliedCase_hh
#define G4BiasingAppliedCase_hh

enum G4BiasingAppliedCase
 {
   BAC_None,               // -- not under biasing
   BAC_NonPhysics,         // -- splitting, killing (not a physics process biasing)
   BAC_DenyInteraction,    // -- physics process occurence biasing, denying production of physics process final state
   BAC_FinalState,         // -- physics process final state biasing only
   BAC_Occurence           // -- physics process occurence biasing; may come together with a final state biasing
 };

#endif
