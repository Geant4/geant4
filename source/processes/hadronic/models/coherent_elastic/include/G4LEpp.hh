 //
 // G4 Low energy model: nn or pp scattering
 // L.G. Greeniaus and F.W. Jones, TRIUMF, May-1999
 // modified by Joe Chuma, Jan 2000
 //  
 // For further comments see G4LEpp.cc
 //
#ifndef G4LEpp_h
#define G4LEpp_h 1
 
#include "globals.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4HadronicInteraction.hh"

class G4LEpp : public G4HadronicInteraction
{
 private:

   enum { NENERGY=22, NANGLE=180 };

 public:

   G4LEpp();

   ~G4LEpp();
 
   G4VParticleChange* ApplyYourself(const G4Track& aTrack,
                                    G4Nucleus& targetNucleus);
 private:
 //
 // The following arrays are declared static to allow the use of initializers.
 // They are initialized in G4LEpp.cc
 //
   static G4float sig[NENERGY][NANGLE];
   static G4float pcm[NENERGY], elab[NENERGY], 
     dsigmax[NENERGY], sigtot[NENERGY];
 };

#endif
