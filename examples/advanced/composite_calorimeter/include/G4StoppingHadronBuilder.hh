#ifndef G4StoppingHadronBuilder_h
#define G4StoppingHadronBuilder_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"


// At rest processes
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4PionMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorption.hh"

class G4StoppingHadronBuilder 
{
  public: 
    G4StoppingHadronBuilder();
    virtual ~G4StoppingHadronBuilder();

  public: 
    virtual void Build();

  private:

   // Pi -
   G4PionMinusAbsorptionAtRest thePionMinusAbsorption;

   // K -
   G4KaonMinusAbsorption theKaonMinusAbsorption;
 
   // anti-proton
   G4AntiProtonAnnihilationAtRest  theAntiProtonAnnihilation;
    
   // anti-neutron
   G4AntiNeutronAnnihilationAtRest  theAntiNeutronAnnihilation;
         
};
// 2002 by J.P. Wellisch


#endif

