#ifndef G4NeutronIsotopeProduction_h
#define G4NeutronIsotopeProduction_h

#include "globals.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronElementIsoCrossSections.hh"

class G4NeutronIsotopeProduction : public G4VDiscreteProcess
{
  public:
  
  G4NeutronIsotopeProduction();
  ~G4NeutronIsotopeProduction();

      G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

      G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );

  private:
    
  G4ParticleChangeForIsotopeProduction theParticleChange;
  G4NeutronElementIsoCrossSections ** theData;
  G4int numberOfElements;
};

#endif
