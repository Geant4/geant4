// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#ifndef G4PreCompoundProton_h
#define G4PreCompoundProton_h 1

#include "G4VPreCompoundNucleon.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4PreCompoundParameters.hh"
#include "Randomize.hh"


class G4PreCompoundProton : public G4VPreCompoundNucleon
{
public:
  // default constructor
  G4PreCompoundProton():G4VPreCompoundNucleon(1,1) {};

  // copy constructor
  G4PreCompoundProton(const G4PreCompoundProton &right):
    G4VPreCompoundNucleon(right) {};

  ~G4PreCompoundProton() {};

  // operators  
  const G4PreCompoundProton & operator=(const G4PreCompoundProton &right) {
    if (&right != this) this->G4VPreCompoundNucleon::operator=(right);
    return *this;
  };

  G4bool operator==(const G4PreCompoundProton &right) const
  {return G4VPreCompoundNucleon::operator==(right);};

  
  G4bool operator!=(const G4PreCompoundProton &right) const
  {return G4VPreCompoundNucleon::operator!=(right);};


  const G4DynamicParticle GetDynamicParticle() const
    {
      G4DynamicParticle theDynamicParticle(G4Proton::ProtonDefinition(),GetMomentum());
      return theDynamicParticle;
    }


public:
  G4double ProbabilityDistributionFunction(const G4double & eKin,
					   const G4Fragment & aFragment);

  // Gives the kinetic energy for fragments in pre-equilibrium decay
  G4double GetKineticEnergy(const G4Fragment & aFragment);

};

#endif
 
