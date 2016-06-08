// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara


#ifndef G4PreCompoundNeutron_h
#define G4PreCompoundNeutron_h 1

#include "G4VPreCompoundNucleon.hh"
#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"
#include "G4PreCompoundParameters.hh"
#include "Randomize.hh"



class G4PreCompoundNeutron : public G4VPreCompoundNucleon
{
public:
  // default constructor
  G4PreCompoundNeutron() : G4VPreCompoundNucleon(1,0) {};

  // copy constructor
  G4PreCompoundNeutron(const G4PreCompoundNeutron &right):
    G4VPreCompoundNucleon(right) {};

  ~G4PreCompoundNeutron() {};

  // operators  
  const G4PreCompoundNeutron & operator=(const G4PreCompoundNeutron &right) {
    if (&right != this) this->G4VPreCompoundNucleon::operator=(right);
    return *this;
  };

  G4bool operator==(const G4PreCompoundNeutron &right) const
  {return G4VPreCompoundNucleon::operator==(right);};
  
  G4bool operator!=(const G4PreCompoundNeutron &right) const
  {return G4VPreCompoundNucleon::operator!=(right);};


  const G4DynamicParticle GetDynamicParticle() const
    {
      G4DynamicParticle theDynamicParticle(G4Neutron::NeutronDefinition(),GetMomentum());
      return theDynamicParticle;
    }



public:
  G4double ProbabilityDistributionFunction(const G4double & eKin,
					   const G4Fragment & aFragment);
  // Gives the kinetic energy for fragments in pre-equilibrium decay
  G4double GetKineticEnergy(const G4Fragment & aFragment);

};

#endif
 
