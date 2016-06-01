// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara


#ifndef G4PreCompoundIon_h
#define G4PreCompoundIon_h 1


#include "G4VPreCompoundFragment.hh"
#include "G4PreCompoundParameters.hh"
#include "Randomize.hh"


class G4VPreCompoundIon : public G4VPreCompoundFragment
{
protected:
  // default constructor
  G4VPreCompoundIon() {};

public:

  // copy constructor
  G4VPreCompoundIon(const G4VPreCompoundIon &right):
    G4VPreCompoundFragment(right) {};

  // constructor  
  G4VPreCompoundIon(const G4double anA, const G4double aZ):
    G4VPreCompoundFragment(anA,aZ) {};

  virtual ~G4VPreCompoundIon() {};

  // operators  
  const G4VPreCompoundIon & operator=(const G4VPreCompoundIon &right) {
    if (&right != this) this->G4VPreCompoundFragment::operator=(right);
    return *this;
  };


  G4bool operator==(const G4VPreCompoundIon &right) const
  {return G4VPreCompoundFragment::operator==(right);};

  G4bool operator!=(const G4VPreCompoundIon &right) const
  {return G4VPreCompoundFragment::operator!=(right);};



public:
  G4double ProbabilityDistributionFunction(const G4double & eKin,
					   const G4Fragment & aFragment);

  // Gives the kinetic energy for fragments in pre-equilibrium decay
  G4double GetKineticEnergy(const G4Fragment & aFragment);


};

#endif
 
