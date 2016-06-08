// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#ifndef G4VPreCompoundNucleon_h
#define G4VPreCompoundNucleon_h 1

#include "G4VPreCompoundFragment.hh"


class G4VPreCompoundNucleon : public G4VPreCompoundFragment
{
protected:
  // copy constructor
  G4VPreCompoundNucleon() {};

public:

  // copy constructor
  G4VPreCompoundNucleon(const G4VPreCompoundNucleon &right):
    G4VPreCompoundFragment(right) {};

  // constructor  
  G4VPreCompoundNucleon(const G4double anA, const G4double aZ):
    G4VPreCompoundFragment(anA,aZ) {};

  virtual ~G4VPreCompoundNucleon() {};

  // operators  
  const G4VPreCompoundNucleon & operator=(const G4VPreCompoundNucleon &right) {
    if (&right != this) this->G4VPreCompoundFragment::operator=(right);
    return *this;
  };

  G4bool operator==(const G4VPreCompoundNucleon &right) const 
  {return G4VPreCompoundFragment::operator==(right); };
  
  G4bool operator!=(const G4VPreCompoundNucleon &right) const 
  {return G4VPreCompoundFragment::operator!=(right); };


  void CalcExcitonLevelDensityRatios(const G4double Excitons, 
				     const G4double Particles)
  {
    // Level density ratios are calculated according to the formula
    // (P!*(N-1)!)/((P-Af)!*(N-1-Af)!*Af!)
    // where  P is number of particles
    //        N is number of excitons
    //        Af atomic number of emitting fragment
    // the next is a simplification for nucleons (Af = 1)

    SetExcitonLevelDensityRatio(Particles*(Excitons-1.0));
  }
   
  void CalcCondensationProbability(const G4double A)
    // This method computes condensation probability to create a fragment
    // consisting from N nucleons inside a nucleus with A nucleons 
    // This value comes from the formula N^3 (N/A)^(N-1) with N = 1 (nucleon)
  {
    SetCondensationProbability(1.0);
  }

};

#endif
 
