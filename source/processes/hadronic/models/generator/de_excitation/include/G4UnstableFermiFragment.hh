// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4UnstableFermiFragment_h
#define G4UnstableFermiFragment_h 1

#include "G4VFermiFragment.hh"
#include "Randomize.hh"

class G4UnstableFermiFragment : public G4VFermiFragment
{
public:
  G4UnstableFermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE):
    G4VFermiFragment(anA,aZ,Pol,ExE)
    {}; 

  ~G4UnstableFermiFragment();
  
protected:
  G4UnstableFermiFragment();
private:
  G4UnstableFermiFragment(const G4UnstableFermiFragment &right);
  
  const G4UnstableFermiFragment & operator=(const G4UnstableFermiFragment &right);
  G4bool operator==(const G4UnstableFermiFragment &right) const;
  G4bool operator!=(const G4UnstableFermiFragment &right) const;
  
public:

  G4RWTPtrOrderedVector<G4LorentzVector> *
  FragmentsMomentum(G4double KinE, const G4int K, const G4double * Masses);

private:

  G4double RNKSI(const G4int K);

  G4ParticleMomentum IsotropicVector(const G4double Magnitude = 1.0);


};


#endif


