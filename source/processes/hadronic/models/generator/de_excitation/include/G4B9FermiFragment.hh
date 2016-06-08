// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4B9FermiFragment_h
#define G4B9FermiFragment_h 1

#include "G4UnstableFermiFragment.hh"
#include "G4IonTable.hh"

class G4B9FermiFragment : public G4UnstableFermiFragment
{
public:
  G4B9FermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE):
    G4UnstableFermiFragment(anA,aZ,Pol,ExE)
    {}; 

  ~G4B9FermiFragment();
  
private:
  G4B9FermiFragment();

  G4B9FermiFragment(const G4B9FermiFragment &right);
  
  const G4B9FermiFragment & operator=(const G4B9FermiFragment &right);
  G4bool operator==(const G4B9FermiFragment &right) const;
  G4bool operator!=(const G4B9FermiFragment &right) const;
  
public:

  G4FragmentVector * GetFragment(const G4LorentzVector & aMomentum);

};


#endif


