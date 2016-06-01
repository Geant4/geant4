// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4Li5FermiFragment_h
#define G4Li5FermiFragment_h 1

#include "G4UnstableFermiFragment.hh"
#include "G4IonTable.hh"

class G4Li5FermiFragment : public G4UnstableFermiFragment
{
public:
  G4Li5FermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE):
    G4UnstableFermiFragment(anA,aZ,Pol,ExE)
    {}; 

  ~G4Li5FermiFragment();
  
private:
  G4Li5FermiFragment();

  G4Li5FermiFragment(const G4Li5FermiFragment &right);
  
  const G4Li5FermiFragment & operator=(const G4Li5FermiFragment &right);
  G4bool operator==(const G4Li5FermiFragment &right) const;
  G4bool operator!=(const G4Li5FermiFragment &right) const;
  
public:

  G4FragmentVector * GetFragment(const G4LorentzVector & aMomentum);


};


#endif


