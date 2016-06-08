// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4StableFermiFragment.hh"


G4StableFermiFragment::G4StableFermiFragment()
{
}

G4StableFermiFragment::G4StableFermiFragment(const G4StableFermiFragment &right)
{
  G4Exception("G4StableFermiFragment::copy_constructor meant to not be accessable");
}


G4StableFermiFragment::~G4StableFermiFragment()
{
}


const G4StableFermiFragment & G4StableFermiFragment::operator=(const G4StableFermiFragment &right)
{
  G4Exception("G4StableFermiFragment::operator= meant to not be accessable");
  return *this;
}


G4bool G4StableFermiFragment::operator==(const G4StableFermiFragment &right) const
{
  return false;
}

G4bool G4StableFermiFragment::operator!=(const G4StableFermiFragment &right) const
{
  return true;
}



G4FragmentVector * G4StableFermiFragment::GetFragment(const G4LorentzVector & aMomentum)
{
  G4FragmentVector * theResult = new G4FragmentVector;

  theResult->insert(new G4Fragment(A,Z,aMomentum));
  
  return theResult;
}
