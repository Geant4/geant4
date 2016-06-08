// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4VFermiFragment_h
#define G4VFermiFragment_h 1

#include "G4FragmentVector.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

class G4VFermiFragment 
{
public:
  G4VFermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE):
    A(anA),
    Z(aZ),
    Polarization(Pol),
    ExcitEnergy(ExE)
    {}

  virtual ~G4VFermiFragment() {};
  
protected:
  G4VFermiFragment() {};

private:

  G4VFermiFragment(const G4VFermiFragment &right);
  
  const G4VFermiFragment & operator=(const G4VFermiFragment &right);
  G4bool operator==(const G4VFermiFragment &right) const;
  G4bool operator!=(const G4VFermiFragment &right) const;
  
public:

  virtual G4FragmentVector * GetFragment(const G4LorentzVector & aMomentum) = 0;

  G4int GetA(void) {return A;}
  G4int GetZ(void) {return Z;}
  G4int GetPolarization(void) {return Polarization;}
  G4double GetExcitationEnergy(void) {return ExcitEnergy;}

  G4double GetFragmentMass(void){
	  return G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z,A) + ExcitEnergy;
	  }

protected:

  G4int A;

  G4int Z;

  G4int Polarization;

  G4double ExcitEnergy;


};


#endif


