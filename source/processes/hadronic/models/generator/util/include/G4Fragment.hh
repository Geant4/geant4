// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Fragment.hh,v 1.10 1998/12/16 11:46:47 larazb Exp $
// GEANT4 tag $Name: geant4-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)

#ifndef G4Fragment_h
#define G4Fragment_h 1

#include "G4ios.hh"
#include <iomanip.h>
#include <rw/tpordvec.h>

#include "globals.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4ThreeVector.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

class G4ParticleDefinition;

class G4Fragment;     // Forward deckaration
typedef RWTPtrOrderedVector<G4Fragment> G4FragmentVector;

class G4Fragment 
{
public:
  G4Fragment();
  ~G4Fragment();
  G4Fragment(const G4Fragment &right);
  G4Fragment(const G4int A, const G4int Z, const G4LorentzVector aMomentum);
  G4Fragment(const G4LorentzVector aMomentum, G4ParticleDefinition * aParticleDefinition);
  G4Fragment(const G4int A, const G4int Z, const G4double anExEnergy, const G4LorentzVector aMomentum,
	     const G4ThreeVector anAngularMomentum, const G4int aNumberOfExcitons, const G4int aNumberOfHoles,
	     const G4int aNumberOfCharged);
  G4Fragment(const G4int A, const G4int Z, const G4double anExEnergy, 
	     const G4LorentzVector aMomentum, 
	     const G4ThreeVector anAngularMomentum, const G4int aNumberOfExcitons, 
	     const G4int aNumberOfHoles, const G4int aNumberOfCharged,
	     G4ParticleDefinition * aParticleDefinition);
    
  const G4Fragment & operator=(const G4Fragment &right);
  G4bool operator==(const G4Fragment &right) const;
  G4bool operator!=(const G4Fragment &right) const;

  friend ostream& operator<<(ostream&, const G4Fragment*);
  friend ostream& operator<<(ostream&, const G4Fragment&);


  inline G4double GetA(void) const;
  void SetA(const G4double value);
  
  G4double GetZ(void) const;
  void SetZ(const G4double value);
  
  G4double GetExcitationEnergy(void) const;
  void SetExcitationEnergy(const G4double value);
  
  const G4LorentzVector GetMomentum(void) const;
  void SetMomentum(const G4LorentzVector value);
  
  const G4ThreeVector GetAngularMomentum(void) const;
  void SetAngularMomentum(const G4ThreeVector value);
  
  G4int GetNumberOfExcitons(void) const;
  void SetNumberOfExcitons(const G4int value);
  
  G4int GetNumberOfHoles(void) const;
  void SetNumberOfHoles(const G4int value);
  
  G4int GetNumberOfCharged(void) const;
  void SetNumberOfCharged(const G4int value);

  G4int GetNumberOfParticles(void) const;


  inline G4ParticleDefinition * GetParticleDefinition(void) const;
  void SetParticleDefinition(G4ParticleDefinition * aParticleDefinition);

  void DumpInfo() const;
  
private:

  G4double theA;
  
  G4double theZ;
  
  G4double theExcitationEnergy;

  G4LorentzVector theMomentum;
  
  G4ThreeVector theAngularMomentum;
  
  G4int numberOfExcitons;
  
  G4int numberOfHoles;
  
  G4int numberOfCharged;

  G4ParticleDefinition * theParticleDefinition;
  
};

// Class G4Fragment 

inline G4double G4Fragment::GetA() const
{
  return theA;
}

inline void G4Fragment::SetA(const G4double value)
{
  theA = value;
}

inline G4double G4Fragment::GetZ()  const
{
  return theZ;
}

inline void G4Fragment::SetZ(const G4double value)
{
  theZ = value;
}

inline G4double G4Fragment::GetExcitationEnergy()  const
{
  return theExcitationEnergy;
}

inline void G4Fragment::SetExcitationEnergy(const G4double value)
{
  theExcitationEnergy = value;
}

inline const G4LorentzVector G4Fragment::GetMomentum()  const
{
  return theMomentum;
}

inline void G4Fragment::SetMomentum(const G4LorentzVector value)
{
  theMomentum = value;
}

inline const G4ThreeVector G4Fragment::GetAngularMomentum()  const
{
  return theAngularMomentum;
}

inline void G4Fragment::SetAngularMomentum(const G4ThreeVector value)
{
  theAngularMomentum = value;
}

inline G4int G4Fragment::GetNumberOfExcitons()  const
{
  return numberOfExcitons;
}

inline void G4Fragment::SetNumberOfExcitons(const G4int value)
{
  numberOfExcitons = value;
}

inline G4int G4Fragment::GetNumberOfHoles()  const
{
  return numberOfHoles;
}

inline void G4Fragment::SetNumberOfHoles(const G4int value)
{
  numberOfHoles = value;
}

inline G4int G4Fragment::GetNumberOfCharged()  const
{
  return numberOfCharged;
}

inline void G4Fragment::SetNumberOfCharged(const G4int value)
{
  numberOfCharged = value;
}

inline G4int G4Fragment::GetNumberOfParticles()  const
{
  return numberOfExcitons - numberOfHoles;
}


inline G4ParticleDefinition * G4Fragment::GetParticleDefinition(void) const
{
  return theParticleDefinition;
}

inline void G4Fragment::SetParticleDefinition(G4ParticleDefinition * aParticleDefinition)
{
  theParticleDefinition = aParticleDefinition;
}



inline void G4Fragment::DumpInfo() const
{
  G4cout << this;
}
#endif


