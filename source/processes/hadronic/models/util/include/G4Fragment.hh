//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4Fragment.hh 104779 2017-06-16 09:20:56Z gcosmo $
//
//---------------------------------------------------------------------
//
// Geant4 header G4Fragment
//
// by V. Lara (May 1998)
//
// Modifications:
// 03.05.2010 V.Ivanchenko General cleanup of inline functions: objects 
//            are accessed by reference; remove double return 
//            tolerance of excitation energy at modent it is computed;
//            safe computation of excitation for exotic fragments
// 18.05.2010 V.Ivanchenko added member theGroundStateMass and inline
//            method which allowing to compute this value once and use 
//            many times
// 26.09.2010 V.Ivanchenko added number of protons, neutrons, proton holes
//            and neutron holes as members of the class and Get/Set methods;
//            removed not needed 'const'; removed old debug staff and unused
//            private methods; add comments and reorder methods for 
//            better reading

#ifndef G4Fragment_h
#define G4Fragment_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4NuclearPolarization.hh"
#include "G4NucleiProperties.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include <vector>

class G4ParticleDefinition;

class G4Fragment;     
typedef std::vector<G4Fragment*> G4FragmentVector;

class G4Fragment 
{
public:

  // ============= CONSTRUCTORS ==================

  // Default constructor - obsolete
  G4Fragment();

  // Destructor
  ~G4Fragment();

  // Copy constructor
  G4Fragment(const G4Fragment &right);

  // A,Z and 4-momentum - main constructor for fragment
  G4Fragment(G4int A, G4int Z, const G4LorentzVector& aMomentum);

  // 4-momentum and pointer to G4particleDefinition (for gammas, e-)
  G4Fragment(const G4LorentzVector& aMomentum, 
	     const G4ParticleDefinition* aParticleDefinition);

  // ============= OPERATORS ==================
    
  G4Fragment & operator=(const G4Fragment &right);
  G4bool operator==(const G4Fragment &right) const;
  G4bool operator!=(const G4Fragment &right) const;

  friend std::ostream& operator<<(std::ostream&, const G4Fragment&);

  //  new/delete operators are overloded to use G4Allocator
  inline void *operator new(size_t);
  inline void operator delete(void *aFragment);

  // ============= GENERAL METHODS ==================

  inline G4int GetZ_asInt() const;
  inline G4int GetA_asInt() const;
  inline void SetZandA_asInt(G4int Znew, G4int Anew);
  
  inline G4double GetExcitationEnergy() const;
  inline void SetExcEnergyAndMomentum(G4double eexc, const G4LorentzVector&);

  inline G4double GetGroundStateMass() const;
   
  inline G4double GetBindingEnergy() const;
  
  inline const G4LorentzVector& GetMomentum() const;
  inline void SetMomentum(const G4LorentzVector& value);
  
  // computation of mass for any Z and A
  inline G4double ComputeGroundStateMass(G4int Z, G4int A) const;

  // extra methods
  inline G4double GetSpin() const;
  inline void SetSpin(G4double value);

  inline G4int GetCreatorModelType() const;
  inline void SetCreatorModelType(G4int value);

  // obsolete methods
  
  inline G4double GetZ() const;
  inline G4double GetA() const;
  inline void SetZ(G4double value);
  inline void SetA(G4double value);
  
  // ============= METHODS FOR PRE-COMPOUND MODEL ===============

  inline G4int GetNumberOfExcitons() const;
  
  inline G4int GetNumberOfParticles() const;
  inline G4int GetNumberOfCharged() const;
  inline void SetNumberOfExcitedParticle(G4int valueTot, G4int valueP);

  inline G4int GetNumberOfHoles() const;
  inline G4int GetNumberOfChargedHoles() const;
  inline void SetNumberOfHoles(G4int valueTot, G4int valueP=0);
  
  // these methods will be removed in future
  inline void SetNumberOfParticles(G4int value);
  inline void SetNumberOfCharged(G4int value);

  // ============= METHODS FOR PHOTON EVAPORATION ===============

  inline G4int GetNumberOfElectrons() const;
  inline void SetNumberOfElectrons(G4int value);

  inline G4int GetFloatingLevelNumber() const;
  inline void SetFloatingLevelNumber(G4int value);

  inline const G4ParticleDefinition * GetParticleDefinition() const;
  inline void SetParticleDefinition(const G4ParticleDefinition * p);

  inline G4double GetCreationTime() const;
  inline void SetCreationTime(G4double time);

  inline G4NuclearPolarization* NuclearPolarization();
  inline G4NuclearPolarization* GetNuclearPolarization() const;
  void SetNuclearPolarization(G4NuclearPolarization*);

  void SetAngularMomentum(const G4ThreeVector&);
  G4ThreeVector GetAngularMomentum() const;

  // ============= PRIVATE METHODS ==============================

private:

  void ExcitationEnergyWarning();

  void NumberOfExitationWarning(const G4String&);

  inline void CalculateExcitationEnergy();

  inline void CalculateGroundStateMass();

  // ============= DATA MEMBERS ==================

  G4int theA;
  
  G4int theZ;
  
  G4double theExcitationEnergy;

  G4double theGroundStateMass;

  G4LorentzVector theMomentum;
  
  // Nuclear polarisation by default is nullptr
  G4NuclearPolarization* thePolarization;

  // creator model type
  G4int creatorModel;

  // Exciton model data members  
  G4int numberOfParticles;  
  G4int numberOfCharged;
  G4int numberOfHoles;
  G4int numberOfChargedHoles;

  // Gamma evaporation data members
  G4int numberOfShellElectrons;
  G4int xLevel;

  const G4ParticleDefinition* theParticleDefinition;
  
  G4double spin;
  G4double theCreationTime;

  static const G4double minFragExcitation;
};

// ============= INLINE METHOD IMPLEMENTATIONS ===================
/*
inline void G4Fragment::SetNuclearPolarization(G4NuclearPolarization* p)
{
  if(p !=  thePolarization) {
    delete thePolarization;
    thePolarization = p;
  }
}
*/
#if defined G4HADRONIC_ALLOC_EXPORT
  extern G4DLLEXPORT G4ThreadLocal G4Allocator<G4Fragment> *pFragmentAllocator;
#else
  extern G4DLLIMPORT G4ThreadLocal G4Allocator<G4Fragment> *pFragmentAllocator;
#endif

inline void * G4Fragment::operator new(size_t)
{
  if (!pFragmentAllocator) { pFragmentAllocator = new G4Allocator<G4Fragment>; }
  return (void*) pFragmentAllocator->MallocSingle();
}

inline void G4Fragment::operator delete(void * aFragment)
{
  ((G4Fragment *)aFragment)->SetNuclearPolarization(nullptr);
  pFragmentAllocator->FreeSingle((G4Fragment *) aFragment);
}

inline void G4Fragment::CalculateExcitationEnergy()
{
  theExcitationEnergy = theMomentum.mag() - theGroundStateMass;
  if(theExcitationEnergy < minFragExcitation) { 
    if(theExcitationEnergy < -minFragExcitation) {  ExcitationEnergyWarning(); }
    theExcitationEnergy = 0.0;
  }
}

inline G4double 
G4Fragment::ComputeGroundStateMass(G4int Z, G4int A) const
{
  return G4NucleiProperties::GetNuclearMass(A, Z); 
}
	 
inline void G4Fragment::CalculateGroundStateMass() 
{
  theGroundStateMass = G4NucleiProperties::GetNuclearMass(theA, theZ);
}

inline G4int G4Fragment::GetA_asInt() const
{
  return theA;
}

inline G4int G4Fragment::GetZ_asInt()  const
{
  return theZ;
}

inline void G4Fragment::SetZandA_asInt(G4int Znew, G4int Anew)
{
  theZ = Znew;
  theA = Anew;
  CalculateGroundStateMass();
}

inline G4double G4Fragment::GetExcitationEnergy()  const
{
  return theExcitationEnergy;
}

inline G4double G4Fragment::GetGroundStateMass() const
{
  return theGroundStateMass; 
}

inline void G4Fragment::SetExcEnergyAndMomentum(G4double eexc, 
						const G4LorentzVector& v)
{
  theExcitationEnergy = eexc;
  theMomentum.set(0.0, 0.0, 0.0, theGroundStateMass + eexc);
  theMomentum.boost(v.boostVector());
}

inline G4double G4Fragment::GetBindingEnergy() const
{
  return (theA-theZ)*CLHEP::neutron_mass_c2 + theZ*CLHEP::proton_mass_c2 
    - theGroundStateMass;
}

inline const G4LorentzVector& G4Fragment::GetMomentum()  const
{
  return theMomentum;
}

inline void G4Fragment::SetMomentum(const G4LorentzVector& value)
{
  theMomentum = value;
  CalculateExcitationEnergy();
}

inline G4double G4Fragment::GetZ()  const
{
  return G4double(theZ);
}

inline G4double G4Fragment::GetA() const
{
  return G4double(theA);
}

inline void G4Fragment::SetZ(const G4double value)
{
  theZ = G4lrint(value);
  CalculateGroundStateMass();
}

inline void G4Fragment::SetA(const G4double value)
{
  theA = G4lrint(value);
  CalculateGroundStateMass();
}

inline G4int G4Fragment::GetNumberOfExcitons()  const
{
  return numberOfParticles + numberOfHoles;
}

inline G4int G4Fragment::GetNumberOfParticles()  const
{
  return numberOfParticles;
}

inline G4int G4Fragment::GetNumberOfCharged()  const
{
  return numberOfCharged;
}

inline 
void G4Fragment::SetNumberOfExcitedParticle(G4int valueTot, G4int valueP)
{
  numberOfParticles = valueTot;
  numberOfCharged = valueP;
  if(valueTot < valueP)  { 
    NumberOfExitationWarning("SetNumberOfExcitedParticle"); 
  }
}

inline G4int G4Fragment::GetNumberOfHoles()  const
{
  return numberOfHoles;
}

inline G4int G4Fragment::GetNumberOfChargedHoles()  const
{
  return numberOfChargedHoles;
}

inline void G4Fragment::SetNumberOfHoles(G4int valueTot, G4int valueP)
{
  numberOfHoles = valueTot;
  numberOfChargedHoles = valueP;
  if(valueTot < valueP)  { 
    NumberOfExitationWarning("SetNumberOfHoles"); 
  }
}

inline void G4Fragment::SetNumberOfParticles(G4int value)
{
  numberOfParticles = value;
}

inline void G4Fragment::SetNumberOfCharged(G4int value)
{
  numberOfCharged = value;
  if(value > numberOfParticles)  { 
    NumberOfExitationWarning("SetNumberOfCharged"); 
  }
}

inline G4int G4Fragment::GetNumberOfElectrons() const
{
  return numberOfShellElectrons;
}

inline void G4Fragment::SetNumberOfElectrons(G4int value)
{
  numberOfShellElectrons = value;
}

inline G4int G4Fragment::GetCreatorModelType() const
{
  return creatorModel;
}

inline void G4Fragment::SetCreatorModelType(G4int value)
{
  creatorModel = value;
}

inline G4double G4Fragment::GetSpin() const
{
  return spin;
}

inline void G4Fragment::SetSpin(G4double value)
{
  spin = value;
}

inline G4int G4Fragment::GetFloatingLevelNumber() const
{
  return xLevel;
}

inline void G4Fragment::SetFloatingLevelNumber(G4int value)
{
  xLevel = value;
}

inline 
const G4ParticleDefinition* G4Fragment::GetParticleDefinition(void) const
{
  return theParticleDefinition;
}

inline void G4Fragment::SetParticleDefinition(const G4ParticleDefinition * p)
{
  theParticleDefinition = p;
}

inline G4double G4Fragment::GetCreationTime() const 
{
  return theCreationTime;
}

inline void G4Fragment::SetCreationTime(G4double time)
{
  theCreationTime = time;
}

inline G4NuclearPolarization* G4Fragment::NuclearPolarization()
{
  return thePolarization;
}

inline G4NuclearPolarization* G4Fragment::GetNuclearPolarization() const
{
  return thePolarization;
}

#endif


