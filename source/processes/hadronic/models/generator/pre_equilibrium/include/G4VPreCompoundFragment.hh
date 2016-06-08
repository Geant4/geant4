// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#ifndef G4VPreCompoundFragment_h
#define G4VPreCompoundFragment_h 1

#include "G4ios.hh"
#include "g4std/iomanip"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
//#include "G4DynamicParticle.hh"

#include "G4Fragment.hh"

class G4DynamicParticle;

class G4VPreCompoundFragment
{
protected:
  // default constructor
  G4VPreCompoundFragment() {};

public:  
  // copy constructor
  G4VPreCompoundFragment(const G4VPreCompoundFragment &right);

  // constructor  
  G4VPreCompoundFragment(const G4double anA, const G4double aZ);

  virtual ~G4VPreCompoundFragment();

  // operators  
  const G4VPreCompoundFragment& operator=(const G4VPreCompoundFragment &right);

  G4int operator==(const G4VPreCompoundFragment &right) const;
  
  G4int operator!=(const G4VPreCompoundFragment &right) const;

  friend G4std::ostream& operator<<(G4std::ostream&, const G4VPreCompoundFragment*);
  friend G4std::ostream& operator<<(G4std::ostream&, const G4VPreCompoundFragment&);


  // methods

  void Init(const G4Fragment & aFragment);
  
  virtual void CalcExcitonLevelDensityRatios(const G4double Excitons,
			                     const G4double Particles) = 0;
  
  virtual G4double GetKineticEnergy(const G4Fragment & aFragment) = 0;

  // Calculates condensation probabilities to create fragment consisting from Nf nucleons 
  // inside a nucleus with A nucleons
  virtual void CalcCondensationProbability(const G4double A) = 0;

  // Calculates the total (integrated over kinetic energy) emission
  // probability of a fragment
  G4double CalcEmissionProbability(const G4Fragment & aFragment);

  void SetA(const G4double value);
  const G4double GetA() const;

  void SetZ(const G4double value);
  const G4double GetZ() const;
  

  void SetRestA(const G4double value);
  const G4double GetRestA() const;

  void SetRestZ(const G4double value);
  const G4double GetRestZ() const;
  

  void SetCoulombBarrier(const G4double value);
  const G4double GetCoulombBarrier() const;

  void SetBindingEnergy(const G4double value);
  const G4double GetBindingEnergy() const;

  void SetMaximalKineticEnergy(const G4double value);
  const G4double GetMaximalKineticEnergy() const;

  void SetExcitonLevelDensityRatio(const G4double value);
  const G4double GetExcitonLevelDensityRatio() const;

  void SetEmissionProbability(const G4double value);
  const G4double GetEmissionProbability() const;

  void SetCondensationProbability(const G4double value);
  const G4double GetCondensationProbability() const;


  const G4double GetNuclearMass() const;
  const G4double GetRestNuclearMass() const;

protected:
  virtual G4double ProbabilityDistributionFunction(const G4double & K,
					   const G4Fragment & aFragment) = 0;
private:
  G4double CalcCoulombBarrier(const G4double & NucRad);
  
   // This method performs integration for probability function over 
   // fragment kinetic energy
  G4double IntegrateEmissionProbability(const G4double & Low, const G4double & Up, 
  					const G4Fragment & aFragment);

public:
  void SetMomentum(const G4LorentzVector value);
  const G4LorentzVector GetMomentum() const;

  virtual const G4DynamicParticle GetDynamicParticle() const = 0; 

private:

  G4double theA;

  G4double theZ;

  G4double theRestNucleusA;

  G4double theRestNucleusZ;

  G4double CoulombBarrier;

  G4double BindingEnergy;

  G4double MaximalKineticEnergy;

  G4double ExcitonLevelDensityRatio;

  G4double EmissionProbability;

  G4double CondensationProbability;

  G4LorentzVector Momentum;

};


inline void G4VPreCompoundFragment::SetA(const G4double value)
{
  theA = value;
}

inline const G4double G4VPreCompoundFragment::GetA() const
{
  return theA;
}

inline void G4VPreCompoundFragment::SetZ(const G4double value)
{
  theZ = value;
}

inline const G4double G4VPreCompoundFragment::GetZ() const
{
  return theZ;
}

inline void G4VPreCompoundFragment::SetRestA(const G4double value)
{
  theRestNucleusA = value - theA;
}

inline const G4double G4VPreCompoundFragment::GetRestA() const
{
  return theRestNucleusA;
}

inline void G4VPreCompoundFragment::SetRestZ(const G4double value)
{
  theRestNucleusZ = value - theZ;
}

inline const G4double G4VPreCompoundFragment::GetRestZ() const
{
  return theRestNucleusZ;
}

inline void G4VPreCompoundFragment::SetCoulombBarrier(const G4double value)
{
  CoulombBarrier = value;
}

inline const G4double G4VPreCompoundFragment::GetCoulombBarrier() const
{
  return CoulombBarrier;
}

inline void G4VPreCompoundFragment::SetBindingEnergy(const G4double value)
{
  BindingEnergy = value;
}

inline const G4double G4VPreCompoundFragment::GetBindingEnergy() const
{
  return BindingEnergy;
}


inline void G4VPreCompoundFragment::SetMaximalKineticEnergy(const G4double value)
{
  MaximalKineticEnergy = value;
}

inline const G4double G4VPreCompoundFragment::GetMaximalKineticEnergy() const
{
  return MaximalKineticEnergy;
}

inline void G4VPreCompoundFragment::SetExcitonLevelDensityRatio(const G4double value)
{
  ExcitonLevelDensityRatio = value;
}

inline const G4double G4VPreCompoundFragment::GetExcitonLevelDensityRatio() const
{
  return ExcitonLevelDensityRatio;
}


inline void G4VPreCompoundFragment::SetEmissionProbability(const G4double value)
{
  EmissionProbability = value;
}

inline const G4double G4VPreCompoundFragment::GetEmissionProbability() const
{
  return EmissionProbability;
}

inline void G4VPreCompoundFragment::SetCondensationProbability(const G4double value)
{
  CondensationProbability = value;
}

inline const G4double G4VPreCompoundFragment::GetCondensationProbability() const
{
  return CondensationProbability;
}


inline const G4double G4VPreCompoundFragment::GetNuclearMass() const
  // Calculate nucleus atomic mass (MeV)
{
  return G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(theZ,theA)/MeV;
}

inline const G4double G4VPreCompoundFragment::GetRestNuclearMass() const
  // Calculate nucleus atomic mass (MeV)
{
  return G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(theRestNucleusZ,theRestNucleusA)/MeV;
}


inline void G4VPreCompoundFragment::SetMomentum(const G4LorentzVector value)
{
  Momentum = value;
}

inline const G4LorentzVector G4VPreCompoundFragment::GetMomentum() const
{
  return Momentum;
}



#endif
