// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
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
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Fragment.hh"

class G4DynamicParticle;

class G4VPreCompoundFragment
{
	// ============================
	// Constructors and destructor
	// ============================

protected:
	// default constructor
	G4VPreCompoundFragment() {};

public:  
	// copy constructor
	G4VPreCompoundFragment(const G4VPreCompoundFragment &right);

	// constructor  
	G4VPreCompoundFragment(const G4double anA, const G4double aZ);

	virtual ~G4VPreCompoundFragment();

	// ==========
	// operators 
	// ========== 
	
	const G4VPreCompoundFragment& operator=(const G4VPreCompoundFragment &right);

	G4int operator==(const G4VPreCompoundFragment &right) const;
  
	G4int operator!=(const G4VPreCompoundFragment &right) const;

	friend G4std::ostream& operator<<(G4std::ostream&, const G4VPreCompoundFragment*);
	friend G4std::ostream& operator<<(G4std::ostream&, const G4VPreCompoundFragment&);

	// =====================
	// Pure Virtual methods
	// =====================
	virtual void CalcExcitonLevelDensityRatios(const G4double Excitons, const G4double Particles) = 0;
  
	virtual G4double GetKineticEnergy(const G4Fragment & aFragment) = 0;

	// Calculates condensation probabilities to create clusters  
	// consisting of N nucleons inside a nucleus with A nucleons
	virtual void CalcCondensationProbability(const G4double A) = 0;
	
	virtual const G4DynamicParticle GetDynamicParticle() const = 0; 	

protected:
  virtual G4double ProbabilityDistributionFunction(const G4double & K, const G4Fragment & aFragment) = 0;

public:

	// =====================
	// Initialization method
	// =====================
	void Init(const G4Fragment & aFragment);
  
	// ================================================
	// Methods for calculating the emission probability
	// ================================================

	// Calculates the total (integrated over kinetic energy) emission
	// probability of a fragment
	G4double CalcEmissionProbability(const G4Fragment & aFragment);

	// See above (in virtual methods) the method ProbabilityDistributionFunction
private:	
	// This method performs integration for probability function over 
	// fragment kinetic energy
	G4double IntegrateEmissionProbability(const G4double & Low, const G4double & Up, 
  													 const G4Fragment & aFragment);	

	// ========================================
	// Method for calculate the Coulomb barrier
	// ========================================
private:
	G4double CalcCoulombBarrier(const G4double NucRad, const G4double aZ);

	// ============================
	// Data members access methods
	// ============================

public:
	const G4double GetA() const { return theA;}

	const G4double GetZ() const { return theZ;}
  
	const G4double GetRestA() const { return theRestNucleusA;}

	const G4double GetRestZ() const { return theRestNucleusZ;}
  
	const G4double GetCoulombBarrier() const {return theCoulombBarrier;}

	const G4double GetBindingEnergy() const { return theBindingEnergy;}

	const G4double GetMaximalKineticEnergy() const { return theMaximalKineticEnergy;}

	const G4double GetExcitonLevelDensityRatio() const { return theExcitonLevelDensityRatio;}
	void SetExcitonLevelDensityRatio(const G4double value) { theExcitonLevelDensityRatio = value;} 
  
	void SetEmissionProbability(const G4double value) { theEmissionProbability = value;}
	
	const G4double GetEmissionProbability() const { return theEmissionProbability;}

	const G4double GetCondensationProbability() const { return theCondensationProbability;}
	void SetCondensationProbability(const G4double value) { theCondensationProbability = value;} 

	const G4double GetNuclearMass() const {
		return G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(theZ,theA);
	}

	const G4double GetRestNuclearMass() const {
		return G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(theRestNucleusZ,theRestNucleusA);
	}



	const G4double GetReducedMass() const {
		return GetRestNuclearMass()*GetNuclearMass()/(GetNuclearMass()+GetRestNuclearMass());
	}

	const G4LorentzVector GetMomentum() const { return theMomentum;}

	void SetMomentum(const G4LorentzVector & value) { theMomentum = value;}

protected:

	virtual G4double GetBarrierPenetrationFactor(const G4double aZ) const
	{return 1.0;}

//	virtual G4double Get 



	// =============
	// Data members
	// =============
	
private:

  G4double theA;

  G4double theZ;

  G4double theRestNucleusA;

  G4double theRestNucleusZ;

  G4double theCoulombBarrier;

  G4double theBindingEnergy;

  G4double theMaximalKineticEnergy;

  G4double theExcitonLevelDensityRatio;

  G4double theEmissionProbability;

  G4double theCondensationProbability;

  G4LorentzVector theMomentum;

};

#endif
