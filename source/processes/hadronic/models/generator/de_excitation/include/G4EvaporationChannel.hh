// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#ifndef G4EvaporationChannel_h
#define G4EvaporationChannel_h 1

#include "G4VEvaporationChannel.hh"
#include "G4VEmissionProbability.hh"
#include "G4EvaporationProbability.hh"
#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4NucleiProperties.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "g4rw/tvvector.h"

class G4EvaporationChannel : public G4VEvaporationChannel
{
public:
  // only available constructor
  G4EvaporationChannel(const G4int theGamma,
		       const G4int theA,
		       const G4int theZ,
		       G4RWTValVector<G4double> * theExcitationEnergies,
		       G4RWTValVector<G4int> * theExcitationSpins);

  // destructor
  ~G4EvaporationChannel();
  
private:
  // default constructor
  G4EvaporationChannel() {};
  // copy constructor
  G4EvaporationChannel(const G4EvaporationChannel & right);

  const G4EvaporationChannel & operator=(const G4EvaporationChannel & right);
public:
  G4bool operator==(const G4EvaporationChannel & right) const;
  G4bool operator!=(const G4EvaporationChannel & right) const;

public:
  void Initialize(const G4Fragment & fragment);

  G4FragmentVector * BreakUp(const G4Fragment & theNucleus);
  
  inline void SetEmissionStrategy(G4VEmissionProbability * aStrategy)
    {
      if (MyOwnEvaporationProbability) delete theEvaporationProbabilityPtr;
      theEvaporationProbabilityPtr = aStrategy;
      MyOwnEvaporationProbability = false;
    }


  inline void SetLevelDensityParameter(G4VLevelDensityParameter * aLevelDensity)
    {
      if (MyOwnLevelDensity) delete theLevelDensityPtr;
      theLevelDensityPtr = aLevelDensity;
      MyOwnLevelDensity = false;
    }

  inline G4double GetLevelDensityParameter(void) const { return LevelDensityParameter;}

private:

  // This data member define the channel. 
  // They are intializated at object creation (constructor) time.

  // Gamma is A_f(2S_f+1) factor, where A_f is fragment atomic number and S_f is fragment spin
  G4int Gamma;

  // Atomic Number
  G4int A;

  // Charge
  G4int Z;

  // 
  G4RWTValVector<G4double> * ExcitationEnergies;

  //
  G4RWTValVector<G4int> * ExcitationSpins;


  // For evaporation probability calcualtion
  G4bool MyOwnEvaporationProbability;
  G4VEmissionProbability * theEvaporationProbabilityPtr;

  // For Level Density calculation
  G4bool MyOwnLevelDensity;
  G4VLevelDensityParameter * theLevelDensityPtr;
  G4double LevelDensityParameter;

  //---------------------------------------------------

  // This values depends on the nucleus that is being evaporated.
  // They are calculated through the Initialize method which takes as parameters 
  // the atomic number, charge and excitation energy of nucleus.

  // Residual Atomic Number
  G4int AResidual;

  // Residual Charge
  G4int ZResidual;

  // Coulomb Barrier
  G4double CoulombBarrier;

  // Binding Energy
  G4double BindingEnergy;

  // Emission Probability
  G4double EmissionProbability;


  // Maximal Kinetic Energy that can be carried by fragment
  G4double MaximalKineticEnergy;


public:
  inline G4int GetGamma(void) const
    {return Gamma;}

  inline G4int GetA(void) const 
    {return A;}

  inline G4int GetZ(void) const 
    {return Z;}

  inline G4double GetCoulombBarrier(void) const 
    {return CoulombBarrier;}

  inline G4double GetBindingEnergy(void) const 
    {return BindingEnergy;}

  inline G4double GetEmissionProbability(void) const 
    {return EmissionProbability;}

  inline G4double GetExcitationEnergy(const G4int i) const 
    {
      if (ExcitationEnergies != 0 && i < ExcitationEnergies->length()) 
	return ExcitationEnergies->operator()(i); 
      else return 0.0;
    }

  inline G4int GetExcitationSpin(const G4int i) const 
    {
      if (ExcitationSpins != 0 && i < ExcitationSpins->length())
	return ExcitationSpins->operator()(i); 
      else return 0;
    }


  inline G4double GetMaximalKineticEnergy(void) const 
    { return MaximalKineticEnergy; }

  // ----------------------

  inline G4int GetResidualA(void) const 
    { return AResidual; }

  inline G4int GetResidualZ(void) const 
    { return ZResidual; }



private: 

  // Coulomb barrier calculation
  G4double CalcCoulombBarrier(const G4int ARes, const G4int ZRes);

  // Calculate Binding Energy for separate fragment from nucleus
  G4double CalcBindingEnergy(const G4int anA, const G4int aZ);

  // Calculate maximal kinetic energy that can be carried by fragment (in MeV)
  G4double CalcMaximalKineticEnergy(const G4double U);

  // Samples fragment kinetic energy (in MeV).
  G4double CalcKineticEnergy(void);




  // This has to be removed and put in Random Generator
  G4ThreeVector IsotropicVector(const G4double Magnitude  = 1.0);

};


#endif
