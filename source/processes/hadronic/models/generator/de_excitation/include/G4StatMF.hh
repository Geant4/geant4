#ifndef G4StatMF_h
#define G4StatMF_h 1

#include <rw/tvordvec.h>

#include "globals.hh"
#include "G4MultiFragmentation.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4StatMFFragment.hh"
#include "G4StatMFParameters.hh"
#include "G4VStatMFCanonical.hh"
#include "G4StatMFMicrocanonical.hh"
#include "G4StatMFMacrocanonical.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"



class G4StatMF : public G4MultiFragmentation
{
public:
  G4StatMF();
  ~G4StatMF();

private:
  G4StatMF(const G4StatMF & right);
  G4StatMF & operator=(const G4StatMF & right);
  G4bool operator==(const G4StatMF & right);
  G4bool operator!=(const G4StatMF & right);

public:
  G4FragmentVector *BreakItUp(const G4Fragment &theNucleus);


private:

  // This finds temperature of breaking channel.
  G4bool FindTemperatureOfBreakingChannel(const G4Fragment & theFragment, 
					  const G4double & Multiplicity,
					  G4double & Temperature,
					  G4double & EnergyCol);

  // Calculate asymptotic fragments momenta 
  void CoulombImpulse(const G4Fragment & theFragment,
		      const G4int & NumberOfChargedFragments,
		      const G4int & Multiplicity,
		      const G4double & Temperature,
		      const G4double & CoulombEnergy,
		      G4ThreeVector * MomentumOfFragments);



  // Randomly samples fragments positions inside prolongated ellipsoid
  void Place(const G4Fragment & theFragment,
	     const G4int & Multiplicity,
	     G4ThreeVector * Position);


  // This method will find a solution of Newton's equation of motion
  // for fragments in the self-consistent time-dependent Coulomb field
  void SolveEqOfMotion(G4ThreeVector * InitialPos,
		       G4ThreeVector * InitialVel,
		       G4ThreeVector * FinalVel,
		       const G4int & Multiplicity,
		       const G4double & CoulombEnergy,
		       const G4double & KineticEnergy);

  // Calculates fragments momentum components at the breakup instant.
  // Fragment kinetic energies will be calculated according to the
  // Boltzamann distribution at given temperature.
  void CalculateFragmentsMomentum(const G4int & INET,
				  const G4int & NFrags,
				  const G4double & T,
				  const G4double & TotKineticE,
				  G4ThreeVector * Momentum);

  // Rotates a 3-vector P to close momentum triangle P + A + B = 0
  G4ThreeVector Rotor(const G4ThreeVector & P,
		      const G4ThreeVector & A,
		      const G4ThreeVector & B);

  G4double CalculateFragmentExcitationEnergy(const G4int & index, const G4double & T);

  // Samples a isotropic random vectorwith a magnitud given by Magnitude.
  // By default Magnitude = 1
  G4ThreeVector IsotropicVector(const G4double Magnitude = 1.0);


private:

  //  G4StatMFMicrocanonical * theMicrocanonicalSim;
  //  G4StatMFMacrocanonical * theMacrocanonicalSim;

  G4VStatMFCanonical * theSim;



};

#endif
