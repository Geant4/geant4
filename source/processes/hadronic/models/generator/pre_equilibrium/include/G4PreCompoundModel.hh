// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

// Class Description
// Model implementation for pre-equilibrium decay models in geant4. 
// To be used in your physics list, in case you neeed this kind of physics.
// Can be used as a stand-allone model, but also in conjunction with an intra-nuclear
// transport, or any of the string-parton models.
// Class Description - End

#ifndef G4PreCompoundModel_h
#define G4PreCompoundModel_h 1

#include "G4VPreCompoundModel.hh"
#include "G4PreCompoundTransitions.hh"
#include "G4PreCompoundEmission.hh"
#include "G4LorentzVector.hh"

#include "G4NucleiProperties.hh"
#include "G4Proton.hh"
#include "G4VPreCompoundFragment.hh"
#include "G4PreCompoundParameters.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"


class G4PreCompoundModel : public G4VPreCompoundModel
{
public:
  
	G4PreCompoundModel(G4ExcitationHandler * const value) : 
		G4VPreCompoundModel(value) {};

	~G4PreCompoundModel() {};

private:
	G4PreCompoundModel() {};
  
	G4PreCompoundModel(const G4PreCompoundModel &right) {};
  
	const G4PreCompoundModel& operator=(const G4PreCompoundModel &right);
	G4bool operator==(const G4PreCompoundModel &right) const;
	G4bool operator!=(const G4PreCompoundModel &right) const;

public:
  G4VParticleChange * ApplyYourself(const G4Track & thePrimary, G4Nucleus & theNucleus);
  
  G4ReactionProductVector* DeExcite(const G4Fragment& aFragment) const;

private:  

  void PerformEquilibriumEmission(const G4Fragment & aFragment, 
				  G4ReactionProductVector * theResult) const;
				  
	G4ParticleChange theResult;
};


#endif


