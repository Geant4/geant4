// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PreCompoundModel.hh,v 1.6 1999/12/15 14:52:38 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// by V. Lara

#ifndef G4PreCompoundModel_h
#define G4PreCompoundModel_h 1

#include "G4VPreCompoundModel.hh"
#include "G4PreCompoundNeutron.hh"
#include "G4PreCompoundProton.hh"
#include "G4PreCompoundDeuteron.hh"
#include "G4PreCompoundTriton.hh"
#include "G4PreCompoundHe3.hh"
#include "G4PreCompoundAlpha.hh"
#include "G4PreCompoundTransitions.hh"
#include "G4LorentzVector.hh"

#include "G4NucleiProperties.hh"
#include "G4Proton.hh"
#include "G4VPreCompoundFragment.hh"
#include "G4PreCompoundParameters.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"



class G4Fragment;



class G4PreCompoundModel : public G4VPreCompoundModel
{
public:
  
  G4PreCompoundModel(G4ExcitationHandler * const value);

  ~G4PreCompoundModel();

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

  
  G4ParticleChange theResult;
  
  
  
  //  static const G4int NumberOfPossibleFragments = 6;
  enum {NumberOfPossibleFragments = 6};
  
  // The possible emitted fragments 
  G4RWTPtrOrderedVector<G4VPreCompoundFragment> theChannels;



  G4ThreeVector IsotropicRandom3Vector(G4double Magnitude = 1.0) const;


  void PerformEquilibriumEmission(const G4Fragment & aFragment, 
				  G4ReactionProductVector * theResult) const;

	G4ParticleMomentum RotateMomentum(G4ParticleMomentum Pa, G4ParticleMomentum V,
				    G4ParticleMomentum P) const;

};


#endif


