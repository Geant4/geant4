// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpBoundaryProcess.hh,v 1.5 1999-12-15 14:53:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// Optical Photon Boundary Process Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpBoundaryProcess.hh
// Description: Discrete Process -- reflection/refraction at
//                                  optical interfaces
// Version:     1.1
// Created:     1997-06-18
// Modified:    1999-10-29 add method and class descriptors
//              1999-10-10 - Fill NewMomentum/NewPolarization in 
//                           DoAbsorption. These members need to be
//                           filled since DoIt calls 
//                           aParticleChange.SetMomentumChange etc.
//                           upon return (thanks to: Clark McGrew)
//
// Author:      Peter Gumplinger
//              adopted from work by Werner Keil - April 2/96
// mail:        gum@triumf.ca
//
// CVS version tag: 
////////////////////////////////////////////////////////////////////////

#ifndef G4OpBoundaryProcess_h
#define G4OpBoundaryProcess_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4OpticalPhoton.hh"
#include "G4TransportationManager.hh"

// Class Description:
// Discrete Process -- reflection/refraction at optical interfaces.
// Class inherits publicly from G4VDiscreteProcess.                  
// Class Description - End:             

/////////////////////
// Class Definition
/////////////////////

class G4OpBoundaryProcess : public G4VDiscreteProcess 
{

enum G4OpBoundaryProcessStatus {  Undefined,
				  FresnelRefraction, FresnelReflection, 
				  TotalInternalReflection,
				  LambertianReflection, LobeReflection,
				  SpikeReflection, BackScattering,
				  Absorption, Detection };

private:

        //////////////
        // Operators
        //////////////

        // G4OpBoundaryProcess& operator=(const G4OpBoundaryProcess &right);

        // G4OpBoundaryProcess(const G4OpBoundaryProcess &right);

public: // Without description

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        G4OpBoundaryProcess(const G4String& processName = "Boundary");

	~G4OpBoundaryProcess();

	////////////
	// Methods
        ////////////

public: // With description

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an optical photon.

	G4double GetMeanFreePath(const G4Track& ,
				 G4double ,
				 G4ForceCondition* condition);
        // Returns infinity; i. e. the process does not limit the step,
        // but sets the 'Forced' condition for the DoIt to be invoked at
        // every step. However, only at a boundary will any action be
        // taken. 

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				       const G4Step&  aStep);
        // This is the method implementing boundary processes.

	G4OpticalSurfaceModel GetModel() const;
        // Returns the optical surface mode.

	void           SetModel(G4OpticalSurfaceModel model);
	// Set the optical surface model to be followed 
        // (glisur || unified). 

private:

	void G4Swap(G4double* a, G4double* b) const;

	void G4Swap(G4Material* a, G4Material* b) const;

	void G4VectorSwap(G4ThreeVector* vec1, G4ThreeVector* vec2) const;

	G4bool G4BooleanRand(const G4double prob) const;

	G4ThreeVector G4IsotropicRand() const;

	G4ThreeVector G4LambertianRand(const G4ThreeVector& normal);

	G4ThreeVector G4PlaneVectorRand(const G4ThreeVector& normal) const; 

	G4ThreeVector GetFacetNormal(const G4ThreeVector& Momentum,
				     const G4ThreeVector&  Normal) const;

	void DielectricMetal();
	void DielectricDielectric();

	void ChooseReflection();
	void DoAbsorption();
	void DoReflection();

private:

	G4double thePhotonMomentum;

	G4ThreeVector OldMomentum;
	G4ThreeVector OldPolarization;

	G4ThreeVector NewMomentum;
	G4ThreeVector NewPolarization;

	G4ThreeVector theGlobalNormal;
	G4ThreeVector theFacetNormal;

	G4Material* Material1;
	G4Material* Material2;

	G4OpticalSurface* OpticalSurface;

	G4double Rindex1;
	G4double Rindex2;

	G4double cost1, cost2, sint1, sint2;

	G4OpBoundaryProcessStatus theStatus;

	G4OpticalSurfaceModel theModel;

	G4OpticalSurfaceFinish theFinish;

	G4double theReflectivity;
	G4double theEfficiency;
	G4double prob_sl, prob_ss, prob_bs;

};

////////////////////
// Inline methods
////////////////////

inline
void G4OpBoundaryProcess::G4Swap(G4double* a, G4double* b) const
{
	// swaps the contents of the objects pointed 
	// to by 'a' and 'b'!

  G4double temp;

  temp = *a;
  *a = *b;
  *b = temp;
}

inline
void G4OpBoundaryProcess::G4Swap(G4Material* a, G4Material* b) const
{
	// ONLY swaps the pointers; i.e. what used to be pointed
	// to by 'a' is now pointed to by 'b' and vice versa!

   G4Material* temp = a;

   a = b;
   b = temp;
}

inline
void G4OpBoundaryProcess::G4VectorSwap(G4ThreeVector* vec1,
				       G4ThreeVector* vec2) const
{
        // swaps the contents of the objects pointed
        // to by 'vec1' and 'vec2'!

  G4ThreeVector temp;

  temp = *vec1;
  *vec1 = *vec2;
  *vec2 = temp;
}

inline
G4bool G4OpBoundaryProcess::G4BooleanRand(const G4double prob) const
{
  /* Returns a random boolean variable with the specified probability */

  return (G4UniformRand() < prob);
}

inline
G4ThreeVector G4OpBoundaryProcess::G4IsotropicRand() const
{
  /* Returns a random isotropic unit vector. */

  G4ThreeVector vect;
  G4double len2;

  do {

    vect.setX(G4UniformRand() - 0.5);
    vect.setY(G4UniformRand() - 0.5);
    vect.setZ(G4UniformRand() - 0.5);

    len2 = vect.mag2();

  } while (len2 < 0.01 || len2 > 0.25);

  return vect.unit();
}

inline
G4ThreeVector G4OpBoundaryProcess::
	      G4LambertianRand(const G4ThreeVector& normal)
{
  /* Returns a random lambertian unit vector. */

  G4ThreeVector vect;
  G4double ndotv;

  do {
    vect = G4IsotropicRand();

    ndotv = normal * vect;

    if (ndotv < 0.0) {
      vect = -vect;
      ndotv = -ndotv;
    }

  } while (!G4BooleanRand(ndotv));
  return vect;
}

inline
G4ThreeVector G4OpBoundaryProcess::
	      G4PlaneVectorRand(const G4ThreeVector& normal) const

  /* This function chooses a random vector within a plane given
     by the unit normal */
{
  G4ThreeVector vec1 = normal.orthogonal();

  G4ThreeVector vec2 = vec1.cross(normal);

  G4double cost = 2.*G4UniformRand() - 1.0;
  G4double sint = sqrt(1.0 - cost * cost);

  return cost * vec1 + sint * vec2;
}

inline
G4bool G4OpBoundaryProcess::IsApplicable(const G4ParticleDefinition& 
					               aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

inline
G4OpticalSurfaceModel G4OpBoundaryProcess::GetModel() const
{
   return theModel;
}

inline
void G4OpBoundaryProcess::SetModel(G4OpticalSurfaceModel model)
{
   theModel = model;
}

inline
void G4OpBoundaryProcess::ChooseReflection()
{
                 G4double rand = G4UniformRand();
                 if ( rand >= 0.0 && rand < prob_ss ) {
                    theStatus = SpikeReflection;
                    theFacetNormal = theGlobalNormal;
                 }
                 else if ( rand >= prob_ss &&
                           rand <= prob_ss+prob_sl) {
                    theStatus = LobeReflection;
                 }
                 else if ( rand > prob_ss+prob_sl &&
                           rand < prob_ss+prob_sl+prob_bs ) {
                    theStatus = BackScattering;
                 }
                 else {
                    theStatus = LambertianReflection;
                 }
}

inline
void G4OpBoundaryProcess::DoAbsorption()
{
              theStatus = Absorption;

              if ( G4BooleanRand(theEfficiency) ) {
                 // EnergyDeposited =/= 0 means: photon has been detected
                 theStatus = Detection;
                 aParticleChange.SetLocalEnergyDeposit(thePhotonMomentum);
              }
              else {
                 aParticleChange.SetLocalEnergyDeposit(0.0);
              }

              NewMomentum = OldMomentum;
              NewPolarization = OldPolarization;

//              aParticleChange.SetEnergyChange(0.0);
              aParticleChange.SetStatusChange(fStopAndKill);
}

inline
void G4OpBoundaryProcess::DoReflection()
{
        if ( theStatus == LambertianReflection ) {

          NewMomentum = G4LambertianRand(theGlobalNormal);
          theFacetNormal = (NewMomentum - OldMomentum).unit();

        }
        else if ( theFinish == ground ) {

	  theStatus = LobeReflection;
          theFacetNormal = GetFacetNormal(OldMomentum,theGlobalNormal);
          G4double PdotN = OldMomentum * theFacetNormal;
          NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;

        }
        else {

          theStatus = SpikeReflection;
          theFacetNormal = theGlobalNormal;
          G4double PdotN = OldMomentum * theFacetNormal;
          NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;

        }
        G4double EdotN = OldPolarization * theFacetNormal;
        NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
}

#endif /* G4OpBoundaryProcess_h */
