// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFSFissionFS.hh,v 1.2 1999-06-29 18:43:57 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFSFissionFS_h
#define G4NeutronHPFSFissionFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4DynamicParticleVector.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPNeutronYield.hh"
#include "G4NeutronHPVector.hh"
#include "G4NeutronHPFissionERelease.hh"
#include "G4NeutronHPEnergyDistribution.hh"
#include "G4NeutronHPPhotonDist.hh"
#include "G4NeutronHPAngular.hh"

class G4NeutronHPFSFissionFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPFSFissionFS();
  ~G4NeutronHPFSFissionFS();
  
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  
  G4DynamicParticleVector * ApplyYourself(G4int Prompt, G4int delayed, G4double *decayconst);
  
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPFSFissionFS * theNew = new G4NeutronHPFSFissionFS;
   return theNew;
  }
  
  inline G4double GetMass(){ return targetMass; }
  
  void SampleNeutronMult(G4int&all, 
	  		 G4int&Prompt, 
			 G4int&delayed, 
			 G4double energy,
			 G4int off);
						 
  inline void SetNeutron(const G4ReactionProduct & aNeutron)
                        { 
                          theNeutron = aNeutron;
                          theNeutronAngularDis.SetNeutron(aNeutron);
                        }
  
  inline void SetTarget(const G4ReactionProduct & aTarget)
                        { 
                          theTarget = aTarget; 
                          theNeutronAngularDis.SetTarget(aTarget);
                        }
    
  G4DynamicParticleVector * GetPhotons();
  
  inline G4NeutronHPFissionERelease * GetEnergyRelease()
  {
    return &theEnergyRelease;
  }
  
  private:

  G4ParticleChange * ApplyYourself(const G4Track & theTrack) { return NULL; }
  
  G4double targetMass;
  
  G4NeutronHPNeutronYield theFinalStateNeutrons;
  G4NeutronHPEnergyDistribution thePromptNeutronEnDis;
  G4NeutronHPEnergyDistribution theDelayedNeutronEnDis;
  G4NeutronHPAngular theNeutronAngularDis;
  
  G4NeutronHPPhotonDist theFinalStatePhotons;
  G4NeutronHPFissionERelease theEnergyRelease;
  
  G4ReactionProduct theNeutron;
  G4ReactionProduct theTarget;
  
  private:
  
  G4NeutronHPNames theNames;
  
};
#endif
