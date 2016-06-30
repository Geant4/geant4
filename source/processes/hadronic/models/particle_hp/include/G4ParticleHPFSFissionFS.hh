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
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPFSFissionFS_h
#define G4ParticleHPFSFissionFS_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ParticleHPNames.hh"
#include "G4ParticleHPParticleYield.hh"
#include "G4ParticleHPVector.hh"
#include "G4ParticleHPFissionERelease.hh"
#include "G4ParticleHPEnergyDistribution.hh"
#include "G4ParticleHPPhotonDist.hh"
#include "G4ParticleHPAngular.hh"
#include "G4Cache.hh"

class G4ParticleHPFSFissionFS : public G4ParticleHPFinalState
{
   struct toBeCached {
      const G4ReactionProduct* theNeutronRP;
      const G4ReactionProduct* theTarget;
      toBeCached() : theNeutronRP(NULL),theTarget(NULL){};
   };

  public:
  
  G4ParticleHPFSFissionFS(){ hasXsec = true; }
  ~G4ParticleHPFSFissionFS(){}
  
  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType, G4ParticleDefinition*);
  
  G4DynamicParticleVector * ApplyYourself(G4int Prompt, G4int delayed, G4double *decayconst);
  
  G4ParticleHPFinalState * New() 
  {
   G4ParticleHPFSFissionFS * theNew = new G4ParticleHPFSFissionFS;
   return theNew;
  }
  
  inline G4double GetMass(){ return theFinalStateNeutrons.GetTargetMass(); }
  
  void SampleNeutronMult(G4int&all, 
	  		 G4int&Prompt, 
			 G4int&delayed, 
			 G4double energy,
			 G4int off);
						 
  inline void SetNeutronRP(const G4ReactionProduct & aNeutron)
                        { 
                          fCache.Get().theNeutronRP = &aNeutron; 
                          theNeutronAngularDis.SetProjectileRP(aNeutron);
                        }
  
  inline void SetTarget(const G4ReactionProduct & aTarget)
                        { 
                          fCache.Get().theTarget = &aTarget; 
                          theNeutronAngularDis.SetTarget(aTarget);
                        }
    
  G4DynamicParticleVector * GetPhotons();
  
  inline G4ParticleHPFissionERelease * GetEnergyRelease()
  {
    return &theEnergyRelease;
  }
  
  private:

  G4HadFinalState * ApplyYourself(const G4HadProjectile & ) { return 0; }  
  //G4double targetMass;
  
  G4ParticleHPParticleYield theFinalStateNeutrons;
  G4ParticleHPEnergyDistribution thePromptNeutronEnDis;
  G4ParticleHPEnergyDistribution theDelayedNeutronEnDis;
  G4ParticleHPAngular theNeutronAngularDis;
  
  G4ParticleHPPhotonDist theFinalStatePhotons;
  G4ParticleHPFissionERelease theEnergyRelease;
  
  //G4ReactionProduct theNeutronRP;
  //G4ReactionProduct theTarget;
   private:
      G4Cache<toBeCached> fCache;
  
  private:
  
  G4ParticleHPNames theNames;
  
};
#endif
