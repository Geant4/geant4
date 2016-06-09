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
#ifndef G4ELEMENTARY_PARTICLE_COLLIDER_HH
#define G4ELEMENTARY_PARTICLE_COLLIDER_HH

#include "G4Collider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4CascadSpecialFunctions.hh"
#include "G4LorentzConvertor.hh"

#include "G4CascadeKplusPChannel.hh"
#include "G4CascadeKplusNChannel.hh"
#include "G4CascadeKzeroPChannel.hh"
#include "G4CascadeKzeroNChannel.hh"
#include "G4CascadeKminusPChannel.hh"
#include "G4CascadeKminusNChannel.hh"
#include "G4CascadeKzeroBarPChannel.hh"
#include "G4CascadeKzeroBarNChannel.hh"
#include "G4CascadeLambdaPChannel.hh"
#include "G4CascadeLambdaNChannel.hh"
#include "G4CascadeSigmaPlusPChannel.hh"
#include "G4CascadeSigmaPlusNChannel.hh"
#include "G4CascadeSigmaZeroPChannel.hh"
#include "G4CascadeSigmaZeroNChannel.hh"
#include "G4CascadeSigmaMinusPChannel.hh"
#include "G4CascadeSigmaMinusNChannel.hh"
#include "G4CascadeXiZeroPChannel.hh"
#include "G4CascadeXiZeroNChannel.hh"
#include "G4CascadeXiMinusPChannel.hh"
#include "G4CascadeXiMinusNChannel.hh"

using namespace G4InuclSpecialFunctions;
using namespace G4CascadSpecialFunctions;

class G4ElementaryParticleCollider : public G4Collider {

public:

  G4ElementaryParticleCollider();

  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				    G4InuclParticle* target);

private:

  G4CascadeKplusPChannel kpp;
  G4CascadeKplusNChannel kpn;
  G4CascadeKzeroPChannel k0p;
  G4CascadeKzeroNChannel k0n;
  G4CascadeKminusPChannel kmp;
  G4CascadeKminusNChannel kmn;
  G4CascadeKzeroBarPChannel k0bp;
  G4CascadeKzeroBarNChannel k0bn;
  G4CascadeLambdaPChannel lp;
  G4CascadeLambdaNChannel ln;
  G4CascadeSigmaPlusPChannel spp;
  G4CascadeSigmaPlusNChannel spn;
  G4CascadeSigmaZeroPChannel s0p;
  G4CascadeSigmaZeroNChannel s0n;
  G4CascadeSigmaMinusPChannel smp;
  G4CascadeSigmaMinusNChannel smn;
  G4CascadeXiZeroPChannel x0p;
  G4CascadeXiZeroNChannel x0n;
  G4CascadeXiMinusPChannel xmp;
  G4CascadeXiMinusNChannel xmn;

  G4int verboseLevel;
  G4int generateMultiplicity(G4int is, 
			     G4double ekin) const;
      
  std::vector<G4InuclElementaryParticle> generateSCMfinalState(G4double ekin, 
							  G4double etot_scm, G4double pscm,	     
							  G4InuclElementaryParticle* particle1,
							  G4InuclElementaryParticle* particle2, 
							  G4LorentzConvertor* toSCM) const; 

  std::vector<G4double> generateMomModules(const std::vector<G4int>& kinds, 
				      G4int mult,
				      G4int is, 
				      G4double ekin, 
				      G4double etot_cm) const; 
      
  G4bool reChargering(G4double ekin, 
		      G4int is) const;


  std::vector<G4double> particleSCMmomentumFor2to2(G4int is, 
			             G4int kw, 
				     G4double ekin,
				     G4double pscm) const; 
    
  G4int getElasticCase(G4int is, 
		       G4int kw, 
		       G4double ekin) const;

  std::vector<G4int> generateOutgoingKindsFor2toMany(G4int is, 
						G4int mult, 
						G4double ekin) const;


  std::vector<G4int> generateStrangeChannelPartTypes(G4int is, 
						G4int mult, 
						G4double ekin) const;


  G4double getMomModuleFor2toMany(G4int is, 
				  G4int mult, 
				  G4int knd, 
				  G4double ekin) const; 

  G4bool satisfyTriangle(const std::vector<G4double>& modules) const; 
	
  std::vector<G4double> particleSCMmomentumFor2to3(G4int is, 
					      G4int knd, 
					      G4double ekin, 
					      G4double pmod) const; 
	
  G4int getIL(G4int is, 
	      G4int mult) const; 

  std::pair<G4double, G4double> adjustIntervalForElastic(G4double ekin, 
						    G4double ak, 
						    G4double ae,
						    G4int k, 
						    G4int l, 
						    const std::vector<G4double>& ssv, 
						    G4double st) const;
 
  std::vector<G4InuclElementaryParticle> 
  generateSCMpionAbsorption(G4double etot_scm,
			    G4InuclElementaryParticle* particle1,
			    G4InuclElementaryParticle* particle2) const; 
    
};        

#endif // G4ELEMENTARY_PARTICLE_COLLIDER_HH


