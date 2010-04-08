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
// $Id: G4ElementaryParticleCollider.hh,v 1.26 2010-04-08 15:48:00 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100315  M. Kelsey -- Remove "using" directive and unnecessary #includes.
// 20100407  M. Kelsey -- Eliminate return-by-value std::vector<> by creating
//		three data buffers for particles, momenta, and particle types.
//		The following functions now return void and are non-const:
//		  ::generateSCMfinalState()
//		  ::generateMomModules() (also remove input vector)
//		  ::generateStrangeChannelPartTypes()
//		  ::generateSCMpionAbsorption()

#ifndef G4ELEMENTARY_PARTICLE_COLLIDER_HH
#define G4ELEMENTARY_PARTICLE_COLLIDER_HH

#include "G4CollisionOutput.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4LorentzVector.hh"
#include "G4NucleonSampler.hh"
#include "G4PionSampler.hh"
#include <vector>

class G4LorentzConvertor;

class G4ElementaryParticleCollider {

public:

  G4ElementaryParticleCollider();

  G4CollisionOutput collide(G4InuclParticle* bullet,
			    G4InuclParticle* target);

private:

  G4int verboseLevel;

  void initializeArrays();

  G4int generateMultiplicity(G4int is, G4double ekin) const;

  void collide(G4InuclElementaryParticle* bullet,
	       G4InuclElementaryParticle* target,
	       G4CollisionOutput& output);

      
  void generateSCMfinalState(G4double ekin, G4double etot_scm, G4double pscm,
			     G4InuclElementaryParticle* particle1,
			     G4InuclElementaryParticle* particle2, 
			     G4LorentzConvertor* toSCM); 

  void generateSCMpionAbsorption(G4double etot_scm,
				 G4InuclElementaryParticle* particle1,
				 G4InuclElementaryParticle* particle2); 

  void generateMomModules(G4int mult, G4int is, G4double ekin,
			  G4double etot_cm); 


  G4LorentzVector
  particleSCMmomentumFor2to2(G4int is, G4int kw, G4double ekin,
			     G4double pscm) const; 


  G4int getElasticCase(G4int is, G4int kw, G4double ekin) const;


  void generateStrangeChannelPartTypes(G4int is, G4int mult, 
				       G4double ekin);


  G4double getMomModuleFor2toMany(G4int is, G4int mult, G4int knd, 
				  G4double ekin) const; 


  G4bool satisfyTriangle(const std::vector<G4double>& modules) const; 
	
  G4LorentzVector
  particleSCMmomentumFor2to3(G4int is, G4int knd, G4double ekin, 
			     G4double pmod) const; 


  std::pair<G4double, G4double> 
  adjustIntervalForElastic(G4double ekin, G4double ak, G4double ae,
                           G4int k, G4int l, const std::vector<G4double>& ssv, 
			   G4double st) const;
 
  G4NucleonSampler nucSampler;
  G4PionSampler piSampler;

  // Internal buffers for lists of secondaries
  std::vector<G4InuclElementaryParticle> particles;
  std::vector<G4double> modules;
  std::vector<G4int> particle_kinds;

  // Parameter arrays

  G4double rmn[14][10][2];    
  G4double ang[4][4][13];
  G4double abn[4][4][4];

};

#endif // G4ELEMENTARY_PARTICLE_COLLIDER_HH


