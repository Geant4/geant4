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
#ifndef G4CASCADE_COLLIDER_BASE_HH
#define G4CASCADE_COLLIDER_BASE_HH
// $Id: G4CascadeColliderBase.hh,v 1.4 2010-12-15 07:39:32 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100714  M. Kelsey -- Move functionality from G4VCascadeCollider, and
//		provide conservation-checking here, with wrapper function
//		and control flag.
// 20100720  M. Kelsey -- Change G4CascadeCheckBalance to pointer member
// 20010925  M. Kelsey -- Add explosion(A,Z,Eex) and explosion(G4Fragment)
//		interfaces

#include "G4VCascadeCollider.hh"

#include "globals.hh"
#include "G4InteractionCase.hh"
#include <vector>

class G4InuclElementaryParticle;
class G4InuclNuclei;
class G4InuclParticle;
class G4CollisionOutput;
class G4CascadeCheckBalance;
class G4Fragment;

class G4CascadeColliderBase : public G4VCascadeCollider {
public:
  G4CascadeColliderBase(const char* name, G4int verbose=0);
  virtual ~G4CascadeColliderBase();

  virtual void setConservationChecks(G4bool doBalance=true) {
    doConservationChecks = doBalance;
  }

protected:
  G4InteractionCase interCase;		// Determine bullet vs. target
  G4bool doConservationChecks;		// Conservation-law validation
  G4CascadeCheckBalance* balance;

  // Decide whether to use G4ElementaryParticleCollider or not
  virtual G4bool useEPCollider(G4InuclParticle* bullet, 
			       G4InuclParticle* target) const;

  // Decide whether to use G4BigBanger or not
  virtual G4bool explosion(G4InuclNuclei* target) const;
  virtual G4bool explosion(G4Fragment* target) const;
  virtual G4bool explosion(G4int A, G4int Z, G4double excitation) const;

  // Decide whether to use G4IntraNuclearCascader or not
  virtual G4bool inelasticInteractionPossible(G4InuclParticle* bullet,
					      G4InuclParticle* target, 
					      G4double ekin) const;

  // ==> Provide same interfaces as G4CascadeCheckBalance itself

  // Validate output for energy, momentum conservation, etc.
  virtual G4bool validateOutput(G4InuclParticle* bullet,
				G4InuclParticle* target,
				G4CollisionOutput& output);

  // This is for use with G4EPCollider and G4BigBanger
  virtual G4bool validateOutput(G4InuclParticle* bullet,
				G4InuclParticle* target,
		const std::vector<G4InuclElementaryParticle>& particles);

  // This is for use with G4Fissioner
  virtual G4bool validateOutput(G4InuclParticle* bullet,
				G4InuclParticle* target,
		const std::vector<G4InuclNuclei>& fragments);
};        

#endif	/* G4CASCADE_COLLIDER_BASE_HH */
