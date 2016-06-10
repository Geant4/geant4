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
// $Id: G4CascadeColliderBase.hh 71942 2013-06-28 19:08:11Z mkelsey $
//
// 20100714  M. Kelsey -- Move functionality from G4VCascadeCollider, and
//		provide conservation-checking here, with wrapper function
//		and control flag.
// 20100720  M. Kelsey -- Change G4CascadeCheckBalance to pointer member
// 20100925  M. Kelsey -- Add explosion(A,Z,Eex) and explosion(G4Fragment)
//		interfaces
// 20110225  M. Kelsey -- Add setVerboseLevel(), calls through to members
// 20110304  M. Kelsey -- Add dummy rescatter() interface here, to enforce
//		consistency in subclass colliders.
// 20110321  M. Kelsey -- Hide names of arguments to rescatter(), to avoid
//		compiler warnings on some GCC versions.
// 20130620  Address Coverity complaint about missing copy actions
// 20130621  Move doConservationChecks to G4CascadeParameters; change
//		explosion to use reference, add validateOutput() w/G4Fragment
// 20130622  Move fragment-handling functions to G4CascadeDeexciteBase
// 20140930  Change name from "const char*" to "const G4String"

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
class G4KineticTrackVector;
class G4V3DNucleus;


class G4CascadeColliderBase : public G4VCascadeCollider {
public:
  G4CascadeColliderBase(const G4String& name, G4int verbose=0);
  virtual ~G4CascadeColliderBase();

  // For use with top-level Propagate to preload a set of secondaries
  virtual void rescatter(G4InuclParticle* /*bullet*/,
			 G4KineticTrackVector* /*theSecondaries*/,
			 G4V3DNucleus* /*theNucleus*/,
			 G4CollisionOutput& /*globalOutput*/) { ; }

  virtual void setVerboseLevel(G4int verbose=0);

protected:
  G4InteractionCase interCase;		// Determine bullet vs. target

  // Decide whether to use G4ElementaryParticleCollider or not
  virtual G4bool useEPCollider(G4InuclParticle* bullet, 
			       G4InuclParticle* target) const;

  // Decide whether to use G4IntraNuclearCascader or not
  virtual G4bool inelasticInteractionPossible(G4InuclParticle* bullet,
					      G4InuclParticle* target, 
					      G4double ekin) const;

  // ==> Provide same interfaces as G4CascadeCheckBalance itself
  G4CascadeCheckBalance* balance;

  // Validate output for energy, momentum conservation, etc.
  virtual G4bool validateOutput(G4InuclParticle* bullet,
				G4InuclParticle* target,
				G4CollisionOutput& output);

  // This is for use after de-excitation
  virtual G4bool validateOutput(const G4Fragment& fragment,
				G4CollisionOutput& output);

  // This is for use with G4EPCollider
  virtual G4bool validateOutput(G4InuclParticle* bullet,
				G4InuclParticle* target,
		const std::vector<G4InuclElementaryParticle>& particles);

private:
  // Copying of modules is forbidden
  G4CascadeColliderBase(const G4CascadeColliderBase&);
  G4CascadeColliderBase& operator=(const G4CascadeColliderBase&);
};        

#endif	/* G4CASCADE_COLLIDER_BASE_HH */
