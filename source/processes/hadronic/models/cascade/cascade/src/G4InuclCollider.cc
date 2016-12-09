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
// $Id: G4InuclCollider.cc 71769 2013-06-21 21:23:50Z mkelsey $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100309  M. Kelsey -- Eliminate some unnecessary std::pow()
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100418  M. Kelsey -- Move lab-frame transformation code to G4CollisonOutput
// 20100429  M. Kelsey -- Change "photon()" to "isPhoton()"
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members, consolidate code
// 20100620  M. Kelsey -- Reorganize top level if-blocks to reduce nesting,
//		use new four-vector conservation check.
// 20100701  M. Kelsey -- Bug fix energy-conservation after equilibrium evap,
//		pass verbosity through to G4CollisionOutput
// 20100714  M. Kelsey -- Move conservation checking to base class, report
//		number of iterations at end
// 20100715  M. Kelsey -- Remove all the "Before xxx" and "After xxx"
//		conservation checks, as the colliders now all do so.  Move
//		local buffers outside while() loop, use new "combined add()"
//		interface for copying local buffers to global.
// 20100716  M. Kelsey -- Drop G4IntraNucleiCascader::setInteractionCase()
// 20100720  M. Kelsey -- Make all the collders pointer members (to reducde
//		external compile dependences).
// 20100915  M. Kelsey -- Move post-cascade colliders to G4CascadeDeexcitation,
//		simplify operational code somewhat
// 20100922  M. Kelsey -- Add functions to select de-excitation method;
//		default is G4CascadeDeexcitation (i.e., built-in modules)
// 20100924  M. Kelsey -- Migrate to integer A and Z
// 20101019  M. Kelsey -- CoVerity report: check dynamic_cast<> for null
// 20110224  M. Kelsey -- Add ::rescatter() function which takes a list of
//		pre-existing secondaries as input.  Add setVerboseLevel().
// 20110301  M. Kelsey -- Pass verbosity to new or changed de-excitation
// 20110304  M. Kelsey -- Modify rescatter to use original Propagate() input
// 20110308  M. Kelsey -- Separate de-excitation block from collide(); check
//		for single-nucleon "fragment", rather than for null fragment
// 20110413  M. Kelsey -- Modify diagnostic messages in ::rescatter() to be
//		equivalent to those from ::collide().
// 20111003  M. Kelsey -- Prepare for gamma-N interactions by checking for
//		final-state tables instead of particle "isPhoton()"
// 20130621  M. Kelsey -- Pass G4Fragment to de-excitation modules directly
// 20140929  M. Kelsey -- Make PreCompound the default de-excitation
// 20141111  M. Kelsey -- Revert default use of PreCompound; replace
//		G4Fragment::GetA() call with GetA_asInt().
// 20150205  M. Kelsey -- New photonuclearOkay() filter to reject events
//		around giant dipole resonance with no hadronic secondaries.
//		Addresses bug #1680.
// 20150220  M. Kelsey -- Improve photonuclearOkay() filter by just checking
//		final-state nucleus vs. target, rather than all secondaries.
// 20150608  M. Kelsey -- Label all while loops as terminating.

#include "G4InuclCollider.hh"
#include "G4CascadeChannelTables.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CascadeDeexcitation.hh"
#include "G4CollisionOutput.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4LorentzConvertor.hh"
#include "G4PreCompoundDeexcitation.hh"


G4InuclCollider::G4InuclCollider()
  : G4CascadeColliderBase("G4InuclCollider"),
    theElementaryParticleCollider(new G4ElementaryParticleCollider),
    theIntraNucleiCascader(new G4IntraNucleiCascader),
    theDeexcitation(new G4PreCompoundDeexcitation) {}

G4InuclCollider::~G4InuclCollider() {
  delete theElementaryParticleCollider;
  delete theIntraNucleiCascader;
  delete theDeexcitation;
}


// Set verbosity and pass on to member objects
void G4InuclCollider::setVerboseLevel(G4int verbose) {
  G4CascadeColliderBase::setVerboseLevel(verbose);

  theElementaryParticleCollider->setVerboseLevel(verboseLevel);
  theIntraNucleiCascader->setVerboseLevel(verboseLevel);
  theDeexcitation->setVerboseLevel(verboseLevel);

  output.setVerboseLevel(verboseLevel);
  DEXoutput.setVerboseLevel(verboseLevel);
}


// Select post-cascade processing (default will be CascadeDeexcitation)

void G4InuclCollider::useCascadeDeexcitation() {
  delete theDeexcitation;
  theDeexcitation = new G4CascadeDeexcitation;
  theDeexcitation->setVerboseLevel(verboseLevel);
}

void G4InuclCollider::usePreCompoundDeexcitation() {
  delete theDeexcitation;
  theDeexcitation = new G4PreCompoundDeexcitation;
  theDeexcitation->setVerboseLevel(verboseLevel);
}


// Main action

void G4InuclCollider::collide(G4InuclParticle* bullet, G4InuclParticle* target,
			      G4CollisionOutput& globalOutput) {
  if (verboseLevel) G4cout << " >>> G4InuclCollider::collide" << G4endl;

  const G4int itry_max = 100;

  // Particle-on-particle collision; no nucleus involved
  if (useEPCollider(bullet,target)) {
    if (verboseLevel > 2)
      G4cout << " InuclCollider -> particle on particle collision" << G4endl;
 
    theElementaryParticleCollider->collide(bullet, target, globalOutput);
    return;
  }
  
  interCase.set(bullet,target);		// Classify collision type
  if (verboseLevel > 2) {
    G4cout << " InuclCollider -> inter case " << interCase.code() << G4endl;
  }

  if (!interCase.valid()) {
    if (verboseLevel > 1)
      G4cerr << " InuclCollider -> no collision possible " << G4endl;

    globalOutput.trivialise(bullet, target);
    return;
  }

  // Target must be a nucleus
  G4InuclNuclei* ntarget = dynamic_cast<G4InuclNuclei*>(interCase.getTarget());
  if (!ntarget) {
    G4cerr << " InuclCollider -> ERROR target is not a nucleus " << G4endl;

    globalOutput.trivialise(bullet, target);
    return;
  }

  G4int btype = 0;
  G4int ab = 0;
  G4int zb = 0;
  
  if (interCase.hadNucleus()) { 	// particle with nuclei
    G4InuclElementaryParticle* pbullet = 
      dynamic_cast<G4InuclElementaryParticle*>(interCase.getBullet());

    if (!pbullet) {
      G4cerr << " InuclCollider -> ERROR bullet is not a hadron " << G4endl;
      globalOutput.trivialise(bullet, target);
      return;
    }

    if (!G4CascadeChannelTables::GetTable(pbullet->type())) {
      G4cerr << " InuclCollider -> ERROR can not collide with "
	     << pbullet->getDefinition()->GetParticleName() << G4endl;
      globalOutput.trivialise(bullet, target);
      return;
    }

    btype = pbullet->type();
  } else { 				// nuclei with nuclei
    G4InuclNuclei* nbullet = 
      dynamic_cast<G4InuclNuclei*>(interCase.getBullet());
    if (!nbullet) {
      G4cerr << " InuclCollider -> ERROR bullet is not a nucleus " << G4endl;
      globalOutput.trivialise(bullet, target);
      return;
    }
    
    ab = nbullet->getA();
    zb = nbullet->getZ();
  }

  G4LorentzConvertor convertToTargetRestFrame(bullet, ntarget);
  G4double ekin = convertToTargetRestFrame.getKinEnergyInTheTRS();
  
  if (verboseLevel > 3) G4cout << " ekin in trs " << ekin << G4endl;

  if (!inelasticInteractionPossible(bullet, target, ekin)) {
    if (verboseLevel > 3)
      G4cout << " InuclCollider -> inelastic interaction is impossible\n"
	     << " due to the coulomb barirer " << G4endl;

    globalOutput.trivialise(bullet, target);
    return;
  }

  // Generate interaction secondaries in rest frame of target nucleus
  convertToTargetRestFrame.toTheTargetRestFrame();
  if (verboseLevel > 3) {
    G4cout << " degenerated? " << convertToTargetRestFrame.trivial()
	   << G4endl;
  }
  
  G4LorentzVector bmom;			// Bullet is along local Z
  bmom.setZ(convertToTargetRestFrame.getTRSMomentum());

  // Need to make copy of bullet with momentum realigned
  G4InuclParticle* zbullet = 0;
  if (interCase.hadNucleus())
    zbullet = new G4InuclElementaryParticle(bmom, btype);
  else
    zbullet = new G4InuclNuclei(bmom, ab, zb);

  G4int itry = 0;
  while (itry < itry_max) {	/* Loop checking 08.06.2015 MHK */
    itry++;
    if (verboseLevel > 2)
      G4cout << " InuclCollider itry " << itry << G4endl;

    globalOutput.reset();		// Clear buffers for this attempt
    output.reset();

    theIntraNucleiCascader->collide(zbullet, target, output);
    
    if (verboseLevel > 1) G4cout << " After Cascade " << G4endl;

    deexcite(output.getRecoilFragment(), output);
    output.removeRecoilFragment();

    //*** TEMPORARY, USE ENVVAR TO ENABLE/DISABLE THIS TEST ***
    if (getenv("G4CASCADE_CHECK_PHOTONUCLEAR"))
      if (!photonuclearOkay(output)) continue;

    if (verboseLevel > 2)
      G4cout << " itry " << itry << " finished, moving to lab frame" << G4endl;

    // convert to the LAB frame and add to final result
    output.boostToLabFrame(convertToTargetRestFrame);

    globalOutput.add(output);

    // Adjust final state particles to balance momentum and energy
    // FIXME:  This should no longer be necessary!
    globalOutput.setOnShell(bullet, target);
    if (globalOutput.acceptable()) {
      if (verboseLevel) 
	G4cout << " InuclCollider output after trials " << itry << G4endl;
      delete zbullet;
      return;
    } else {
      if (verboseLevel>2)
	G4cerr << " InuclCollider setOnShell failed." << G4endl;
    }
  }	// while (itry < itry_max)
  
  if (verboseLevel) {
    G4cout << " InuclCollider -> can not generate acceptable inter. after " 
	   << itry_max << " attempts " << G4endl;
  }
  
  globalOutput.trivialise(bullet, target);

  delete zbullet;
  return;
}


// For use with Propagate to preload a set of secondaries

void G4InuclCollider::rescatter(G4InuclParticle* bullet,
				G4KineticTrackVector* theSecondaries,
				G4V3DNucleus* theNucleus,
				G4CollisionOutput& globalOutput) {
  if (verboseLevel) G4cout << " >>> G4InuclCollider::rescatter" << G4endl;

  G4int itry=1;		// For diagnostic post-processing only
  if (verboseLevel > 2) G4cout << " InuclCollider itry " << itry << G4endl;

  globalOutput.reset();		// Clear buffers for this attempt
  output.reset();

  theIntraNucleiCascader->rescatter(bullet, theSecondaries, theNucleus, 
				    output);

  if (verboseLevel > 1) G4cout << " After Rescatter" << G4endl;

  deexcite(output.getRecoilFragment(), output);
  output.removeRecoilFragment();

  globalOutput.add(output);	// Add local results to global output

  if (verboseLevel) 
    G4cout << " InuclCollider output after trials " << itry << G4endl;
}


// De-excite nuclear fragment to ground state

void G4InuclCollider::deexcite(const G4Fragment& fragment,
			       G4CollisionOutput& globalOutput) {
  if (fragment.GetA_asInt() <= 1) return;	// Nothing real to be de-excited

  if (verboseLevel) G4cout << " >>> G4InuclCollider::deexcite" << G4endl;

  const G4int itry_max = 10;		// Maximum number of attempts
  G4int itry = 0;
  do {					/* Loop checking 08.06.2015 MHK */
    if (verboseLevel > 2) G4cout << " deexcite itry " << itry << G4endl;

    DEXoutput.reset();
    theDeexcitation->deExcite(fragment, DEXoutput);
  } while (!validateOutput(fragment, DEXoutput) && (++itry < itry_max));

  // Add de-excitation products to output buffer
  globalOutput.add(DEXoutput);
}


// Looks for non-gamma final state in photonuclear or leptonuclear

G4bool G4InuclCollider::photonuclearOkay(G4CollisionOutput& checkOutput) const {
  if (interCase.twoNuclei()) return true;	// A-A is not photonuclear

  G4InuclElementaryParticle* bullet =
    dynamic_cast<G4InuclElementaryParticle*>(interCase.getBullet());
  if (!bullet || !(bullet->isPhoton() || bullet->isElectron())) return true;

  if (verboseLevel>1)
    G4cout << " >>> G4InuclCollider::photonuclearOkay" << G4endl;

  if (bullet->getKineticEnergy() > 0.050) return true;

  if (verboseLevel>2) {
    if (checkOutput.numberOfOutgoingNuclei() > 0) {
      G4cout << " comparing final nucleus with initial target:\n"
             << checkOutput.getOutgoingNuclei()[0] << G4endl
             << *(interCase.getTarget()) << G4endl;
    } else {
      G4cout << " no final nucleus remains when target was "
             << *(interCase.getTarget()) << G4endl;
    }
  }

  // Hadron production changes target nucleus
  G4double mfinalNuc = 0.0;
  if (checkOutput.numberOfOutgoingNuclei() > 0)
    mfinalNuc = checkOutput.getOutgoingNuclei()[0].getMass();
  G4double mtargetNuc = interCase.getTarget()->getMass();
  if (mfinalNuc != mtargetNuc) return true;	// Mass from G4Ions is fixed
  
  if (verboseLevel>2)
    G4cout << " photonuclear produced only gammas.  Try again." << G4endl;

  return false;		// Final state is entirely de-excitation photons
}
