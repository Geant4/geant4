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
// Author:  Michael Kelsey (SLAC)
// Date:    15 April 2013
//
// Description: Subclass of models/util G4VHadDecayAlgorithm which uses
//		old INUCL parametrizations for momentum and angular
//		distributions.
//
// 20130509  BUG FIX:  Two-body momentum vector should be rotated into
//		collision axis; three-body "final" vector needs to be rotated
//		into axis of rest of system.  Tweak some diagnostic messages
//		to match old G4EPCollider version.
// 20130612  BUG FIX:  Create separate FillDirThreeBody(), FillDirManyBody()
//		in order to reporoduce old method: N-body states are filled
//		from first to last, while three-body starts with the last.
// 20130702  M. Kelsey: Copy phase-space algorithm from Kopylov; use if
//		runtime envvar G4CASCADE_USE_PHASESPACE is set
// 20140627  BUG FIX:  Use ".c_str()" in diagnostics to avoid IBM XL error.
// 20150608  M. Kelsey -- Label all while loops as terminating.
// 20150619  M. Kelsey -- Replace std::exp with G4Exp

#include "G4CascadeFinalStateAlgorithm.hh"
#include "G4CascadeParameters.hh"
#include "G4Exp.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4LorentzConvertor.hh"
#include "G4MultiBodyMomentumDist.hh"
#include "G4Pow.hh"
#include "G4TwoBodyAngularDist.hh"
#include "G4VMultiBodyMomDst.hh"
#include "G4VTwoBodyAngDst.hh"
#include "Randomize.hh"
#include <vector>
#include <numeric>
#include <cmath>

using namespace G4InuclSpecialFunctions;


// Cut-offs and iteration limits for generation

const G4double G4CascadeFinalStateAlgorithm::maxCosTheta = 0.9999;
const G4double G4CascadeFinalStateAlgorithm::oneOverE = 0.3678794;   
const G4double G4CascadeFinalStateAlgorithm::small = 1.e-10;
const G4int G4CascadeFinalStateAlgorithm::itry_max = 10;


// Constructor and destructor

G4CascadeFinalStateAlgorithm::G4CascadeFinalStateAlgorithm()
  : G4VHadDecayAlgorithm("G4CascadeFinalStateAlgorithm"),
    momDist(0), angDist(0), multiplicity(0), bullet_ekin(0.) {;}

G4CascadeFinalStateAlgorithm::~G4CascadeFinalStateAlgorithm() {;}

void G4CascadeFinalStateAlgorithm::SetVerboseLevel(G4int verbose) {
  G4VHadDecayAlgorithm::SetVerboseLevel(verbose);
  G4MultiBodyMomentumDist::setVerboseLevel(verbose);
  G4TwoBodyAngularDist::setVerboseLevel(verbose);
  toSCM.setVerbose(verbose);
}


// Select distributions to be used for next interaction

void G4CascadeFinalStateAlgorithm::
Configure(G4InuclElementaryParticle* bullet,
	  G4InuclElementaryParticle* target,
	  const std::vector<G4int>& particle_kinds) {
  if (GetVerboseLevel()>1)
    G4cout << " >>> " << GetName() << "::Configure" << G4endl;

  // Identify initial and final state (if two-body) for algorithm selection
  multiplicity = (G4int)particle_kinds.size();
  G4int is = bullet->type() * target->type();
  G4int fs = (multiplicity==2) ? particle_kinds[0]*particle_kinds[1] : 0;

  ChooseGenerators(is, fs);

  // Save kinematics for use with distributions
  SaveKinematics(bullet, target);

  // Save particle types for use with distributions
  kinds = particle_kinds;
}

// Save kinematics for use with generators

void G4CascadeFinalStateAlgorithm::
SaveKinematics(G4InuclElementaryParticle* bullet,
	       G4InuclElementaryParticle* target) {
  if (GetVerboseLevel()>1) 
    G4cout << " >>> " << GetName() << "::SaveKinematics" << G4endl;

  if (target->nucleon()) {	// Which particle originated in nucleus?
    toSCM.setBullet(bullet);
    toSCM.setTarget(target);
  } else {
    toSCM.setBullet(target);
    toSCM.setTarget(bullet);
  }

  toSCM.toTheCenterOfMass();

  bullet_ekin = toSCM.getKinEnergyInTheTRS();
}


// Select generator based on initial and final state

void G4CascadeFinalStateAlgorithm::ChooseGenerators(G4int is, G4int fs) {
  if (GetVerboseLevel()>1) 
    G4cout << " >>> " << GetName() << "::ChooseGenerators"
	   << " is " << is << " fs " << fs << G4endl;

  // Get generators for momentum and angle
  if (G4CascadeParameters::usePhaseSpace()) momDist = 0;
  else momDist = G4MultiBodyMomentumDist::GetDist(is, multiplicity);

  if (fs > 0 && multiplicity == 2) {
    G4int kw = (fs==is) ? 1 : 2;
    angDist = G4TwoBodyAngularDist::GetDist(is, fs, kw);
  } else if (multiplicity == 3) {
    angDist = G4TwoBodyAngularDist::GetDist(is);
  } else {
    angDist = 0;
  }

  if (GetVerboseLevel()>1) {
    G4cout << " " << (momDist?momDist->GetName().c_str():"") << " "
	   << (angDist?angDist->GetName().c_str():"") << G4endl;
  }
}


// Two-body generation uses angular-distribution function

void G4CascadeFinalStateAlgorithm::
GenerateTwoBody(G4double initialMass, const std::vector<G4double>& masses,
		std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()>1) 
    G4cout << " >>> " << GetName() << "::GenerateTwoBody" << G4endl;

  finalState.clear();		// Initialization and sanity checks

  if (multiplicity != 2) return;

  // Generate momentum vector in CMS for back-to-back particles
  G4double pscm = TwoBodyMomentum(initialMass, masses[0], masses[1]);

  G4double costh = angDist ? angDist->GetCosTheta(bullet_ekin, pscm)
                           : (2.*G4UniformRand() - 1.);

  mom.setRThetaPhi(pscm, std::acos(costh), UniformPhi());

  if (GetVerboseLevel()>3) {		// Copied from old G4EPCollider
    G4cout << " Particle kinds = " << kinds[0] << " , " << kinds[1]
	   << "\n pmod " << pscm
	   << "\n before rotation px " << mom.x() << " py " << mom.y()
	   << " pz " << mom.z() << G4endl;
  }

  finalState.resize(2);				// Allows filling by index

  finalState[0].setVectM(mom, masses[0]);
  finalState[0] = toSCM.rotate(finalState[0]);

  if (GetVerboseLevel()>3) {		// Copied from old G4EPCollider
    G4cout << " after rotation px " << finalState[0].x() << " py "
	   << finalState[0].y() << " pz " << finalState[0].z() << G4endl;
  }

  finalState[1].setVectM(-finalState[0].vect(), masses[1]);
}


// N-body generation uses momentum-modulus distribution, computed angles

void G4CascadeFinalStateAlgorithm::
GenerateMultiBody(G4double initialMass, const std::vector<G4double>& masses,
		  std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()>1) 
    G4cout << " >>> " << GetName() << "::GenerateMultiBody" << G4endl;

  if (G4CascadeParameters::usePhaseSpace()) {
    FillUsingKopylov(initialMass, masses, finalState);
    return;
  }

  finalState.clear();		// Initialization and sanity checks
  if (multiplicity < 3) return;
  if (!momDist) return;

  G4int itry = -1;		/* Loop checking 08.06.2015 MHK */
  while ((G4int)finalState.size() != multiplicity && ++itry < itry_max) {
    FillMagnitudes(initialMass, masses);
    FillDirections(initialMass, masses, finalState);
  }
}


void G4CascadeFinalStateAlgorithm::
FillMagnitudes(G4double initialMass, const std::vector<G4double>& masses) {
  if (GetVerboseLevel()>1) 
    G4cout << " >>> " << GetName() << "::FillMagnitudes" << G4endl;

  modules.clear();		// Initialization and sanity checks
  if (!momDist) return;

  modules.resize(multiplicity,0.);	// Pre-allocate to avoid resizing

  G4double mass_last = masses.back();
  G4double pmod = 0.;

  if (GetVerboseLevel() > 3){
    G4cout << " knd_last " << kinds.back() << " mass_last " 
           << mass_last << G4endl;
  }

  G4int itry = -1;
  while (++itry < itry_max) {		/* Loop checking 08.06.2015 MHK */
    if (GetVerboseLevel() > 3) {
      G4cout << " itry in fillMagnitudes " << itry << G4endl;
    }

    G4double eleft = initialMass;

    G4int i;	// For access outside of loop
    for (i=0; i < multiplicity-1; i++) {
      pmod = momDist->GetMomentum(kinds[i], bullet_ekin);

      if (pmod < small) break;
      eleft -= std::sqrt(pmod*pmod + masses[i]*masses[i]);

      if (GetVerboseLevel() > 3) {
	G4cout << " kp " << kinds[i] << " pmod " << pmod
	       << " mass2 " << masses[i]*masses[i] << " eleft " << eleft
	       << "\n x1 " << eleft - mass_last << G4endl;
      }

      if (eleft <= mass_last) break;

      modules[i] = pmod;
    }

    if (i < multiplicity-1) continue;	// Failed to generate full kinematics

    G4double plast = eleft * eleft - masses.back()*masses.back();
    if (GetVerboseLevel() > 2) G4cout << " plast ** 2 " << plast << G4endl;
    
    if (plast <= small) continue;	// Not enough momentum left over

    plast = std::sqrt(plast);		// Final momentum is what's left over
    modules.back() = plast;
    
    if (multiplicity > 3 || satisfyTriangle(modules)) break;	// Successful
  }	// while (itry < itry_max)

  if (itry >= itry_max) {		// Too many attempts
    if (GetVerboseLevel() > 2)
      G4cerr << " Unable to generate momenta for multiplicity "
	     << multiplicity << G4endl;

    modules.clear();		// Something went wrong, throw away partial
  }
}

// For three-body states, check kinematics of momentum moduli

G4bool G4CascadeFinalStateAlgorithm::
satisfyTriangle(const std::vector<G4double>& pmod) const {
  if (GetVerboseLevel() > 3) 
    G4cout << " >>> " << GetName() << "::satisfyTriangle" << G4endl;

  return ( (pmod.size() != 3) ||
	   !(pmod[0] < std::fabs(pmod[1] - pmod[2]) ||
	     pmod[0] > pmod[1] + pmod[2] ||
	     pmod[1] < std::fabs(pmod[0] - pmod[2]) ||
	     pmod[1] > pmod[0] + pmod[2] ||
	     pmod[2] < std::fabs(pmod[0] - pmod[1]) ||
	     pmod[2] > pmod[1] + pmod[0])
	   );
}

// Generate momentum directions into final state

void G4CascadeFinalStateAlgorithm::
FillDirections(G4double initialMass, const std::vector<G4double>& masses,
	       std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()>1) 
    G4cout << " >>> " << GetName() << "::FillDirections" << G4endl;

  finalState.clear();			// Initialization and sanity check
  if ((G4int)modules.size() != multiplicity) return;

  // Different order of processing for three vs. N body
  if (multiplicity == 3)
    FillDirThreeBody(initialMass, masses, finalState);
  else
    FillDirManyBody(initialMass, masses, finalState);
}

void G4CascadeFinalStateAlgorithm::
FillDirThreeBody(G4double initialMass, const std::vector<G4double>& masses,
		 std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()>1) 
    G4cout << " >>> " << GetName() << "::FillDirThreeBody" << G4endl;

  finalState.resize(3);

  G4double costh = GenerateCosTheta(kinds[2], modules[2]);
  finalState[2] = generateWithFixedTheta(costh, modules[2], masses[2]);
  finalState[2] = toSCM.rotate(finalState[2]);	// Align target axis

  // Generate direction of first particle
  costh = -0.5 * (modules[2]*modules[2] + modules[0]*modules[0] -
		  modules[1]*modules[1]) / modules[2] / modules[0];

  if (std::fabs(costh) >= maxCosTheta) {  // Bad kinematics; abort generation
    finalState.clear();
    return;
  }

  // Report success
  if (GetVerboseLevel()>2) G4cout << " ok for mult 3" << G4endl;

  // First particle is at fixed angle to recoil system
  finalState[0] = generateWithFixedTheta(costh, modules[0], masses[0]);
  finalState[0] = toSCM.rotate(finalState[2], finalState[0]);

  // Remaining particle is constrained to recoil from entire rest of system
  finalState[1].set(0.,0.,0.,initialMass);
  finalState[1] -= finalState[0] + finalState[2];
}

void G4CascadeFinalStateAlgorithm::
FillDirManyBody(G4double initialMass, const std::vector<G4double>& masses,
		std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()>1) 
    G4cout << " >>> " << GetName() << "::FillDirManyBody" << G4endl;

  // Fill all but the last two particles randomly
  G4double costh = 0.;

  finalState.resize(multiplicity);

  for (G4int i=0; i<multiplicity-2; i++) {
    costh = GenerateCosTheta(kinds[i], modules[i]);
    finalState[i] = generateWithFixedTheta(costh, modules[i], masses[i]);
    finalState[i] = toSCM.rotate(finalState[i]);	// Align target axis
  }

  // Total momentum so far, to compute recoil of last two particles
  G4LorentzVector psum =
    std::accumulate(finalState.begin(), finalState.end()-2, G4LorentzVector());
  G4double pmod = psum.rho();

  costh = -0.5 * (pmod*pmod +
		  modules[multiplicity-2]*modules[multiplicity-2] -
		  modules[multiplicity-1]*modules[multiplicity-1])
    / pmod / modules[multiplicity-2];

  if (GetVerboseLevel() > 2) G4cout << " ct last " << costh << G4endl;

  if (std::fabs(costh) >= maxCosTheta) {  // Bad kinematics; abort generation
    finalState.clear();
    return;
  }

  // Report success
  if (GetVerboseLevel()>2) G4cout << " ok for mult " << multiplicity << G4endl;

  // First particle is at fixed angle to recoil system
  finalState[multiplicity-2] =
    generateWithFixedTheta(costh, modules[multiplicity-2],
			   masses[multiplicity-2]);
  finalState[multiplicity-2] = toSCM.rotate(psum, finalState[multiplicity-2]);

  // Remaining particle is constrained to recoil from entire rest of system
  finalState[multiplicity-1].set(0.,0.,0.,initialMass);
  finalState[multiplicity-1] -= psum + finalState[multiplicity-2];
}


// Generate polar angle for three- and multi-body systems

G4double G4CascadeFinalStateAlgorithm::
GenerateCosTheta(G4int ptype, G4double pmod) const {
  if (GetVerboseLevel() > 2) {
    G4cout << " >>> " << GetName() << "::GenerateCosTheta " << ptype
	   << " " << pmod << G4endl;
  }

  if (multiplicity == 3) {		// Use distribution for three-body
    return angDist->GetCosTheta(bullet_ekin, ptype);
  }

  // Throw multi-body distribution
  G4double p0 = ptype<3 ? 0.36 : 0.25;	// Nucleon vs. everything else
  G4double alf = 1.0 / p0 / (p0 - (pmod+p0)*G4Exp(-pmod / p0));

  G4double sinth = 2.0;

  G4int itry1 = -1;		/* Loop checking 08.06.2015 MHK */
  while (std::fabs(sinth) > maxCosTheta && ++itry1 < itry_max) {
    G4double s1 = pmod * inuclRndm();
    G4double s2 = alf * oneOverE * p0 * inuclRndm();
    G4double salf = s1 * alf * G4Exp(-s1 / p0);
    if (GetVerboseLevel() > 3) {
      G4cout << " s1 * alf * G4Exp(-s1 / p0) " << salf
	     << " s2 " << s2 << G4endl;
    }
    
    if (salf > s2) sinth = s1 / pmod;
  }
  
  if (GetVerboseLevel() > 3)
    G4cout << " itry1 " << itry1 << " sinth " << sinth << G4endl;
  
  if (itry1 == itry_max) {
    if (GetVerboseLevel() > 2)
      G4cout << " high energy angles generation: itry1 " << itry1 << G4endl;
    
    sinth = 0.5 * inuclRndm();
  }

  // Convert generated sin(theta) to cos(theta) with random sign
  G4double costh = std::sqrt(1.0 - sinth * sinth);
  if (inuclRndm() > 0.5) costh = -costh;

  return costh;
}


// SPECIAL:  Generate N-body phase space using Kopylov algorithm
//	     ==> Code is copied verbatim from G4HadPhaseSpaceKopylov

void G4CascadeFinalStateAlgorithm::
FillUsingKopylov(G4double initialMass,
		 const std::vector<G4double>& masses,
		 std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()>2)
    G4cout << " >>> " << GetName() << "::FillUsingKopylov" << G4endl;

  finalState.clear();

  std::size_t N = masses.size();
  finalState.resize(N);

  G4double mtot = std::accumulate(masses.begin(), masses.end(), 0.0);
  G4double mu = mtot;
  G4double Mass = initialMass;
  G4double T = Mass-mtot;
  G4double recoilMass = 0.0;
  G4ThreeVector momV, boostV;		// Buffers to reduce memory churn
  G4LorentzVector recoil(0.0,0.0,0.0,Mass);

  for (std::size_t k=N-1; k>0; --k) {
    mu -= masses[k];
    T *= (k>1) ? BetaKopylov((G4int)k) : 0.;
    
    recoilMass = mu + T;

    boostV = recoil.boostVector();	// Previous system's rest frame
    
    // Create momentum with a random direction isotropically distributed
    // FIXME:  Should theta distribution should use Bertini fit function?
    momV.setRThetaPhi(TwoBodyMomentum(Mass,masses[k],recoilMass),
		      UniformTheta(), UniformPhi());
    
    finalState[k].setVectM(momV,masses[k]);
    recoil.setVectM(-momV,recoilMass);

    finalState[k].boost(boostV);
    recoil.boost(boostV);
    Mass = recoilMass;
  }
  
  finalState[0] = recoil;
}

G4double G4CascadeFinalStateAlgorithm::BetaKopylov(G4int K) const {
  G4Pow* g4pow = G4Pow::GetInstance();

  G4int N = 3*K - 5;
  G4double xN = G4double(N);
  G4double Fmax = std::sqrt(g4pow->powN(xN/(xN+1.),N)/(xN+1.)); 

  G4double F, chi;
  do {					/* Loop checking 08.06.2015 MHK */
    chi = G4UniformRand();
    F = std::sqrt(g4pow->powN(chi,N)*(1.-chi));      
  } while ( Fmax*G4UniformRand() > F);
  return chi;
}
