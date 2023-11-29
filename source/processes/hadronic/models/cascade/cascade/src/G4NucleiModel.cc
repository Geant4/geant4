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
// For the best approximation to a physical-units model, set the following:
//	setenv G4NUCMODEL_XSEC_SCALE   0.1
//	setenv G4NUCMODEL_RAD_SCALE    1.0
//	setenv G4NUCMODEL_RAD_2PAR     1
//	setenv G4NUCMODEL_RAD_SMALL    1.992
//	setenv G4NUCMODEL_RAD_ALPHA    0.84
//	setenv G4NUCMODEL_FERMI_SCALE  0.685
//	setenv G4NUCMODEL_RAD_TRAILING 0.70
//
// 20100112  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100114  M. Kelsey -- Use G4ThreeVector for position
// 20100309  M. Kelsey -- Use new generateWithRandomAngles for theta,phi stuff;
//		eliminate some unnecessary std::pow(), move ctor here
// 20100407  M. Kelsey -- Replace std::vector<>::resize(0) with ::clear().
//		Move output vectors from generateParticleFate() and 
//		::generateInteractionPartners() to be data members; return via
//		const-ref instead of by value.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100418  M. Kelsey -- Reference output particle lists via const-ref
// 20100421  M. Kelsey -- Replace hardwired p/n masses with G4PartDef's
// 20100517  M. Kelsey -- Use G4CascadeINterpolator for cross-section
//		calculations.  use G4Cascade*Channel for total xsec in pi-N
//		N-N channels.  Move absorptionCrossSection() from SpecialFunc.
// 20100610  M. Kelsey -- Replace another random-angle code block; add some
//		diagnostic output for partner-list production.
// 20100617  M. Kelsey -- Replace preprocessor flag CHC_CHECK with
//		G4CASCADE_DEBUG_CHARGE
// 20100620  M. Kelsey -- Improve error message on empty partners list, add
//		four-momentum checking after EPCollider
// 20100621  M. Kelsey -- In boundaryTransition() account for momentum transfer
//		to secondary by boosting into recoil nucleus "rest" frame.
//		Don't need temporaries to create dummy "partners" for list.
// 20100622  M. Kelsey -- Restore use of "bindingEnergy()" function name, which
//		is now a wrapper for G4NucleiProperties::GetBindingEnergy().
// 20100623  M. Kelsey -- Eliminate some temporaries terminating partner-list,
//		discard recoil boost for now. Initialize all data
//		members in ctors.  Allow generateModel() to be called
//		mutliple times, by clearing vectors each time through;
//		avoid extra work by returning if A and Z are same as
//		before.
// 20100628  M. Kelsey -- Two momentum-recoil bugs; don't subtract energies!
// 20100715  M. Kelsey -- Make G4InuclNuclei in generateModel(), use for
//		balance checking (only verbose>2) in generateParticleFate().
// 20100721  M. Kelsey -- Use new G4CASCADE_CHECK_ECONS for conservation checks
// 20100723  M. Kelsey -- Move G4CollisionOutput buffer to .hh for reuse
// 20100726  M. Kelsey -- Preallocate arrays with number_of_zones dimension.
// 20100902  M. Kelsey -- Remove resize(3) directives from qdeutron/acsecs
// 20100907  M. Kelsey -- Limit interaction targets based on current nucleon
//		counts (e.g., no pp if protonNumberCurrent < 2!)
// 20100914  M. Kelsey -- Migrate to integer A and Z
// 20100928  M. Kelsey -- worthToPropagate() use nuclear potential for hadrons.
// 20101005  M. Kelsey -- Annotate hardwired numbers for upcoming rationalizing,
//		move hardwired numbers out to static data members, factor out
//		all the model construction pieces to separate functions, fix
//		wrong value for 4/3 pi.
// 20101019  M. Kelsey -- CoVerity reports: unitialized constructor, dtor leak
// 20101020  M. Kelsey -- Bug fixes to refactoring changes (5 Oct).  Back out
//		worthToPropagate() changes for better regression testing.
// 20101020  M. Kelsey -- Re-activate worthToPropagate() changes.
// 20101119  M. Kelsey -- Hide "negative path" and "no partners" messages in
//		verbosity.
// 20110218  M. Kelsey -- Add crossSectionUnits and radiusUnits scale factors,
//		use "theoretical" numbers for radii etc., multipled by scale
//		factor; set scale factors using environment variables
// 20110303  M. Kelsey -- Add comments why using fabs() with B.E. differences?
// 20110321  M. Kelsey -- Replace strtof() with strtod() for envvar conversion
// 20110321  M. Kelsey -- Use fm and fm^2 as default units, Per D. Wright
//		(NOTE: Restored from original 20110318 commit)
// 20110324  D. Wright -- Implement trailing effect
// 20110324  M. Kelsey -- Move ::reset() here, as it has more code.
// 20110519  M. Kelsey -- Used "rho" after assignment, instead of recomputing
// 20110525  M. Kelsey -- Revert scale factor changes (undo 20110321 changes)
// 20110617  M. Kelsey -- Apply scale factor to trailing-effect radius, make
//		latter runtime adjustable (G4NUCMODEL_RAD_TRAILING)
// 20110720  M. Kelsey -- Follow interface change for cross-section tables,
//		eliminating switch blocks.
// 20110806  M. Kelsey -- Reduce memory churn by pre-allocating buffers
// 20110823  M. Kelsey -- Remove local cross-section tables entirely
// 20110825  M. Kelsey -- Add comments regarding Fermi momentum scale, set of
//		"best guess" parameter values
// 20110831  M. Kelsey -- Make "best guess" parameters the defaults
// 20110922  M. Kelsey -- Follow migrations G4InuclParticle::print(ostream&)
//		and G4CascadParticle::print(ostream&)
// 20111018  M. Kelsey -- Correct kaon potential to be positive, not negative
// 20111107  M. Kelsey -- *** REVERT TO OLD NON-PHYSICAL PARAMETERS FOR 9.5 ***
// 20120306  M. Kelsey -- For incident projectile (CParticle incoming to
//		nucleus from outside) don't use pw cut; force interaction.
// 20120307  M. Kelsey -- Compute zone volumes during setup, not during each
//		interaction.  Include photons in absorption by dibaryons.
//		Discard unused code in totalCrossSection().  Encapsulate
//		interaction path calculations in generateInteractionLength().
// 20120308  M. Kelsey -- Add binned photon absorption cross-section.
// 20120320  M. Kelsey -- Bug fix: fill zone_volumes for single-zone case
// 20120321  M. Kelsey -- Add check in isProjectile() for single-zone case.
// 20120608  M. Kelsey -- Fix variable-name "shadowing" compiler warnings.
// 20120822  M. Kelsey -- Move envvars to G4CascadeParameters.
// 20121205  M. Kelsey -- Bug fix in generateParticleFate(): daughters should
//		have generation count incremented from parent.  This allows
//		isProjectile() to simply check generation number == 0.
// 20130221  M. Kelsey -- For incident photons, move to random point along
//		trajectory through nucleus for first (forced) interaction.
// 20130223  M. Kelsey -- Fix bugs in weighting and selection of traj points
// 20130302  M. Kelsey -- Replace "isProjectile()" with "forceFirst()" in
//		path-length calculation.
// 20130314  M. Kelsey -- Restore null initializer and if-block for _TLS_.
// 20130508  D. Wright -- Support muon-nuclear interactions.
// 20130510  M. Kelsey -- Check for zero-momentum in choosePointAlongTraj().
// 20130511  M. Kelsey -- Move choosePointAlongTraj() to initializeCascad(),
//		instead of after generateIP().  Change "spath<path" to <= to
//		handle at-rest particles.  Skip mu-/neutron interactions.
//		Hide "zone 0 transition" error message behind verbosity.
// 20130523  M. Kelsey -- Bug fix: Initialize dummy_convertor in inverseMFP.
//		For capture-at-rest, replace exterior nucleus with outer zone.
// 20130524  M. Kelsey -- Move "large" and "small" cutoffs to static const.
//		Check zone argument to inverseMFP() to ensure within range.
// 20130611  M. Kelsey -- Undo "spath<=path" change (20130511), replace with
//		explicit check on spath==0.
// 20130619  A. Ribon  -- Fixed reproducibility problem in the method
//              generateParticleFate
// 20130627  M. Kelsey -- Check "path==0.", rather than spath.
// 20130628  M. Kelsey -- Print deuteron type code, rather than array index,
//		Extend useQuasiDeuteron() to check good absorption
// 20130701  M. Kelsey -- Don't average 1/MFP for total interaction; just sum;
//		inverseMFP() returns "large" value instead of input path.
// 20130806  M. Kelsey -- Move static G4InuclEP's to file scope to avoid
//		thread collisions
// 20130924  M. Kelsey -- Use G4Log, G4Exp, G4Pow for CPU speedup
// 20130925  M. Kelsey -- Eliminate unnecessary use of pow() in integrals
// 20131001  M. Kelsey -- Move QDinterp object to data member, thread local
// 20140116  M. Kelsey -- Move all envvar parameters to const data members
// 20141001  M. Kelsey -- Change sign of "dv" in boundaryTransition
// 20150608  M. Kelsey -- Label all while loops as terminating.
// 20150622  M. Kelsey -- Use G4AutoDelete for _TLS_ buffers.
// 20180227  A. Ribon  -- Replaced obsolete std::bind2nd with std::bind

#include "G4NucleiModel.hh"
#include "G4AutoDelete.hh"
#include "G4CascadeChannel.hh"
#include "G4CascadeChannelTables.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CascadeInterpolator.hh"
#include "G4CascadeParameters.hh"
#include "G4CollisionOutput.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4Exp.hh"
#include "Randomize.hh"
#include "G4HadTmpUtil.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticleNames.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4Log.hh"
#include "G4LorentzConvertor.hh"
#include "G4Neutron.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleLargerBeta.hh"
#include "G4PhysicalConstants.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>

using namespace G4InuclParticleNames;
using namespace G4InuclSpecialFunctions;

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

// Cut-off limits for extreme values in calculations
const G4double G4NucleiModel::small = 1.0e-9;
const G4double G4NucleiModel::large = 1000.;

// For convenience in computing densities
const G4double G4NucleiModel::piTimes4thirds = pi*4./3.;

// Zone boundaries as fraction of nuclear radius (from outside in)
const G4double G4NucleiModel::alfa3[3] = { 0.7, 0.3, 0.01 };
const G4double G4NucleiModel::alfa6[6] = { 0.9, 0.6, 0.4, 0.2, 0.1, 0.05 };

// Flat nuclear potentials for mesons and hyperons (GeV)
const G4double G4NucleiModel::pion_vp       = 0.007;
const G4double G4NucleiModel::pion_vp_small = 0.007;
const G4double G4NucleiModel::kaon_vp       = 0.015;
const G4double G4NucleiModel::hyperon_vp    = 0.030;

namespace {
  // Quasideuteron absorption cross section is taken to be the 
  // deuteron photo-disintegration cross section (gam + D -> p + n)
  // Values taken from smooth curve drawn through data from 2.4 - 500 MeV, 
  // plus angle integration of JLab data from 0.5 - 3 GeV
  //    M. Mirazita et al., Phys. Rev. C 70, 014005 (2004)
  // Points above 3 GeV are extrapolated.
  const G4double kebins[] = {
      0.0,  0.0024, 0.0032, 0.0042, 0.0056, 0.0075, 0.01,  0.024, 0.075, 0.1,
      0.13, 0.18,   0.24,   0.32,   0.42,   0.56,   0.75,  1.0,   1.3,   1.8,
      2.4,  3.2,    4.2,    5.6,    7.5,   10.0,   13.0,  18.0,  24.0,  32.0 };

  const G4double gammaQDxsec[30] = {
     0.0,     0.7,    2.0,    2.2,    2.1,    1.8,    1.3,     0.4,     0.098,   0.071,
     0.055,   0.055,  0.065,  0.045,  0.017,  0.007,  2.37e-3, 6.14e-4, 1.72e-4, 4.2e-5,
     1.05e-5, 3.0e-6, 7.0e-7, 1.3e-7, 2.3e-8, 3.2e-9, 4.9e-10, 0.0,     0.0,     0.0 };
}


// Constructors and destructor
G4NucleiModel::G4NucleiModel()
  : verboseLevel(0), nuclei_radius(0.), nuclei_volume(0.), number_of_zones(0),
    A(0), Z(0), theNucleus(0), neutronNumber(0), protonNumber(0),
    neutronNumberCurrent(0), protonNumberCurrent(0), current_nucl1(0),
    current_nucl2(0), gammaQDinterp(kebins),
    crossSectionUnits(G4CascadeParameters::xsecScale()),
    radiusUnits(G4CascadeParameters::radiusScale()),
    skinDepth(0.611207*radiusUnits),
    radiusScale((G4CascadeParameters::useTwoParam()?1.16:1.2)*radiusUnits),
    radiusScale2((G4CascadeParameters::useTwoParam()?-1.3456:0.)*radiusUnits),
    radiusForSmall(G4CascadeParameters::radiusSmall()),
    radScaleAlpha(G4CascadeParameters::radiusAlpha()),
    fermiMomentum(G4CascadeParameters::fermiScale()),
    R_nucleon(G4CascadeParameters::radiusTrailing()),
    gammaQDscale(G4CascadeParameters::gammaQDScale()),
    potentialThickness(1.0),
    neutronEP(neutron), protonEP(proton) {}

G4NucleiModel::G4NucleiModel(G4int a, G4int z)
  : verboseLevel(0), nuclei_radius(0.), nuclei_volume(0.), number_of_zones(0),
    A(0), Z(0), theNucleus(0), neutronNumber(0), protonNumber(0),
    neutronNumberCurrent(0), protonNumberCurrent(0), current_nucl1(0),
    current_nucl2(0), gammaQDinterp(kebins),
    crossSectionUnits(G4CascadeParameters::xsecScale()),
    radiusUnits(G4CascadeParameters::radiusScale()),
    skinDepth(0.611207*radiusUnits),
    radiusScale((G4CascadeParameters::useTwoParam()?1.16:1.2)*radiusUnits),
    radiusScale2((G4CascadeParameters::useTwoParam()?-1.3456:0.)*radiusUnits),
    radiusForSmall(G4CascadeParameters::radiusSmall()),
    radScaleAlpha(G4CascadeParameters::radiusAlpha()),
    fermiMomentum(G4CascadeParameters::fermiScale()),
    R_nucleon(G4CascadeParameters::radiusTrailing()),
    gammaQDscale(G4CascadeParameters::gammaQDScale()),
    potentialThickness(1.0),
    neutronEP(neutron), protonEP(proton) {
  generateModel(a,z);
}

G4NucleiModel::G4NucleiModel(G4InuclNuclei* nuclei)
  : verboseLevel(0), nuclei_radius(0.), nuclei_volume(0.), number_of_zones(0),
    A(0), Z(0), theNucleus(0), neutronNumber(0), protonNumber(0),
    neutronNumberCurrent(0), protonNumberCurrent(0), current_nucl1(0),
    current_nucl2(0), gammaQDinterp(kebins),
    crossSectionUnits(G4CascadeParameters::xsecScale()),
    radiusUnits(G4CascadeParameters::radiusScale()),
    skinDepth(0.611207*radiusUnits),
    radiusScale((G4CascadeParameters::useTwoParam()?1.16:1.2)*radiusUnits),
    radiusScale2((G4CascadeParameters::useTwoParam()?-1.3456:0.)*radiusUnits),
    radiusForSmall(G4CascadeParameters::radiusSmall()),
    radScaleAlpha(G4CascadeParameters::radiusAlpha()),
    fermiMomentum(G4CascadeParameters::fermiScale()),
    R_nucleon(G4CascadeParameters::radiusTrailing()),
    gammaQDscale(G4CascadeParameters::gammaQDScale()),
    potentialThickness(1.0),
    neutronEP(neutron), protonEP(proton) {
  generateModel(nuclei);
}

G4NucleiModel::~G4NucleiModel() {
  delete theNucleus;
  theNucleus = 0;
}


// Initialize model state for new cascade

void G4NucleiModel::reset(G4int nHitNeutrons, G4int nHitProtons,
			  const std::vector<G4ThreeVector>* hitPoints) {
  neutronNumberCurrent = neutronNumber - nHitNeutrons;
  protonNumberCurrent  = protonNumber - nHitProtons;
  
  // zero or copy collision point array for trailing effect
  if (!hitPoints || !hitPoints->empty()) collisionPts.clear();
  else collisionPts = *hitPoints;
}


// Generate nuclear model parameters for given nucleus

void G4NucleiModel::generateModel(G4InuclNuclei* nuclei) {
  generateModel(nuclei->getA(), nuclei->getZ());
}

void G4NucleiModel::generateModel(G4int a, G4int z) {
  if (verboseLevel) {
    G4cout << " >>> G4NucleiModel::generateModel A " << a << " Z " << z
	   << G4endl;
  }

  // If model already built, just return; otherwise initialize everything
  if (a == A && z == Z) {
    if (verboseLevel > 1) G4cout << " model already generated" << z << G4endl;
    reset();
    return;
  }

  A = a;
  Z = z;
  delete theNucleus;
  theNucleus = new G4InuclNuclei(A,Z);		// For conservation checking

  neutronNumber = A - Z;
  protonNumber = Z;
  reset();

  if (verboseLevel > 3) {
    G4cout << "  crossSectionUnits = " << crossSectionUnits << G4endl
	   << "  radiusUnits = " << radiusUnits << G4endl
	   << "  skinDepth = " << skinDepth << G4endl
	   << "  radiusScale = " << radiusScale << G4endl
	   << "  radiusScale2 = " << radiusScale2 << G4endl
	   << "  radiusForSmall = " << radiusForSmall << G4endl
	   << "  radScaleAlpha  = " << radScaleAlpha << G4endl
	   << "  fermiMomentum = " << fermiMomentum << G4endl
	   << "  piTimes4thirds = " << piTimes4thirds << G4endl;
  }

  G4double nuclearRadius;		// Nuclear radius computed from A
  if (A>4) nuclearRadius = radiusScale*G4cbrt(A) + radiusScale2/G4cbrt(A);
  else     nuclearRadius = radiusForSmall * (A==4 ? radScaleAlpha : 1.);

  // This will be used to pre-allocate lots of arrays below
  number_of_zones = (A < 5) ? 1 : (A < 100) ? 3 : 6;

  // Clear all parameters arrays for reloading
  binding_energies.clear();
  nucleon_densities.clear();
  zone_potentials.clear();
  fermi_momenta.clear();
  zone_radii.clear();
  zone_volumes.clear();

  fillBindingEnergies();
  fillZoneRadii(nuclearRadius);

  G4double tot_vol = fillZoneVolumes(nuclearRadius);	// Woods-Saxon integral

  fillPotentials(proton, tot_vol);		// Protons
  fillPotentials(neutron, tot_vol);		// Neutrons

  // Additional flat zone potentials for other hadrons
  const std::vector<G4double> vp(number_of_zones, (A>4)?pion_vp:pion_vp_small);
  const std::vector<G4double> kp(number_of_zones, kaon_vp);
  const std::vector<G4double> hp(number_of_zones, hyperon_vp);

  zone_potentials.push_back(vp);
  zone_potentials.push_back(kp);
  zone_potentials.push_back(hp);

  nuclei_radius = zone_radii.back();
  nuclei_volume = std::accumulate(zone_volumes.begin(),zone_volumes.end(),0.);

  if (verboseLevel > 3) printModel();
}


// Load binding energy array for current nucleus

void G4NucleiModel::fillBindingEnergies() {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::fillBindingEnergies" << G4endl;

  G4double dm = bindingEnergy(A,Z);

  // Binding energy differences for proton and neutron loss, respectively
  // FIXME:  Why is fabs() used here instead of the signed difference?
  binding_energies.push_back(std::fabs(bindingEnergy(A-1,Z-1)-dm)/GeV);
  binding_energies.push_back(std::fabs(bindingEnergy(A-1,Z)-dm)/GeV);
}

// Load zone boundary radius array for current nucleus

void G4NucleiModel::fillZoneRadii(G4double nuclearRadius) {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::fillZoneRadii" << G4endl;

  G4double skinRatio = nuclearRadius/skinDepth;
  G4double skinDecay = G4Exp(-skinRatio);    

  if (A < 5) {			// Light ions treated as simple balls
    zone_radii.push_back(nuclearRadius);
    ur[0] = 0.;
    ur[1] = 1.;
  } else if (A < 12) {		// Small nuclei have Gaussian potential
    G4double rSq = nuclearRadius * nuclearRadius;
    G4double gaussRadius = std::sqrt(rSq * (1.0 - 1.0/A) + 6.4);
    
    ur[0] = 0.0;
    for (G4int i = 0; i < number_of_zones; i++) {
      G4double y = std::sqrt(-G4Log(alfa3[i]));
      zone_radii.push_back(gaussRadius * y);
      ur[i+1] = y;
    }
  } else if (A < 100) {		// Intermediate nuclei have three-zone W-S
    ur[0] = -skinRatio;
    for (G4int i = 0; i < number_of_zones; i++) {
      G4double y = G4Log((1.0 + skinDecay)/alfa3[i] - 1.0);
      zone_radii.push_back(nuclearRadius + skinDepth * y);
      ur[i+1] = y;
    }
  } else {			// Heavy nuclei have six-zone W-S
    ur[0] = -skinRatio;
    for (G4int i = 0; i < number_of_zones; i++) {
      G4double y = G4Log((1.0 + skinDecay)/alfa6[i] - 1.0);
      zone_radii.push_back(nuclearRadius + skinDepth * y);
      ur[i+1] = y;
    }
  }
}

// Compute zone-by-zone density integrals

G4double G4NucleiModel::fillZoneVolumes(G4double nuclearRadius) {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::fillZoneVolumes" << G4endl;

  G4double tot_vol = 0.;	// Return value omitting 4pi/3 for sphere!

  if (A < 5) {			// Light ions are treated as simple balls
    v[0] = v1[0] = 1.;
    tot_vol = zone_radii[0]*zone_radii[0]*zone_radii[0];
    zone_volumes.push_back(tot_vol*piTimes4thirds);	// True volume of sphere

    return tot_vol;
  }

  PotentialType usePotential = (A < 12) ? Gaussian : WoodsSaxon;

  for (G4int i = 0; i < number_of_zones; i++) {
    if (usePotential == WoodsSaxon) {
      v[i] = zoneIntegralWoodsSaxon(ur[i], ur[i+1], nuclearRadius);
    } else {
      v[i] = zoneIntegralGaussian(ur[i], ur[i+1], nuclearRadius);
    }
    
    tot_vol += v[i];
    
    v1[i] = zone_radii[i]*zone_radii[i]*zone_radii[i];
    if (i > 0) v1[i] -= zone_radii[i-1]*zone_radii[i-1]*zone_radii[i-1];

    zone_volumes.push_back(v1[i]*piTimes4thirds);	// True volume of shell
  }

  return tot_vol;	// Sum of zone integrals, not geometric volume
}

// Load potentials for nucleons (using G4InuclParticle codes)
void G4NucleiModel::fillPotentials(G4int type, G4double tot_vol) {
  if (verboseLevel > 1) 
    G4cout << " >>> G4NucleiModel::fillZoneVolumes(" << type << ")" << G4endl;

  if (type != proton && type != neutron) return;

  const G4double mass = G4InuclElementaryParticle::getParticleMass(type);

  // FIXME:  This is the fabs() binding energy difference, not signed
  const G4double dm = binding_energies[type-1];

  rod.clear(); rod.reserve(number_of_zones);
  pf.clear();  pf.reserve(number_of_zones);
  vz.clear();  vz.reserve(number_of_zones);

  G4int nNucleons = (type==proton) ? protonNumber : neutronNumber;
  G4double dd0 = nNucleons / tot_vol / piTimes4thirds;

  for (G4int i = 0; i < number_of_zones; i++) {
    G4double rd = dd0 * v[i] / v1[i];
    rod.push_back(rd);
    G4double pff = fermiMomentum * G4cbrt(rd);
    pf.push_back(pff);
    vz.push_back(0.5 * pff * pff / mass + dm);
  }
  
  nucleon_densities.push_back(rod);
  fermi_momenta.push_back(pf);
  zone_potentials.push_back(vz);
}

// Zone integral of Woods-Saxon density function
G4double G4NucleiModel::zoneIntegralWoodsSaxon(G4double r1, G4double r2, 
					       G4double nuclearRadius) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::zoneIntegralWoodsSaxon" << G4endl;
  }

  const G4double epsilon = 1.0e-3;
  const G4int itry_max = 1000;

  G4double skinRatio = nuclearRadius / skinDepth;

  G4double d2 = 2.0 * skinRatio;
  G4double dr = r2 - r1;
  G4double fr1 = r1 * (r1 + d2) / (1.0 + G4Exp(r1));
  G4double fr2 = r2 * (r2 + d2) / (1.0 + G4Exp(r2));
  G4double fi = (fr1 + fr2) / 2.;
  G4double fun1 = fi * dr;
  G4double fun;
  G4int jc = 1;
  G4double dr1 = dr;
  G4int itry = 0;

  while (itry < itry_max) {	/* Loop checking 08.06.2015 MHK */
    dr /= 2.;
    itry++;

    G4double r = r1 - dr;
    fi = 0.0;

    for (G4int i = 0; i < jc; i++) { 
      r += dr1; 
      fi += r * (r + d2) / (1.0 + G4Exp(r));
    }

    fun = 0.5 * fun1 + fi * dr;

    if (std::fabs((fun - fun1) / fun) <= epsilon) break;

    jc *= 2;
    dr1 = dr;
    fun1 = fun;
  }	// while (itry < itry_max)

  if (verboseLevel > 2 && itry == itry_max)
    G4cout << " zoneIntegralWoodsSaxon-> n iter " << itry_max << G4endl;

  G4double skinDepth3 = skinDepth*skinDepth*skinDepth;

  return skinDepth3 * (fun + skinRatio*skinRatio*G4Log((1.0 + G4Exp(-r1)) / (1.0 + G4Exp(-r2))));
}


// Zone integral of Gaussian density function
G4double G4NucleiModel::zoneIntegralGaussian(G4double r1, G4double r2, 
					     G4double nucRad) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::zoneIntegralGaussian" << G4endl;
  }

  G4double gaussRadius = std::sqrt(nucRad*nucRad * (1.0 - 1.0/A) + 6.4);

  const G4double epsilon = 1.0e-3;
  const G4int itry_max = 1000;

  G4double dr = r2 - r1;
  G4double fr1 = r1 * r1 * G4Exp(-r1 * r1);
  G4double fr2 = r2 * r2 * G4Exp(-r2 * r2);
  G4double fi = (fr1 + fr2) / 2.;
  G4double fun1 = fi * dr;
  G4double fun;
  G4int jc = 1;
  G4double dr1 = dr;
  G4int itry = 0;

  while (itry < itry_max) {	/* Loop checking 08.06.2015 MHK */
    dr /= 2.;
    itry++;
    G4double r = r1 - dr;
    fi = 0.0;

    for (G4int i = 0; i < jc; i++) { 
      r += dr1; 
      fi += r * r * G4Exp(-r * r);
    }

    fun = 0.5 * fun1 + fi * dr;  

    if (std::fabs((fun - fun1) / fun) <= epsilon) break;

    jc *= 2;
    dr1 = dr;
    fun1 = fun;
  }	// while (itry < itry_max)

  if (verboseLevel > 2 && itry == itry_max)
    G4cerr << " zoneIntegralGaussian-> n iter " << itry_max << G4endl;

  return gaussRadius*gaussRadius*gaussRadius * fun;
}


void G4NucleiModel::printModel() const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::printModel" << G4endl;
  }

  G4cout << " nuclei model for A " << A << " Z " << Z << G4endl
	 << " proton binding energy " << binding_energies[0]
	 << " neutron binding energy " << binding_energies[1] << G4endl
	 << " Nuclei radius " << nuclei_radius << " volume " << nuclei_volume
	 << " number of zones " << number_of_zones << G4endl;

  for (G4int i = 0; i < number_of_zones; i++)
    G4cout << " zone " << i+1 << " radius " << zone_radii[i] 
	   << " volume " << zone_volumes[i] << G4endl
	   << " protons: density " << getDensity(1,i) << " PF " << 
      getFermiMomentum(1,i) << " VP " << getPotential(1,i) << G4endl
	   << " neutrons: density " << getDensity(2,i) << " PF " << 
      getFermiMomentum(2,i) << " VP " << getPotential(2,i) << G4endl
	   << " pions: VP " << getPotential(3,i) << G4endl;
}


G4double G4NucleiModel::getFermiKinetic(G4int ip, G4int izone) const {
  G4double ekin = 0.0;
  
  if (ip < 3 && izone < number_of_zones) {	// ip for proton/neutron only
    G4double pfermi = fermi_momenta[ip - 1][izone]; 
    G4double mass = G4InuclElementaryParticle::getParticleMass(ip);
    
    ekin = std::sqrt(pfermi*pfermi + mass*mass) - mass;
  }  
  return ekin;
}


G4LorentzVector 
G4NucleiModel::generateNucleonMomentum(G4int type, G4int zone) const {
  G4double pmod = getFermiMomentum(type, zone) * G4cbrt(inuclRndm());
  G4double mass = G4InuclElementaryParticle::getParticleMass(type);

  return generateWithRandomAngles(pmod, mass);
}


G4InuclElementaryParticle 
G4NucleiModel::generateNucleon(G4int type, G4int zone) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::generateNucleon" << G4endl;
  }

  G4LorentzVector mom = generateNucleonMomentum(type, zone);
  return G4InuclElementaryParticle(mom, type);
}


G4InuclElementaryParticle
G4NucleiModel::generateQuasiDeuteron(G4int type1, G4int type2,
				    G4int zone) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::generateQuasiDeuteron" << G4endl;
  }

  // Quasideuteron consists of an unbound but associated nucleon pair
  
  // FIXME:  Why generate two separate nucleon momenta (randomly!) and
  //         add them, instead of just throwing a net momentum for the
  //	     dinulceon state?  And why do I have to capture the two
  //	     return values into local variables?
  G4LorentzVector mom1 = generateNucleonMomentum(type1, zone);
  G4LorentzVector mom2 = generateNucleonMomentum(type2, zone);
  G4LorentzVector dmom = mom1+mom2;

  G4int dtype = 0;
       if (type1*type2 == pro*pro) dtype = 111;
  else if (type1*type2 == pro*neu) dtype = 112;
  else if (type1*type2 == neu*neu) dtype = 122;

  return G4InuclElementaryParticle(dmom, dtype);
}


void
G4NucleiModel::generateInteractionPartners(G4CascadParticle& cparticle) {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::generateInteractionPartners" << G4endl;
  }

  thePartners.clear();		// Reset buffer for next cycle

  G4int ptype = cparticle.getParticle().type();
  G4int zone = cparticle.getCurrentZone();

  G4double r_in;		// Destination radius if inbound
  G4double r_out;		// Destination radius if outbound

  if (zone == number_of_zones) {
    r_in = nuclei_radius;
    r_out = 0.0;
  } else if (zone == 0) { 			// particle is outside core
    r_in = 0.0;
    r_out = zone_radii[0];
  } else {
    r_in = zone_radii[zone - 1];
    r_out = zone_radii[zone];
  }  

  G4double path = cparticle.getPathToTheNextZone(r_in, r_out);

  if (verboseLevel > 2) {
    if (isProjectile(cparticle)) G4cout << " incident particle: ";
    G4cout << " r_in " << r_in << " r_out " << r_out << " path " << path
	   << G4endl;
  }

  if (path < -small) { 			// something wrong
    if (verboseLevel)
      G4cerr << " generateInteractionPartners-> negative path length" << G4endl;
    return;
  }

  if (std::fabs(path) < small) {	// Not moving, or just at boundary
    if (cparticle.getMomentum().vect().mag() > small) {
      if (verboseLevel > 3)
	G4cout << " generateInteractionPartners-> zero path" << G4endl;
      
      thePartners.push_back(partner());	// Dummy list terminator with zero path
      return;
    }

    if (zone >= number_of_zones)	// Place captured-at-rest in nucleus
      zone = number_of_zones-1;
  }

  G4double invmfp = 0.;			// Buffers for interaction probability
  G4double spath  = 0.;
  for (G4int ip = 1; ip < 3; ip++) {
    // Only process nucleons which remain active in target
    if (ip==proton  && protonNumberCurrent < 1) continue;
    if (ip==neutron && neutronNumberCurrent < 1) continue;
    if (ip==neutron && ptype==muonMinus) continue;	// mu-/n forbidden

    // All nucleons are assumed to be at rest when colliding
    G4InuclElementaryParticle particle = generateNucleon(ip, zone);
    invmfp = inverseMeanFreePath(cparticle, particle);
    spath  = generateInteractionLength(cparticle, path, invmfp);

    if (path<small || spath < path) {
      if (verboseLevel > 3) {
	G4cout << " adding partner[" << thePartners.size() << "]: "
	       << particle << G4endl;
      }
      thePartners.push_back(partner(particle, spath));
    }
  }	// for (G4int ip...
  
  if (verboseLevel > 2) {
    G4cout << " after nucleons " << thePartners.size() << " path " << path << G4endl;
  }

  // Absorption possible for pions or photons interacting with dibaryons
  if (useQuasiDeuteron(cparticle.getParticle().type())) {
    if (verboseLevel > 2) {
      G4cout << " trying quasi-deuterons with bullet: "
	     << cparticle.getParticle() << G4endl;
    }
  
    // Initialize buffers for quasi-deuteron results
    qdeutrons.clear();
    acsecs.clear();
  
    G4double tot_invmfp = 0.0;		// Total inv. mean-free-path for all QDs

    // Proton-proton state interacts with pi-, mu- or neutrals
    if (protonNumberCurrent >= 2 && ptype != pip) {
      G4InuclElementaryParticle ppd = generateQuasiDeuteron(pro, pro, zone);
      if (verboseLevel > 2)
	G4cout << " ptype=" << ptype << " using pp target\n" << ppd << G4endl;
      
      invmfp = inverseMeanFreePath(cparticle, ppd);      
      tot_invmfp += invmfp;
      acsecs.push_back(invmfp);
      qdeutrons.push_back(ppd);
    }
    
    // Proton-neutron state interacts with any pion type or photon
    if (protonNumberCurrent >= 1 && neutronNumberCurrent >= 1) {
      G4InuclElementaryParticle npd = generateQuasiDeuteron(pro, neu, zone);
      if (verboseLevel > 2) 
	G4cout << " ptype=" << ptype << " using np target\n" << npd << G4endl;
      
      invmfp = inverseMeanFreePath(cparticle, npd);
      tot_invmfp += invmfp;
      acsecs.push_back(invmfp);
      qdeutrons.push_back(npd);
    }
    
    // Neutron-neutron state interacts with pi+ or neutrals
    if (neutronNumberCurrent >= 2 && ptype != pim && ptype != mum) {
      G4InuclElementaryParticle nnd = generateQuasiDeuteron(neu, neu, zone);
      if (verboseLevel > 2)
	G4cout << " ptype=" << ptype << " using nn target\n" << nnd << G4endl;
      
      invmfp = inverseMeanFreePath(cparticle, nnd);
      tot_invmfp += invmfp;
      acsecs.push_back(invmfp);
      qdeutrons.push_back(nnd);
    } 
    
    // Select quasideuteron interaction from non-zero cross-section choices
    if (verboseLevel > 2) {
      for (std::size_t i=0; i<qdeutrons.size(); i++) {
	G4cout << " acsecs[" << qdeutrons[i].getDefinition()->GetParticleName()
	       << "] " << acsecs[i];
      }
      G4cout << G4endl;
    }
  
    if (tot_invmfp > small) {		// Must have absorption cross-section
      G4double apath = generateInteractionLength(cparticle, path, tot_invmfp);
      
      if (path<small || apath < path) {		// choose the qdeutron
	G4double sl = inuclRndm() * tot_invmfp;
	G4double as = 0.0;
	
	for (std::size_t i = 0; i < qdeutrons.size(); i++) {
	  as += acsecs[i];
	  if (sl < as) { 
	    if (verboseLevel > 2)
	      G4cout << " deut type " << qdeutrons[i] << G4endl; 

	    thePartners.push_back(partner(qdeutrons[i], apath));
	    break;
	  }
	}	// for (qdeutrons...
      }		// if (apath < path)
    }		// if (tot_invmfp > small)
  }		// if (useQuasiDeuteron) [pion, muon or photon]
  
  if (verboseLevel > 2) {
    G4cout << " after deuterons " << thePartners.size() << " partners"
	   << G4endl;
  }
  
  if (thePartners.size() > 1) {		// Sort list by path length
    std::sort(thePartners.begin(), thePartners.end(), sortPartners);
  }
  
  G4InuclElementaryParticle particle;		// Total path at end of list
  thePartners.push_back(partner(particle, path));
}


void G4NucleiModel::
generateParticleFate(G4CascadParticle& cparticle,
		     G4ElementaryParticleCollider* theEPCollider,
		     std::vector<G4CascadParticle>& outgoing_cparticles) {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::generateParticleFate" << G4endl;

  if (verboseLevel > 2) G4cout << " cparticle: " << cparticle << G4endl;

  // Create four-vector checking
#ifdef G4CASCADE_CHECK_ECONS
  G4CascadeCheckBalance balance(0.005, 0.01, "G4NucleiModel");	// Second arg is in GeV
  balance.setVerboseLevel(verboseLevel);
#endif

  outgoing_cparticles.clear();		// Clear return buffer for this event
  generateInteractionPartners(cparticle);	// Fills "thePartners" data

  if (thePartners.empty()) { // smth. is wrong -> needs special treatment
    if (verboseLevel)
      G4cerr << " generateParticleFate-> got empty interaction-partners list "
	     << G4endl;
    return;
  }

  std::size_t npart = thePartners.size();	// Last item is a total-path placeholder

  if (npart == 1) { 		// cparticle is on the next zone entry
    if (verboseLevel > 1) 
      G4cout << " no interactions; moving to next zone" << G4endl;
    
    cparticle.propagateAlongThePath(thePartners[0].second);
    cparticle.incrementCurrentPath(thePartners[0].second);
    boundaryTransition(cparticle);
    outgoing_cparticles.push_back(cparticle);

    if (verboseLevel > 2) G4cout << " next zone \n" << cparticle << G4endl;

    //A.R. 19-Jun-2013: Fixed rare cases of non-reproducibility.
    current_nucl1 = 0;
    current_nucl2 = 0;

    return;
  }	// if (npart == 1)

  if (verboseLevel > 1)
    G4cout << " processing " << npart-1 << " possible interactions" << G4endl;
  
  G4ThreeVector old_position = cparticle.getPosition();
  G4InuclElementaryParticle& bullet = cparticle.getParticle();
  G4bool no_interaction = true;
  G4int zone = cparticle.getCurrentZone();
  
  for (std::size_t i=0; i<npart-1; ++i) {	// Last item is a total-path placeholder
    if (i > 0) cparticle.updatePosition(old_position); 
    
    G4InuclElementaryParticle& target = thePartners[i].first; 
    
    if (verboseLevel > 3) {
      if (target.quasi_deutron()) G4cout << " try absorption: ";
      G4cout << " target " << target.type() << " bullet " << bullet.type()
	     << G4endl;
    }
    
    EPCoutput.reset();
    // Pass current (A,Z) configuration for possible recoils
    G4int massNumberCurrent = protonNumberCurrent+neutronNumberCurrent;
    theEPCollider->setNucleusState(massNumberCurrent, protonNumberCurrent);
    theEPCollider->collide(&bullet, &target, EPCoutput);
    
    // If collision failed, exit loop over partners
    if (EPCoutput.numberOfOutgoingParticles() == 0) break;
    
    if (verboseLevel > 2) {
      EPCoutput.printCollisionOutput();
#ifdef G4CASCADE_CHECK_ECONS
      balance.collide(&bullet, &target, EPCoutput);
      balance.okay();		// Do checks, but ignore result
#endif
    }
    
    // Get list of outgoing particles for evaluation
    std::vector<G4InuclElementaryParticle>& outgoing_particles = 
      EPCoutput.getOutgoingParticles();
    
    if (!passFermi(outgoing_particles, zone)) continue; // Interaction fails
    
    // Trailing effect: reject interaction at previously hit nucleon
    cparticle.propagateAlongThePath(thePartners[i].second);
    const G4ThreeVector& new_position = cparticle.getPosition();
    
    if (!passTrailing(new_position)) continue;
    collisionPts.push_back(new_position);
    
    // Sort particles according to beta (fastest first)
    std::sort(outgoing_particles.begin(), outgoing_particles.end(),
	      G4ParticleLargerBeta() );
    
    if (verboseLevel > 2)
      G4cout << " adding " << outgoing_particles.size()
	     << " output particles" << G4endl;
    
    // NOTE:  Embedded temporary is optimized away (no copying gets done)
    G4int nextGen = cparticle.getGeneration()+1;
    for (std::size_t ip = 0; ip < outgoing_particles.size(); ++ip) { 
      outgoing_cparticles.push_back(G4CascadParticle(outgoing_particles[ip],
						     new_position, zone,
						     0.0, nextGen));
    }
    
    no_interaction = false;
    current_nucl1 = 0;
    current_nucl2 = 0;
    
#ifdef G4CASCADE_DEBUG_CHARGE
    {
      G4double out_charge = 0.0;
      
      for (G4int ip = 0; ip < G4int(outgoing_particles.size()); ip++) 
	out_charge += outgoing_particles[ip].getCharge();
      
      G4cout << " multiplicity " << outgoing_particles.size()
	     << " bul type " << bullet.type()
	     << " targ type " << target.type()
	     << "\n initial charge "
	     << bullet.getCharge() + target.getCharge() 
	     << " out charge " << out_charge << G4endl;  
    }
#endif
    
    if (verboseLevel > 2)
      G4cout << " partner type " << target.type() << G4endl;
    
    if (target.nucleon()) {
      current_nucl1 = target.type();
    } else {
      if (verboseLevel > 2) G4cout << " good absorption " << G4endl;
      current_nucl1 = (target.type() - 100) / 10;
      current_nucl2 = target.type() - 100 - 10 * current_nucl1;
    }   
    
    if (current_nucl1 == 1) {
      if (verboseLevel > 3) G4cout << " decrement proton count" << G4endl;
      protonNumberCurrent--;
    } else {
      if (verboseLevel > 3) G4cout << " decrement neutron count" << G4endl;
      neutronNumberCurrent--;
    } 
    
    if (current_nucl2 == 1) {
      if (verboseLevel > 3) G4cout << " decrement proton count" << G4endl;
      protonNumberCurrent--;
    } else if (current_nucl2 == 2) {
      if (verboseLevel > 3) G4cout << " decrement neutron count" << G4endl;
      neutronNumberCurrent--;
    }
    
    break;
  }		// loop over partners
  
  if (no_interaction) { 		// still no interactions
    if (verboseLevel > 1) G4cout << " no interaction " << G4endl;
    
    // For conservation checking (below), get particle before updating
    static G4ThreadLocal G4InuclElementaryParticle *prescatCP_G4MT_TLS_ = 0;
    if (!prescatCP_G4MT_TLS_) {
      prescatCP_G4MT_TLS_ = new G4InuclElementaryParticle;
      G4AutoDelete::Register(prescatCP_G4MT_TLS_);
    }
    G4InuclElementaryParticle &prescatCP = *prescatCP_G4MT_TLS_;	// Avoid memory churn
    prescatCP = cparticle.getParticle();
    
    // Last "partner" is just a total-path placeholder
    cparticle.updatePosition(old_position); 
    cparticle.propagateAlongThePath(thePartners[npart-1].second);
    cparticle.incrementCurrentPath(thePartners[npart-1].second);
    boundaryTransition(cparticle);
    outgoing_cparticles.push_back(cparticle);
    
    // Check conservation for simple scattering (ignore target nucleus!)
#ifdef G4CASCADE_CHECK_ECONS
    if (verboseLevel > 2) {
      balance.collide(&prescatCP, 0, outgoing_cparticles);
      balance.okay();		// Report violations, but don't act on them
    }
#endif
  }	// if (no_interaction)

  return;
}

// Test if particle is suitable for absorption with dibaryons

G4bool G4NucleiModel::useQuasiDeuteron(G4int ptype, G4int qdtype) {
  if (qdtype==pn || qdtype==0)		// All absorptive particles
    return (ptype==pi0 || ptype==pip || ptype==pim || ptype==gam || ptype==mum);
  else if (qdtype==pp)			// Negative or neutral only
    return (ptype==pi0 || ptype==pim || ptype==gam || ptype==mum);
  else if (qdtype==nn)			// Positive or neutral only
    return (ptype==pi0 || ptype==pip || ptype==gam);

  return false;		// Not a quasideuteron target
}


G4bool G4NucleiModel::passFermi(const std::vector<G4InuclElementaryParticle>& particles, 
				G4int zone) {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::passFermi" << G4endl;
  }

  // Only check Fermi momenta for nucleons
  for (G4int i = 0; i < G4int(particles.size()); i++) {
    if (!particles[i].nucleon()) continue;

    G4int type   = particles[i].type();
    G4double mom = particles[i].getMomModule();
    G4double pfermi  = fermi_momenta[type-1][zone];

    if (verboseLevel > 2)
      G4cout << " type " << type << " p " << mom << " pf " << pfermi << G4endl;
    
    if (mom < pfermi) {
	if (verboseLevel > 2) G4cout << " rejected by Fermi" << G4endl;
	return false;
    }
  }
  return true; 
}


// Test here for trailing effect: loop over all previous collision
// locations and test for d > R_nucleon

G4bool G4NucleiModel::passTrailing(const G4ThreeVector& hit_position) {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::passTrailing " << hit_position << G4endl;

  G4double dist;
  for (G4int i = 0; i < G4int(collisionPts.size() ); i++) {
    dist = (collisionPts[i] - hit_position).mag();
    if (verboseLevel > 2) G4cout << " dist " << dist << G4endl;
    if (dist < R_nucleon) {
      if (verboseLevel > 2) G4cout << " rejected by Trailing" << G4endl;
      return false;
    }
  }
  return true;		// New point far enough away to be used
}


void G4NucleiModel::boundaryTransition(G4CascadParticle& cparticle) {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::boundaryTransition" << G4endl;
  }

  G4int zone = cparticle.getCurrentZone();

  if (cparticle.movingInsideNuclei() && zone == 0) {
    if (verboseLevel) G4cerr << " boundaryTransition-> in zone 0 " << G4endl;
    return;
  }

  G4LorentzVector mom = cparticle.getMomentum();
  G4ThreeVector pos = cparticle.getPosition();
  
  G4int type = cparticle.getParticle().type();
  
  G4double r = pos.mag();
  G4double p = mom.vect().mag();           // NAT
  G4double pr = pos.dot(mom.vect()) / r;
  G4double pperp2 = p*p - pr*pr;           // NAT

  G4int next_zone = cparticle.movingInsideNuclei() ? zone - 1 : zone + 1;

  // NOTE:  dv is the height of the wall seen by the particle  
  G4double dv = getPotential(type,next_zone) - getPotential(type,zone);
  if (verboseLevel > 3) {
    G4cout << "Potentials for type " << type << " = " 
	   << getPotential(type,zone) << " , "
  	   << getPotential(type,next_zone) << G4endl;
  }

  G4double qv = dv*dv + 2.0*dv*mom.e() + pr*pr;
  //  G4double qv = dv*dv + 2.0*dv*mom.m() + pr*pr;    // more correct? NAT
  G4double p1r = 0.; 

  // Perpendicular contribution to pr^2 after penetrating     // NAT
  // potential, to leading order in thickness                 // NAT
  G4double qperp = 2.0*pperp2*potentialThickness/r;           // NAT 

  if (verboseLevel > 3) {
    G4cout << " type " << type << " zone " << zone << " next " << next_zone
	   << " qv " << qv << " dv " << dv << G4endl;
  }

  G4bool adjustpperp = false;                                   // NAT
  G4double smallish = 0.001;                                    // NAT

//  if (qv <= 0.0) { 	// reflection
  if (qv <= 0.0 && qv+qperp <=0 ) {     // reflection         // NAT
    if (verboseLevel > 3) G4cout << " reflects off boundary" << G4endl;
    p1r = -pr;
    cparticle.incrementReflectionCounter();
//  } else {            // transition

  } else if (qv > 0.0) {		// standard transition  // NAT
    if (verboseLevel > 3) G4cout << " passes thru boundary" << G4endl;
    p1r = std::sqrt(qv);
    if (pr < 0.0) p1r = -p1r;
    cparticle.updateZone(next_zone);
    cparticle.resetReflection();

  } else {   // transition via transverse kinetic energy (allowed for thick walls)  // NAT
    if (verboseLevel > 3) G4cout << " passes thru boundary due to angular momentum" << G4endl;
    p1r = smallish * pr; // don't want exactly tangent momentum
    adjustpperp = true;

    cparticle.updateZone(next_zone);
    cparticle.resetReflection();
  }
  
  G4double prr = (p1r - pr)/r;  // change to radial momentum, divided by r
  
  if (verboseLevel > 3) {
    G4cout << " prr " << prr << " delta px " << prr*pos.x() << " py "
	   << prr*pos.y()  << " pz " << prr*pos.z() << " mag "
	   << std::fabs(prr*r) << G4endl;
  }

  if (adjustpperp) {                       // NAT
    G4ThreeVector old_pperp = mom.vect() - pos*(pr/r);
    G4double new_pperp_mag = std::sqrt(std::max(0.0, pperp2 + qv - p1r*p1r) );
    // new total momentum found by rescaling p_perp
    mom.setVect(old_pperp * new_pperp_mag/std::sqrt(pperp2));
    // add a small radial component to make sure that we propagate into new zone
    mom.setVect(mom.vect() + pos*p1r/r);
  } else {
    mom.setVect(mom.vect() + pos*prr);
  }

  cparticle.updateParticleMomentum(mom);
}


// Select random point along full trajectory through nucleus
// NOTE:  Intended for projectile photons for initial interaction

void G4NucleiModel::choosePointAlongTraj(G4CascadParticle& cparticle) {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::choosePointAlongTraj" << G4endl;

  // Get trajectory through nucleus by computing exit point of line,
  // assuming that current position is on surface

  // FIXME:  These need to be reusable (static) buffers
  G4ThreeVector pos  = cparticle.getPosition();
  G4ThreeVector rhat = pos.unit();

  G4ThreeVector phat = cparticle.getMomentum().vect().unit();
  if (cparticle.getMomentum().vect().mag() < small) phat.set(0.,0.,1.);

  if (verboseLevel > 3)
    G4cout << " pos " << pos << " phat " << phat << " rhat " << rhat << G4endl;

  G4ThreeVector posout = pos;
  G4double prang = rhat.angle(-phat);

  if (prang < 1e-6) posout = -pos;		// Check for radial incidence
  else {
    G4double posrot = 2.*prang - pi;
    posout.rotate(posrot, phat.cross(rhat));
    if (verboseLevel > 3) G4cout << " posrot " << posrot/deg << " deg";
  }

  if (verboseLevel > 3) G4cout << " posout " << posout << G4endl;

  // Get list of zone crossings along trajectory
  G4ThreeVector posmid = (pos+posout)/2.;	// Midpoint of trajectory
  G4double r2mid = posmid.mag2();
  G4double lenmid = (posout-pos).mag()/2.;	// Half-length of trajectory

  G4int zoneout = number_of_zones-1;
  G4int zonemid = getZone(posmid.mag());	// Middle zone traversed

  // Every zone is entered then exited, so preallocate vector
  G4int ncross = (number_of_zones-zonemid)*2;

  if (verboseLevel > 3) {
    G4cout << " posmid " << posmid << " lenmid " << lenmid
	   << " zoneout " << zoneout << " zonemid " << zonemid
	   << " ncross " << ncross << G4endl;
  }

  // FIXME:  These need to be reusable (static) buffers
  std::vector<G4double> wtlen(ncross,0.);	// CDF from entry point
  std::vector<G4double> len(ncross,0.);		// Distance from entry point

  // Work from outside in, to accumulate CDF steps properly
  G4int i;				// Loop variable, used multiple times
  for (i=0; i<ncross/2; i++) {
    G4int iz = zoneout-i;
    G4double ds = std::sqrt(zone_radii[iz]*zone_radii[iz]-r2mid);

    len[i] = lenmid - ds;		// Distance from entry to crossing
    len[ncross-1-i] = lenmid + ds;	// Distance to outbound crossing

    if (verboseLevel > 3) {
      G4cout << " i " << i << " iz " << iz << " ds " << ds
	     << " len " << len[i] << G4endl;
    }
  }

  // Compute weights for each zone along trajectory and integrate
  for (i=1; i<ncross; i++) {
    G4int iz = (i<ncross/2) ? zoneout-i+1 : zoneout-ncross+i+1;

    G4double dlen = len[i]-len[i-1];	// Distance from previous crossing

    // Uniform probability across entire zone
    //*** G4double wt = dlen*(getDensity(1,iz)+getDensity(2,iz));

    // Probability based on interaction length through zone
    G4double invmfp = (inverseMeanFreePath(cparticle, neutronEP, iz)
		       + inverseMeanFreePath(cparticle, protonEP, iz));

    // Integral of exp(-len/mfp) from start of zone to end
    G4double wt = (G4Exp(-len[i-1]*invmfp)-G4Exp(-len[i]*invmfp)) / invmfp;

    wtlen[i] = wtlen[i-1] + wt;

    if (verboseLevel > 3) {
      G4cout << " i " << i << " iz " << iz << " avg.mfp " << 1./invmfp
	     << " dlen " << dlen  << " wt " << wt << " wtlen " << wtlen[i]
	     << G4endl;
    }
  }

  // Normalize CDF to unit integral
  std::transform(wtlen.begin(), wtlen.end(), wtlen.begin(),
		 std::bind(std::divides<G4double>(), std::placeholders::_1, wtlen.back()));
  
  if (verboseLevel > 3) {
    G4cout << " weights";
    for (i=0; i<ncross; i++) G4cout << " " << wtlen[i];
    G4cout << G4endl;
  }
  
  // Choose random point along trajectory, weighted by density
  G4double rand = G4UniformRand();
  G4long ir = std::upper_bound(wtlen.begin(),wtlen.end(),rand) - wtlen.begin();

  G4double frac = (rand-wtlen[ir-1]) / (wtlen[ir]-wtlen[ir-1]);
  G4double drand = (1.-frac)*len[ir-1] + frac*len[ir];

  if (verboseLevel > 3) {
    G4cout << " rand " << rand << " ir " << ir << " frac " << frac
	   << " drand " << drand << G4endl;
  }

  pos += drand * phat;		// Distance from entry point along trajectory

  // Update input particle with new location
  cparticle.updatePosition(pos);
  cparticle.updateZone(getZone(pos.mag()));

  if (verboseLevel > 2) {
    G4cout << " moved particle to zone " << cparticle.getCurrentZone() 
	   << " @ " << pos << G4endl;
  }
}


// Returns true if particle should interact immediately
G4bool G4NucleiModel::forceFirst(const G4CascadParticle& cparticle) const {
  return (isProjectile(cparticle) && 
	  (cparticle.getParticle().isPhoton() ||
	   cparticle.getParticle().isMuon())
	  );
}

G4bool G4NucleiModel::isProjectile(const G4CascadParticle& cparticle) const {
  return (cparticle.getGeneration() == 0);	// Only initial state particles
}

G4bool G4NucleiModel::worthToPropagate(const G4CascadParticle& cparticle) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::worthToPropagate" << G4endl;
  }

  const G4double ekin_scale = 2.0;

  G4bool worth = true;

  if (cparticle.reflectedNow()) {	// Just reflected -- keep going?
    G4int zone = cparticle.getCurrentZone();
    G4int ip = cparticle.getParticle().type();

    // NOTE:  Temporarily backing out use of potential for non-nucleons
    G4double ekin_cut = (cparticle.getParticle().nucleon()) ?
      getFermiKinetic(ip, zone) : 0.; //*** getPotential(ip, zone);

    worth = cparticle.getParticle().getKineticEnergy()/ekin_scale > ekin_cut;

    if (verboseLevel > 3) {
      G4cout << " type=" << ip
	     << " ekin=" << cparticle.getParticle().getKineticEnergy()
	     << " potential=" << ekin_cut
	     << " : worth? " << worth << G4endl;
    }
  }

  return worth;
}


G4double G4NucleiModel::getRatio(G4int ip) const {
  if (verboseLevel > 4) {
    G4cout << " >>> G4NucleiModel::getRatio " << ip << G4endl;
  }

  switch (ip) {
  case proton:    return G4double(protonNumberCurrent)/G4double(protonNumber);
  case neutron:   return G4double(neutronNumberCurrent)/G4double(neutronNumber);
  case diproton:  return getRatio(proton)*getRatio(proton);
  case unboundPN: return getRatio(proton)*getRatio(neutron);
  case dineutron: return getRatio(neutron)*getRatio(neutron);
  default:        return 0.;
  }

  return 0.;
}

G4double G4NucleiModel::getCurrentDensity(G4int ip, G4int izone) const {
  const G4double pn_spec = 1.0;		// Scale factor for pn vs. pp/nn
  //const G4double pn_spec = 0.5;

  G4double dens = 0.;

  if (ip < 100) dens = getDensity(ip,izone);
  else {	// For dibaryons, remove extra 1/volume term in density product
    switch (ip) {
    case diproton:  
      dens = getDensity(proton,izone) * getDensity(proton,izone);
      break;
    case unboundPN: 
      dens = getDensity(proton,izone) * getDensity(neutron,izone) * pn_spec;
      break;
    case dineutron:
      dens = getDensity(neutron,izone) * getDensity(neutron,izone);
      break;
    default: dens = 0.;
    }
    dens *= getVolume(izone);
  }

  return getRatio(ip) * dens;
}


G4CascadParticle 
G4NucleiModel::initializeCascad(G4InuclElementaryParticle* particle) {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::initializeCascad(particle)" << G4endl;
  }

  // FIXME:  Previous version generated random sin(theta), then used -cos(theta)
  //         Using generateWithRandomAngles changes result!
  // G4ThreeVector pos = generateWithRandomAngles(nuclei_radius).vect();
  G4double costh = std::sqrt(1.0 - inuclRndm());
  G4ThreeVector pos = generateWithFixedTheta(-costh, nuclei_radius);

  // Start particle outside nucleus, unless capture-at-rest
  G4int zone = number_of_zones;
  if (particle->getKineticEnergy() < small) zone--;

  G4CascadParticle cpart(*particle, pos, zone, large, 0);

  // SPECIAL CASE:  Inbound photons are emplanted along through-path
  if (forceFirst(cpart)) choosePointAlongTraj(cpart);

  if (verboseLevel > 2) G4cout << cpart << G4endl;

  return cpart;
}

void G4NucleiModel::initializeCascad(G4InuclNuclei* bullet, 
				     G4InuclNuclei* target,
				     modelLists& output) {
  if (verboseLevel) {
    G4cout << " >>> G4NucleiModel::initializeCascad(bullet,target,output)"
	   << G4endl;
  }

  const G4double max_a_for_cascad = 5.0;
  const G4double ekin_cut = 2.0;
  const G4double small_ekin = 1.0e-6;
  const G4double r_large2for3 = 62.0;
  const G4double r0forAeq3 = 3.92;
  const G4double s3max = 6.5;
  const G4double r_large2for4 = 69.14;
  const G4double r0forAeq4 = 4.16;
  const G4double s4max = 7.0;
  const G4int itry_max = 100;

  // Convenient references to input buffer contents
  std::vector<G4CascadParticle>& casparticles = output.first;
  std::vector<G4InuclElementaryParticle>& particles = output.second;

  casparticles.clear();		// Reset input buffer to fill this event
  particles.clear();

  // first decide whether it will be cascad or compound final nuclei
  G4int ab = bullet->getA();
  G4int zb = bullet->getZ();
  G4int at = target->getA();
  G4int zt = target->getZ();

  G4double massb = bullet->getMass();	// For creating LorentzVectors below

  if (ab < max_a_for_cascad) {

    G4double benb = bindingEnergy(ab,zb)/GeV / G4double(ab);
    G4double bent = bindingEnergy(at,zt)/GeV / G4double(at);
    G4double ben = benb < bent ? bent : benb;

    if (bullet->getKineticEnergy()/ab > ekin_cut*ben) {
      G4int itryg = 0;

	/* Loop checking 08.06.2015 MHK */
      while (casparticles.size() == 0 && itryg < itry_max) {      
	itryg++;
	particles.clear();
      
	//    nucleons coordinates and momenta in nuclei rest frame
	coordinates.clear();
	momentums.clear();
     
	if (ab < 3) { // deuteron, simplest case
	  G4double r = 2.214 - 3.4208 * G4Log(1.0 - 0.981 * inuclRndm());
	  G4ThreeVector coord1 = generateWithRandomAngles(r).vect();
	  coordinates.push_back(coord1);
	  coordinates.push_back(-coord1);

	  G4double p = 0.0;
	  G4bool bad = true;
	  G4int itry = 0;

	  while (bad && itry < itry_max) {  /* Loop checking 08.06.2015 MHK */
	    itry++;
	    p = 456.0 * inuclRndm();

	    if (p * p / (p * p + 2079.36) / (p * p + 2079.36) > 1.2023e-4 * inuclRndm() &&
	       p * r > 312.0) bad = false;
	  }

	  if (itry == itry_max)
	    if (verboseLevel > 2){ 
	      G4cout << " deutron bullet generation-> itry = " << itry_max << G4endl;	
	    }

	  p = 0.0005 * p;

	  if (verboseLevel > 2){ 
	    G4cout << " p nuc " << p << G4endl;
	  }

	  G4LorentzVector mom = generateWithRandomAngles(p, massb);

	  momentums.push_back(mom);
	  mom.setVect(-mom.vect());
	  momentums.push_back(-mom);
	} else {
	  G4ThreeVector coord1;

	  G4bool badco = true;

	  G4int itry = 0;
        
	  if (ab == 3) {
	    while (badco && itry < itry_max) {/* Loop checking 08.06.2015 MHK */
	      if (itry > 0) coordinates.clear();
	      itry++;	
	      G4int i(0);    

	      for (i = 0; i < 2; i++) {
		G4int itry1 = 0;
		G4double ss, u, rho; 
		G4double fmax = G4Exp(-0.5) / std::sqrt(0.5);

		while (itry1 < itry_max) {  /* Loop checking 08.06.2015 MHK */
		  itry1++;
		  ss = -G4Log(inuclRndm());
		  u = fmax * inuclRndm();
		  rho = std::sqrt(ss) * G4Exp(-ss);

		  if (rho > u && ss < s3max) {
		    ss = r0forAeq3 * std::sqrt(ss);
		    coord1 = generateWithRandomAngles(ss).vect();
		    coordinates.push_back(coord1);

		    if (verboseLevel > 2){
		      G4cout << " i " << i << " r " << coord1.mag() << G4endl;
		    }
		    break;
		  }
		}

		if (itry1 == itry_max) { // bad case
		  coord1.set(10000.,10000.,10000.);
		  coordinates.push_back(coord1);
		  break;
		}
	      }

	      coord1 = -coordinates[0] - coordinates[1]; 
	      if (verboseLevel > 2) {
		G4cout << " 3  r " << coord1.mag() << G4endl;
	      }

	      coordinates.push_back(coord1);  	    
	    
	      G4bool large_dist = false;

	      for (i = 0; i < 2; i++) {
		for (G4int j = i+1; j < 3; j++) {
		  G4double r2 = (coordinates[i]-coordinates[j]).mag2();

		  if (verboseLevel > 2) {
		    G4cout << " i " << i << " j " << j << " r2 " << r2 << G4endl;
		  }

		  if (r2 > r_large2for3) {
		    large_dist = true;

		    break; 
		  }      
		}

		if (large_dist) break;
	      } 

	      if(!large_dist) badco = false;

	    }

	  } else { // a >= 4
	    G4double b = 3./(ab - 2.0);
	    G4double b1 = 1.0 - b / 2.0;
	    G4double u = b1 + std::sqrt(b1 * b1 + b);
	    G4double fmax = (1.0 + u/b) * u * G4Exp(-u);
	  
	    while (badco && itry < itry_max) {/* Loop checking 08.06.2015 MHK */

	      if (itry > 0) coordinates.clear();
	      itry++;
	      G4int i(0);
	    
	      for (i = 0; i < ab-1; i++) {
		G4int itry1 = 0;
		G4double ss; 

		while (itry1 < itry_max) {  /* Loop checking 08.06.2015 MHK */
		  itry1++;
		  ss = -G4Log(inuclRndm());
		  u = fmax * inuclRndm();

		  if (std::sqrt(ss) * G4Exp(-ss) * (1.0 + ss/b) > u
		      && ss < s4max) {
		    ss = r0forAeq4 * std::sqrt(ss);
		    coord1 = generateWithRandomAngles(ss).vect();
		    coordinates.push_back(coord1);

		    if (verboseLevel > 2) {
		      G4cout << " i " << i << " r " << coord1.mag() << G4endl;
		    }

		    break;
		  }
		}

		if (itry1 == itry_max) { // bad case
		  coord1.set(10000.,10000.,10000.);
		  coordinates.push_back(coord1);
		  break;
		}
	      }

	      coord1 *= 0.0;	// Cheap way to reset
	      for(G4int j = 0; j < ab -1; j++) coord1 -= coordinates[j];

	      coordinates.push_back(coord1);   

	      if (verboseLevel > 2){
		G4cout << " last r " << coord1.mag() << G4endl;
	      }
	    
	      G4bool large_dist = false;

	      for (i = 0; i < ab-1; i++) {
		for (G4int j = i+1; j < ab; j++) {
	     
		  G4double r2 = (coordinates[i]-coordinates[j]).mag2();

		  if (verboseLevel > 2){
		    G4cout << " i " << i << " j " << j << " r2 " << r2 << G4endl;
		  }

		  if (r2 > r_large2for4) {
		    large_dist = true;

		    break; 
		  }      
		}

		if (large_dist) break;
	      } 

	      if (!large_dist) badco = false;
	    }
	  } 

	  if(badco) {
	    G4cout << " can not generate the nucleons coordinates for a "
		   << ab << G4endl;	

	    casparticles.clear();	// Return empty buffer on error
	    particles.clear();
	    return;

	  } else { // momentums
	    G4double p, u, x;
	    G4LorentzVector mom;
	    //G4bool badp = True;

	    for (G4int i = 0; i < ab - 1; i++) {
	      G4int itry2 = 0;

	      while(itry2 < itry_max) {	/* Loop checking 08.06.2015 MHK */
		itry2++;
		u = -G4Log(0.879853 - 0.8798502 * inuclRndm());
		x = u * G4Exp(-u);

		if(x > inuclRndm()) {
		  p = std::sqrt(0.01953 * u);
		  mom = generateWithRandomAngles(p, massb);
		  momentums.push_back(mom);

		  break;
		}
	      }	// while (itry2 < itry_max)

	      if(itry2 == itry_max) {
		G4cout << " can not generate proper momentum for a "
		       << ab << G4endl;

		casparticles.clear();	// Return empty buffer on error
		particles.clear();
		return;
	      } 
	    }	// for (i=0 ...
	    // last momentum

	    mom *= 0.;		// Cheap way to reset
	    mom.setE(bullet->getEnergy()+target->getEnergy());

	    for(G4int j=0; j< ab-1; j++) mom -= momentums[j]; 

	    momentums.push_back(mom);
	  } 
	}
 
	// Coordinates and momenta at rest are generated, now back to the lab
	G4double rb = 0.0;
	G4int i(0);

	for(i = 0; i < G4int(coordinates.size()); i++) {      
	  G4double rp = coordinates[i].mag2();

	  if(rp > rb) rb = rp;
	}

	// nuclei i.p. as a whole
	G4double s1 = std::sqrt(inuclRndm()); 
	G4double phi = randomPHI();
	G4double rz = (nuclei_radius + rb) * s1;
	G4ThreeVector global_pos(rz*std::cos(phi), rz*std::sin(phi),
				 -(nuclei_radius+rb)*std::sqrt(1.0-s1*s1));

	for (i = 0; i < G4int(coordinates.size()); i++) {
	  coordinates[i] += global_pos;
	}  

	// all nucleons at rest
	raw_particles.clear();

	for (G4int ipa = 0; ipa < ab; ipa++) {
	  G4int knd = ipa < zb ? 1 : 2;
	  raw_particles.push_back(G4InuclElementaryParticle(momentums[ipa], knd));
	} 
      
	G4InuclElementaryParticle dummy(small_ekin, 1);
	G4LorentzConvertor toTheBulletRestFrame(&dummy, bullet);
	toTheBulletRestFrame.toTheTargetRestFrame();

	particleIterator ipart;

	for (ipart = raw_particles.begin(); ipart != raw_particles.end(); ipart++) {
	  ipart->setMomentum(toTheBulletRestFrame.backToTheLab(ipart->getMomentum())); 
	}

	// fill cascad particles and outgoing particles

	for(G4int ip = 0; ip < G4int(raw_particles.size()); ip++) {
	  G4LorentzVector mom = raw_particles[ip].getMomentum();
	  G4double pmod = mom.rho();
	  G4double t0 = -mom.vect().dot(coordinates[ip]) / pmod;
	  G4double det = t0 * t0 + nuclei_radius * nuclei_radius
	               - coordinates[ip].mag2();
	  G4double tr = -1.0;

	  if(det > 0.0) {
	    G4double t1 = t0 + std::sqrt(det);
	    G4double t2 = t0 - std::sqrt(det);

	    if(std::fabs(t1) <= std::fabs(t2)) {	 
	      if(t1 > 0.0) {
		if(coordinates[ip].z() + mom.z() * t1 / pmod <= 0.0) tr = t1;
	      }

	      if(tr < 0.0 && t2 > 0.0) {

		if(coordinates[ip].z() + mom.z() * t2 / pmod <= 0.0) tr = t2;
	      }

	    } else {
	      if(t2 > 0.0) {

		if(coordinates[ip].z() + mom.z() * t2 / pmod <= 0.0) tr = t2;
	      }

	      if(tr < 0.0 && t1 > 0.0) {
		if(coordinates[ip].z() + mom.z() * t1 / pmod <= 0.0) tr = t1;
	      }
	    } 

	  }

	  if(tr >= 0.0) { // cascad particle
	    coordinates[ip] += mom.vect()*tr / pmod;
	    casparticles.push_back(G4CascadParticle(raw_particles[ip], 
						    coordinates[ip], 
						    number_of_zones, large, 0));

	  } else {
	    particles.push_back(raw_particles[ip]); 
	  } 
	}
      }    

      if(casparticles.size() == 0) {
	particles.clear();

	G4cout << " can not generate proper distribution for " << itry_max
	       << " steps " << G4endl;
      }    
    }
  }

  if(verboseLevel > 2){
    G4int ip(0);

    G4cout << " cascad particles: " << casparticles.size() << G4endl;
    for(ip = 0; ip < G4int(casparticles.size()); ip++)
      G4cout << casparticles[ip] << G4endl;

    G4cout << " outgoing particles: " << particles.size() << G4endl;
    for(ip = 0; ip < G4int(particles.size()); ip++)
      G4cout << particles[ip] << G4endl;
  }

  return;	// Buffer has been filled
}


// Compute mean free path for inclusive interaction of projectile and target
G4double 
G4NucleiModel::inverseMeanFreePath(const G4CascadParticle& cparticle,
				   const G4InuclElementaryParticle& target,
				   G4int zone) {
  G4int ptype = cparticle.getParticle().type();
  G4int ip = target.type();

  // Ensure that zone specified is within nucleus, for array lookups
  if (zone<0) zone = cparticle.getCurrentZone();
  if (zone>=number_of_zones) zone = number_of_zones-1;

  // Special cases: neutrinos, and muon-on-neutron, have infinite path
  if (cparticle.getParticle().isNeutrino()) return 0.;
  if (ptype == muonMinus && ip == neutron) return 0.;

  // Configure frame transformation to get kinetic energy for lookups
  dummy_convertor.setBullet(cparticle.getParticle());
  dummy_convertor.setTarget(&target);
  dummy_convertor.toTheCenterOfMass();		// Fill internal kinematics
  G4double ekin = dummy_convertor.getKinEnergyInTheTRS();

  // Get cross section for interacting with target (dibaryons are absorptive)
  G4double csec = (ip < 100) ? totalCrossSection(ekin, ptype*ip)
                             : absorptionCrossSection(ekin, ptype);

  if (verboseLevel > 2) {
    G4cout << " ip " << ip << " zone " << zone << " ekin " << ekin
	   << " dens " << getCurrentDensity(ip, zone)
	   << " csec " << csec << G4endl;
  }

  if (csec <= 0.) return 0.;	// No interaction, avoid unnecessary work

  // Get nuclear density and compute mean free path
  return csec * getCurrentDensity(ip, zone);
}

// Throw random distance for interaction of particle using path and MFP

G4double 
G4NucleiModel::generateInteractionLength(const G4CascadParticle& cparticle,
					 G4double path, G4double invmfp) const {
  // Delay interactions of newly formed secondaries (minimum int. length)
  const G4double young_cut = std::sqrt(10.0) * 0.25;
  const G4double huge_num = 50.0;	// Argument to exponential

  G4double spath = large;		// Buffer for return value

  if (invmfp < small) return spath;	// No interaction, avoid unnecessary work

  G4double pw = -path * invmfp;		// Ratio of path in zone to MFP
  if (pw < -huge_num) pw = -huge_num;
  pw = 1.0 - G4Exp(pw);
    
  if (verboseLevel > 2) 
    G4cout << " mfp " << 1./invmfp << " pw " << pw << G4endl;
    
  // Primary particle(s) should always interact at least once
  if (forceFirst(cparticle) || (inuclRndm() < pw)) {
    spath = -G4Log(1.0 - pw * inuclRndm()) / invmfp;
    if (cparticle.young(young_cut, spath)) spath = large;
    
    if (verboseLevel > 2) 
      G4cout << " spath " << spath << " path " << path << G4endl;
  }

  return spath;
}


// Parametrized cross section for pion and photon absorption

G4double G4NucleiModel::absorptionCrossSection(G4double ke, G4int type) const {
  if (!useQuasiDeuteron(type)) {
    G4cerr << "absorptionCrossSection() only valid for incident pions or gammas" 
           << G4endl;
    return 0.;
  }

  G4double csec = 0.;

  // Pion absorption is parametrized for low vs. medium energy
  // ... use for muon capture as well
  if (type == pionPlus || type == pionMinus || type == pionZero ||
      type == muonMinus) {
    if (ke < 0.3) csec = (0.1106 / std::sqrt(ke) - 0.8
			  + 0.08 / ((ke-0.123)*(ke-0.123) + 0.0056) );
    else if (ke < 1.0) csec = 3.6735 * (1.0-ke)*(1.0-ke);     
  }

  if (type == photon) {
    csec = gammaQDinterp.interpolate(ke, gammaQDxsec) * gammaQDscale;
  }

  if (csec < 0.0) csec = 0.0;

  if (verboseLevel > 2) {
    G4cout << " ekin " << ke << " abs. csec " << csec << " mb" << G4endl;   
  }

  return crossSectionUnits * csec;
}


G4double G4NucleiModel::totalCrossSection(G4double ke, G4int rtype) const
{
  // All scattering cross-sections are available from tables
  const G4CascadeChannel* xsecTable = G4CascadeChannelTables::GetTable(rtype);
  if (!xsecTable) {
    G4cerr << " unknown collison type = " << rtype << G4endl;
    return 0.;
  }

  return (crossSectionUnits * xsecTable->getCrossSection(ke));
}
