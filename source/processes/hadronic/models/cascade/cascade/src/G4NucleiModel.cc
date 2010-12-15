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
// $Id: G4NucleiModel.cc,v 1.96 2010-12-15 07:41:17 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
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
// 20101119  M. Kelsey -- Hide "negative path" and "no partners" messages in
//		verbosity.

#include "G4NucleiModel.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CascadeInterpolator.hh"
#include "G4CascadeNNChannel.hh"
#include "G4CascadeNPChannel.hh"
#include "G4CascadePPChannel.hh"
#include "G4CascadePiMinusNChannel.hh"
#include "G4CascadePiMinusPChannel.hh"
#include "G4CascadePiPlusNChannel.hh"
#include "G4CascadePiPlusPChannel.hh"
#include "G4CascadePiZeroNChannel.hh"
#include "G4CascadePiZeroPChannel.hh"
#include "G4CollisionOutput.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4HadTmpUtil.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticleNames.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4LorentzConvertor.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"

using namespace G4InuclParticleNames;
using namespace G4InuclSpecialFunctions;

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

// Parameters for nuclear structure (actual vs. "should be" in fm and GeV)
const G4double G4NucleiModel::skinDepth = 1.7234;	// 0.61 fm
const G4double G4NucleiModel::radiusScale = 3.3836;	// 1.2 fm [*cbrt(A)]
const G4double G4NucleiModel::radiusForSmall = 8.0;	// R_p = 0.8768 fm
const G4double G4NucleiModel::fermiMomentum = 1.932;	// 0.61 GeV fm

// Zone boundaries as fraction of nuclear radius (from outside in)
const G4double G4NucleiModel::alfa3[3] = { 0.7, 0.3, 0.01 };
const G4double G4NucleiModel::alfa6[6] = { 0.9, 0.6, 0.4, 0.2, 0.1, 0.05 };

// Flat nuclear potentials for mesons and hyperons (GeV)
const G4double G4NucleiModel::pion_vp = 0.007;
const G4double G4NucleiModel::pion_vp_small = 0.007;
const G4double G4NucleiModel::kaon_vp = -0.015;		// WHY NEGATIVE?
const G4double G4NucleiModel::hyperon_vp = 0.030;

// FIXME:  We should not be using this!
const G4double G4NucleiModel::piTimes4thirds = 4.189; // 4.188790204786;


// Constructors
G4NucleiModel::G4NucleiModel()
  : verboseLevel(0), nuclei_radius(0.), number_of_zones(0),
    A(0), Z(0), theNucleus(0), neutronNumber(0), protonNumber(0),
    neutronNumberCurrent(0), protonNumberCurrent(0), current_nucl1(0),
    current_nucl2(0) {}

G4NucleiModel::G4NucleiModel(G4int a, G4int z)
  : verboseLevel(0), nuclei_radius(0.), number_of_zones(0),
    A(0), Z(0), theNucleus(0), neutronNumber(0), protonNumber(0),
    neutronNumberCurrent(0), protonNumberCurrent(0), current_nucl1(0),
    current_nucl2(0) {
  generateModel(a,z);
}

G4NucleiModel::G4NucleiModel(G4InuclNuclei* nuclei)
  : verboseLevel(0), nuclei_radius(0.), number_of_zones(0),
    A(0), Z(0), theNucleus(0), neutronNumber(0), protonNumber(0),
    neutronNumberCurrent(0), protonNumberCurrent(0), current_nucl1(0),
    current_nucl2(0) {
  generateModel(nuclei);
}

G4NucleiModel::~G4NucleiModel() {
  delete theNucleus;
  theNucleus = 0;
}


void G4NucleiModel::generateModel(G4InuclNuclei* nuclei) {
  generateModel(nuclei->getA(), nuclei->getZ());
}

void G4NucleiModel::generateModel(G4int a, G4int z) {
  if (verboseLevel) {
    G4cout << " >>> G4NucleiModel::generateModel A " << a << " Z " << z
	   << G4endl;
  }

  // If model already built, just return; otherwise intialize everything
  if (a == A && z == Z) {
    if (verboseLevel > 1)
      G4cout << " model already generated for A=" << a << ", Z=" << z << G4endl;

    reset();		// Zeros out neutron/proton evaporates
    return;
  }

  A = a;
  Z = z;
  delete theNucleus;
  theNucleus = new G4InuclNuclei(A,Z);		// For conservation checking

  neutronNumber = A - Z;
  protonNumber = Z;
  reset();

  G4double nuclearRadius = radiusScale*G4cbrt(A);	// Nuclear radius

  // This will be used to pre-allocate lots of arrays below
  number_of_zones = (A < 5) ? 1 : (A < 100) ? 3 : 6;

  // Clear all parameters arrays for reloading
  binding_energies.clear();
  nucleon_densities.clear();
  zone_potentials.clear();
  fermi_momenta.clear();
  zone_radii.clear();

  fillBindingEnergies();
  fillZoneRadii(nuclearRadius);

  G4double tot_vol = fillZoneVolumes(nuclearRadius);

  fillPotentials(proton, tot_vol);		// Protons
  fillPotentials(neutron, tot_vol);		// Neutrons

  // Additional flat zone potentials for other hadrons
  const std::vector<G4double> vp(number_of_zones, (A>4)?pion_vp:pion_vp_small);
  const std::vector<G4double> kp(number_of_zones, kaon_vp);
  const std::vector<G4double> hp(number_of_zones, hyperon_vp);

  zone_potentials.push_back(vp);
  zone_potentials.push_back(kp);
  zone_potentials.push_back(hp);

  nuclei_radius = zone_radii[number_of_zones-1];
}


// Load binding energy array for current nucleus

void G4NucleiModel::fillBindingEnergies() {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::fillBindingEnergies" << G4endl;

  G4double dm = bindingEnergy(A,Z);

  // Binding energy differences for proton and neutron loss, respectively
  binding_energies.push_back(std::fabs(bindingEnergy(A-1,Z-1)-dm)/GeV);
  binding_energies.push_back(std::fabs(bindingEnergy(A-1,Z)-dm)/GeV);
}

// Load zone boundary radius array for current nucleus
void 
G4NucleiModel::fillZoneRadii(G4double nuclearRadius) {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::fillZoneRadii" << G4endl;

  G4double skinRatio = nuclearRadius/skinDepth;
  G4double skinDecay = std::exp(-skinRatio);    

  if (A < 5) {			// Light ions treated as simple balls
    G4double smallRad = radiusForSmall;
    if (A == 4) smallRad *= 0.7;	// Alpha is most tightly bound
    zone_radii.push_back(smallRad);
    ur[0] = 0.;
    ur[1] = 1.;
  } else if (A < 12) {		// Small nuclei have Gaussian potential
    G4double rSq = nuclearRadius * nuclearRadius;
    G4double gaussRadius = std::sqrt(rSq * (1.0 - 1.0/A) + 6.4);
    
    ur[0] = 0.0;
    for (G4int i = 0; i < number_of_zones; i++) {
      G4double y = std::sqrt(-std::log(alfa3[i]));
      zone_radii.push_back(gaussRadius * y);
      ur[i+1] = y;
    }
  } else if (A < 100) {		// Intermediate nuclei have three-zone W-S
    ur[0] = -skinRatio;
    for (G4int i = 0; i < number_of_zones; i++) {
      G4double y = std::log((1.0 + skinDecay)/alfa3[i] - 1.0);
      zone_radii.push_back(nuclearRadius + skinDepth * y);
      ur[i+1] = y;
    }
  } else {			// Heavy nuclei have six-zone W-S
    ur[0] = -skinRatio;
    for (G4int i = 0; i < number_of_zones; i++) {
      G4double y = std::log((1.0 + skinDecay)/alfa6[i] - 1.0);
      zone_radii.push_back(nuclearRadius + skinDepth * y);
      ur[i+1] = y;
    }
  }
}

// Compute zone-by-zone density integrals
G4double 
G4NucleiModel::fillZoneVolumes(G4double nuclearRadius) {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::fillZoneVolumes" << G4endl;

  G4double tot_vol = 0.;	// Return value

  if (A < 5) {			// Light ions are treated as simple balls
    v[0] = v1[0] = 1.;
    tot_vol = zone_radii[0]*zone_radii[0]*zone_radii[0];
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
  }

  return tot_vol;
}

// Load potentials for nucleons (using G4InuclParticle codes)
void G4NucleiModel::fillPotentials(G4int type, G4double tot_vol) {
  if (verboseLevel > 1) 
    G4cout << " >>> G4NucleiModel::fillZoneVolumes(" << type << ")" << G4endl;

  if (type != proton && type != neutron) return;

  const G4double mass = G4InuclElementaryParticle::getParticleMass(type);
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
G4double
G4NucleiModel::zoneIntegralWoodsSaxon(G4double r1, G4double r2, 
				      G4double nuclearRadius) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::zoneIntegralWoodsSaxon" << G4endl;
  }

  const G4double epsilon = 1.0e-3;
  const G4int itry_max = 1000;

  G4double skinRatio = nuclearRadius / skinDepth;

  G4double d2 = 2.0 * skinRatio;
  G4double dr = r2 - r1;
  G4double fr1 = r1 * (r1 + d2) / (1.0 + std::exp(r1));
  G4double fr2 = r2 * (r2 + d2) / (1.0 + std::exp(r2));
  G4double fi = (fr1 + fr2) / 2.;
  G4double fun1 = fi * dr;
  G4double fun;
  G4double jc = 1;
  G4double dr1 = dr;
  G4int itry = 0;

  while (itry < itry_max) {
    dr /= 2.;
    itry++;

    G4double r = r1 - dr;
    fi = 0.0;
    G4int jc1 = G4int(std::pow(2.0, jc - 1) + 0.1);

    for (G4int i = 0; i < jc1; i++) { 
      r += dr1; 
      fi += r * (r + d2) / (1.0 + std::exp(r));
    }

    fun = 0.5 * fun1 + fi * dr;

    if (std::fabs((fun - fun1) / fun) <= epsilon) break;

    jc++;
    dr1 = dr;
    fun1 = fun;
  }	// while (itry < itry_max)

  if (verboseLevel > 2 && itry == itry_max)
    G4cout << " zoneIntegralWoodsSaxon-> n iter " << itry_max << G4endl;

  G4double skinDepth3 = 5.11864; //*** skinDepth*skinDepth*skinDepth;

  return skinDepth3 * (fun + skinRatio*skinRatio*std::log((1.0 + std::exp(-r1)) / (1.0 + std::exp(-r2))));
}


// Zone integral of Gaussian density function
G4double
G4NucleiModel::zoneIntegralGaussian(G4double r1, G4double r2, 
				    G4double nucRad) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::zoneIntegralGaussian" << G4endl;
  }

  G4double gaussRadius = std::sqrt(nucRad*nucRad * (1.0 - 1.0/A) + 6.4);

  const G4double epsilon = 1.0e-3;
  const G4int itry_max = 1000;

  G4double dr = r2 - r1;
  G4double fr1 = r1 * r1 * std::exp(-r1 * r1);
  G4double fr2 = r2 * r2 * std::exp(-r2 * r2);
  G4double fi = (fr1 + fr2) / 2.;
  G4double fun1 = fi * dr;
  G4double fun;
  G4double jc = 1;
  G4double dr1 = dr;
  G4int itry = 0;

  while (itry < itry_max) {
    dr /= 2.;
    itry++;
    G4double r = r1 - dr;
    fi = 0.0;
    G4int jc1 = int(std::pow(2.0, jc - 1) + 0.1);

    for (G4int i = 0; i < jc1; i++) { 
      r += dr1; 
      fi += r * r * std::exp(-r * r);
    }

    fun = 0.5 * fun1 + fi * dr;  

    if (std::fabs((fun - fun1) / fun) <= epsilon) break;

    jc++;
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
	 << " Nuclei radius " << nuclei_radius << " number of zones "
	 << number_of_zones << G4endl;

  for (G4int i = 0; i < number_of_zones; i++)
    G4cout << " zone " << i+1 << " radius " << zone_radii[i] << G4endl
	   << " protons: density " << getDensity(1,i) << " PF " << 
      getFermiMomentum(1,i) << " VP " << getPotential(1,i) << G4endl
	   << " neutrons: density " << getDensity(2,i) << " PF " << 
      getFermiMomentum(2,i) << " VP " << getPotential(2,i) << G4endl
	   << " pions: VP " << getPotential(3,i) << G4endl;
}


G4double G4NucleiModel::getFermiKinetic(G4int ip, G4int izone) const {
  G4double ekin = 0.0;
  
  if (ip < 3 && izone < number_of_zones) {	// ip for proton/neutron only
    G4double pf = fermi_momenta[ip - 1][izone]; 
    G4double mass = G4InuclElementaryParticle::getParticleMass(ip);
    
    ekin = std::sqrt(pf * pf + mass * mass) - mass;
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
G4NucleiModel::generateQuasiDeutron(G4int type1, G4int type2,
				    G4int zone) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::generateQuasiDeutron" << G4endl;
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

  const G4double small = 1.0e-10;
  const G4double huge_num = 50.0;
  const G4double pn_spec = 1.0;

  //const G4double pn_spec = 0.5;

  //const G4double young_cut = std::sqrt(10.0) * 0.5;
  //const G4double young_cut = std::sqrt(10.0) * 0.45;
  const G4double young_cut = std::sqrt(10.0) * 0.25;
  //const G4double young_cut = std::sqrt(10.0) * 0.2;
  //const G4double young_cut = std::sqrt(10.0) * 0.1;
  //const G4double young_cut = 0.0;

  thePartners.clear();		// Reset buffer for next cycle

  G4int ptype = cparticle.getParticle().type();
  G4int zone = cparticle.getCurrentZone();
  G4double pmass = cparticle.getParticle().getMass();
  G4LorentzVector pmom = cparticle.getParticle().getMomentum();
  G4double r_in;
  G4double r_out;

  if (zone == number_of_zones) { // particle is outside 
    r_in = nuclei_radius;
    r_out = 0.0;

  } else if (zone == 0) { // particle is outside core
    r_in = 0.0;
    r_out = zone_radii[0];

  } else {
    r_in = zone_radii[zone - 1];
    r_out = zone_radii[zone];
  }  

  G4double path = cparticle.getPathToTheNextZone(r_in, r_out);

  if (verboseLevel > 2) {
    G4cout << " r_in " << r_in << " r_out " << r_out << " path " << path << G4endl;
  }

  if (path < -small) { 			// something wrong
    if (verboseLevel)
      G4cerr << " generateInteractionPartners-> negative path length" << G4endl;

    return;
  }

  if (std::fabs(path) < small) { 	// just on the boundary
    if (verboseLevel > 3)
      G4cout << " generateInteractionPartners-> zero path" << G4endl;

    thePartners.push_back(partner());	// Dummy list terminator with zero path
    return;
  }

  G4LorentzConvertor dummy_convertor;
  dummy_convertor.setBullet(pmom, pmass);

  for (G4int ip = 1; ip < 3; ip++) {
    // Only process nucleons which remain active in target
    if (ip==proton  && protonNumberCurrent < 1) continue;
    if (ip==neutron && neutronNumberCurrent < 1) continue;

    // All nucleons are assumed to be at rest when colliding
    G4InuclElementaryParticle particle = generateNucleon(ip, zone);
    dummy_convertor.setTarget(particle.getMomentum(), particle.getMass());
    G4double ekin = dummy_convertor.getKinEnergyInTheTRS();
    
    G4double csec = totalCrossSection(ekin, ptype * ip);
    
    if(verboseLevel > 2) {
      G4cout << " ip " << ip << " ekin " << ekin << " csec " << csec << G4endl;
    }
    
    G4double dens = nucleon_densities[ip - 1][zone];
    G4double rat = getRatio(ip);
    G4double pw = -path * dens * csec * rat;
    
    if (pw < -huge_num) pw = -huge_num;
    pw = 1.0 - std::exp(pw);
    
    if (verboseLevel > 2) {
      G4cout << " pw " << pw << " rat " << rat << G4endl;
    }
    
    G4double spath = path;
    
    if (inuclRndm() < pw) {
      spath = -1.0 / dens / csec / rat * std::log(1.0 - pw * inuclRndm());
      if (cparticle.young(young_cut, spath)) spath = path;
      
      if (verboseLevel > 2) {
	G4cout << " ip " << ip << " spath " << spath << G4endl;
      }
    }

    if (spath < path) {
      if (verboseLevel > 3) {
	G4cout << " adding partner[" << thePartners.size() << "]: ";
	particle.printParticle();
      }
      thePartners.push_back(partner(particle, spath));
    }
  }	// for (G4int ip...
  
  if (verboseLevel > 2) {
    G4cout << " after nucleons " << thePartners.size() << " path " << path << G4endl;
  }
  
  if (cparticle.getParticle().pion()) { // absorption possible
    if (verboseLevel > 2) {
      G4cout << " trying quasi-deuterons with bullet: ";
      cparticle.getParticle().printParticle();
    }

    // Initialize buffers for quasi-deuteron results
    qdeutrons.clear();
    acsecs.clear();

    G4double tot_abs_csec = 0.0;
    G4double abs_sec;
    G4double vol = zone_radii[zone]*zone_radii[zone]*zone_radii[zone];
    
    if (zone > 0) vol -= zone_radii[zone-1]*zone_radii[zone-1]*zone_radii[zone-1];
    vol *= piTimes4thirds; 
    
    G4double rat  = getRatio(1); 
    G4double rat1 = getRatio(2); 

    // FIXME:  Shouldn't be creating zero-cross-section states!

    // Proton-proton state interacts with pi-, pi0 only
    G4InuclElementaryParticle ppd = generateQuasiDeutron(1, 1, zone);
    
    if (protonNumberCurrent < 2 || !(ptype == pi0 || ptype == pim)) {
      abs_sec = 0.0;
    } else {
      dummy_convertor.setTarget(ppd.getMomentum(), ppd.getMass());
      
      G4double ekin = dummy_convertor.getKinEnergyInTheTRS();
      
      if (verboseLevel > 2) {
	G4cout << " ptype=" << ptype << " using pp target" << G4endl;
	ppd.printParticle();
      }
      
      abs_sec = absorptionCrossSection(ekin, ptype);
      abs_sec *= nucleon_densities[0][zone] * nucleon_densities[0][zone]*
	rat * rat * vol; 
    }
    
    tot_abs_csec += abs_sec;
    acsecs.push_back(abs_sec);
    qdeutrons.push_back(ppd);
    
    // Proton-neutron state interacts with any pion type
    G4InuclElementaryParticle npd = generateQuasiDeutron(1, 2, zone);

    if (protonNumberCurrent < 1 || neutronNumberCurrent < 1) {
      abs_sec = 0.0;
    } else {
      dummy_convertor.setTarget(npd.getMomentum(), npd.getMass());
      
      G4double ekin = dummy_convertor.getKinEnergyInTheTRS();
      
      if (verboseLevel > 2) {
	G4cout << " using np target" << G4endl;
	npd.printParticle();
      }
      
      abs_sec = absorptionCrossSection(ekin, ptype); 
      abs_sec *= pn_spec * nucleon_densities[0][zone] * nucleon_densities[1][zone] *
	rat * rat1 * vol;
    }

    tot_abs_csec += abs_sec;
    acsecs.push_back(abs_sec);
    qdeutrons.push_back(npd);

    // Neutron-neutron state interacts with pi+, pi0 only
    G4InuclElementaryParticle nnd = generateQuasiDeutron(2, 2, zone);
    
    if (neutronNumberCurrent < 2 || !(ptype == pi0 || ptype == pip)) {
      abs_sec = 0.0;
    } else {
      dummy_convertor.setTarget(nnd.getMomentum(), nnd.getMass());
      
      G4double ekin = dummy_convertor.getKinEnergyInTheTRS();
      
      if (verboseLevel > 2) {
	G4cout << " ptype=" << ptype << " using nn target" << G4endl;
	nnd.printParticle();
      }
      
      abs_sec = absorptionCrossSection(ekin, ptype); 
      abs_sec *= nucleon_densities[1][zone] * nucleon_densities[1][zone] *
	rat1 * rat1 * vol; 
    } 

    tot_abs_csec += abs_sec;
    acsecs.push_back(abs_sec);
    qdeutrons.push_back(nnd);

    // Select quasideuteron interaction from non-zero cross-section choices
    if (verboseLevel > 2){
      G4cout << " rod1 " << acsecs[0] << " rod2 " << acsecs[1]  
	     << " rod3 " << acsecs[2] << G4endl;
    }
    
    if (tot_abs_csec > small) {
      G4double pw = -path * tot_abs_csec;
      
      if (pw < -huge_num) pw = -huge_num;
      pw = 1.0 - std::exp(pw);
      
      if (verboseLevel > 2){
	G4cout << " pw " << pw << G4endl;
      }
      
      G4double apath = path;
      
      if (inuclRndm() < pw) 
	apath = -1.0 / tot_abs_csec * std::log(1.0 - pw * inuclRndm());
      
      if (cparticle.young(young_cut, apath)) apath = path;  
      
      if(verboseLevel > 2){
	G4cout << " apath " << apath << " path " << path << G4endl;
      }
      
      if (apath < path) {	// choose the qdeutron
	G4double sl = inuclRndm() * tot_abs_csec;
	G4double as = 0.0;
	
	for (G4int i = 0; i < 3; i++) {
	  as += acsecs[i];
	  if (sl < as) { 
	    if (verboseLevel > 2) G4cout << " deut type " << i << G4endl; 
	    thePartners.push_back(partner(qdeutrons[i], apath));
	    break;
	  }
	}
      }		// if (apath < path)
    }		// if (tot_abs_csec > small)
  }		// if (cparticle.getParticle().pion())

  
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


const std::vector<G4CascadParticle>&
G4NucleiModel::generateParticleFate(G4CascadParticle& cparticle,
                                    G4ElementaryParticleCollider* theElementaryParticleCollider) {
  if (verboseLevel > 1)
    G4cout << " >>> G4NucleiModel::generateParticleFate" << G4endl;

  if (verboseLevel > 2) {
    G4cout << " cparticle: ";
    cparticle.print();
  }

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

    return outgoing_cparticles;
  }

  G4int npart = thePartners.size();	// Last item is a total-path placeholder

  if (npart == 1) { 		// cparticle is on the next zone entry
    cparticle.propagateAlongThePath(thePartners[0].second);
    cparticle.incrementCurrentPath(thePartners[0].second);
    boundaryTransition(cparticle);
    outgoing_cparticles.push_back(cparticle);
    
    if (verboseLevel > 2) {
      G4cout << " next zone " << G4endl;
      cparticle.print();
    }
  } else {			// there are possible interactions
    if (verboseLevel > 1)
      G4cout << " processing " << npart-1 << " possible interactions" << G4endl;

    G4ThreeVector old_position = cparticle.getPosition();
    G4InuclElementaryParticle& bullet = cparticle.getParticle();
    G4bool no_interaction = true;
    G4int zone = cparticle.getCurrentZone();
    
    for (G4int i=0; i<npart-1; i++) {	// Last item is a total-path placeholder
      if (i > 0) cparticle.updatePosition(old_position); 
      
      G4InuclElementaryParticle& target = thePartners[i].first; 

      if (verboseLevel > 3) {
	if (target.quasi_deutron()) G4cout << " try absorption: ";
	G4cout << " target " << target.type() << " bullet " << bullet.type()
	       << G4endl;
      }

      EPCoutput.reset();
      theElementaryParticleCollider->collide(&bullet, &target, EPCoutput);
      
      if (verboseLevel > 2) {
	EPCoutput.printCollisionOutput();
#ifdef G4CASCADE_CHECK_ECONS
	balance.collide(&bullet, &target, EPCoutput);
	balance.okay();		// Do checks, but ignore result
#endif
      }

      // Don't need to copy list, as "output" isn't changed again below
      const std::vector<G4InuclElementaryParticle>& outgoing_particles = 
	EPCoutput.getOutgoingParticles();
      
      if (!passFermi(outgoing_particles, zone)) continue; // Interaction fails

      // Successful interaction, add results to output list
      cparticle.propagateAlongThePath(thePartners[i].second);
      G4ThreeVector new_position = cparticle.getPosition();

      if (verboseLevel > 2)
	G4cout << " adding " << outgoing_particles.size()
	       << " output particles" << G4endl;
      
      for (G4int ip = 0; ip < G4int(outgoing_particles.size()); ip++) { 
	G4CascadParticle temp(outgoing_particles[ip], new_position, zone, 0.0, 0);
	outgoing_cparticles.push_back(temp);
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
      static G4InuclElementaryParticle prescatCP;	// Avoid memory churn
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
    }
  }	// if (npart == 1) [else]

  return outgoing_cparticles;
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
    G4double pf  = fermi_momenta[type-1][zone];

    if (verboseLevel > 2)
      G4cout << " type " << type << " p " << mom << " pf " << pf << G4endl;
    
    if (mom < pf) {
	if (verboseLevel > 2) G4cout << " rejected by Fermi" << G4endl;
	return false;
    }
  }
  return true; 
}

void G4NucleiModel::boundaryTransition(G4CascadParticle& cparticle) {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::boundaryTransition" << G4endl;
  }

  G4int zone = cparticle.getCurrentZone();

  if (cparticle.movingInsideNuclei() && zone == 0) {
    G4cerr << " boundaryTransition-> in zone 0 " << G4endl;
    return;
  }

  G4LorentzVector mom = cparticle.getMomentum();
  G4ThreeVector pos = cparticle.getPosition();
  
  G4int type = cparticle.getParticle().type();
  
  G4double pr = pos.dot(mom.vect());
  G4double r = pos.mag();
  
  pr /= r;
  
  G4int next_zone = cparticle.movingInsideNuclei() ? zone - 1 : zone + 1;
  
  G4double dv = getPotential(type,zone) - getPotential(type, next_zone);
  //    G4cout << "Potentials for type " << type << " = " 
  //           << getPotential(type,zone) << " , "
  //	   << getPotential(type,next_zone) << G4endl;
  
  G4double qv = dv * dv - 2.0 * dv * mom.e() + pr * pr;
  
  G4double p1r;
  
  if (verboseLevel > 3) {
    G4cout << " type " << type << " zone " << zone << " next " << next_zone
	   << " qv " << qv << " dv " << dv << G4endl;
  }

  if (qv <= 0.0) { 	// reflection
    if (verboseLevel > 3) G4cout << " reflects off boundary" << G4endl;
    p1r = -pr;
    cparticle.incrementReflectionCounter();
  } else {		// transition
    if (verboseLevel > 3) G4cout << " passes thru boundary" << G4endl;
    p1r = std::sqrt(qv);
    if(pr < 0.0) p1r = -p1r;
    cparticle.updateZone(next_zone);
    cparticle.resetReflection();
  }
  
  G4double prr = (p1r - pr) / r;  
  
  if (verboseLevel > 3) {
    G4cout << " prr " << prr << " delta px " << prr*pos.x() << " py "
	   << prr*pos.y()  << " pz " << prr*pos.z() << " mag "
	   << std::fabs(prr*r) << G4endl;
  }

  mom.setVect(mom.vect() + pos*prr);
  cparticle.updateParticleMomentum(mom);
}

G4bool G4NucleiModel::worthToPropagate(const G4CascadParticle& cparticle) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::worthToPropagate" << G4endl;
  }

  const G4double ekin_scale = 2.0;

  G4bool worth = true;

  // NOTE:  Temporarily backing out changes flagged "***" below

  if (cparticle.reflectedNow()) {	// Just reflected -- keep going?
    G4int zone = cparticle.getCurrentZone();
    G4int ip = cparticle.getParticle().type();

    G4double ekin_cut = (cparticle.getParticle().nucleon()) ?
      getFermiKinetic(ip, zone) : 0.;	//*** getPotential(ip, zone);

    worth = cparticle.getParticle().getKineticEnergy()/ekin_scale > ekin_cut;

    if (verboseLevel > 3) {
      G4cout //*** << " type=" << ip
	     << " ekin=" << cparticle.getParticle().getKineticEnergy()
	     << " potential=" << ekin_cut
	     << " : worth? " << worth << G4endl;
    }
  }

  return worth;
}


G4double G4NucleiModel::getRatio(G4int ip) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::getRatio" << G4endl;
  }

  if (ip == 1) {
    if (verboseLevel > 2) {
      G4cout << " current " << protonNumberCurrent << " inp " << protonNumber
	     << G4endl;
    }
    return G4double(protonNumberCurrent)/G4double(protonNumber);
  }

  if (ip == 2) {
    if (verboseLevel > 2){
      G4cout << " current " << neutronNumberCurrent << " inp " << neutronNumber
	     << G4endl;
    }
    return G4double(neutronNumberCurrent)/G4double(neutronNumber);
  }

  return 0.;
}

G4CascadParticle 
G4NucleiModel::initializeCascad(G4InuclElementaryParticle* particle) {
  if (verboseLevel > 1) {
    G4cout << " >>> G4NucleiModel::initializeCascad(particle)" << G4endl;
  }

  const G4double large = 1000.0;

  // FIXME:  Previous version generated random sin(theta), then used -cos(theta)
  //         Using generateWithRandomAngles changes result!
  // G4ThreeVector pos = generateWithRandomAngles(nuclei_radius).vect();
  G4double costh = std::sqrt(1.0 - inuclRndm());
  G4ThreeVector pos = generateWithFixedTheta(-costh, nuclei_radius);

  G4CascadParticle cpart(*particle, pos, number_of_zones, large, 0);

  if (verboseLevel > 2) cpart.print();

  return cpart;
}

void G4NucleiModel::initializeCascad(G4InuclNuclei* bullet, 
				     G4InuclNuclei* target,
				     modelLists& output) {
  if (verboseLevel) {
    G4cout << " >>> G4NucleiModel::initializeCascad(bullet,target,output)"
	   << G4endl;
  }

  const G4double large = 1000.0;
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

      while (casparticles.size() == 0 && itryg < itry_max) {      
	itryg++;
	particles.clear();
      
	//    nucleons coordinates and momenta in nuclei rest frame
	coordinates.clear();
	momentums.clear();
     
	if (ab < 3) { // deuteron, simplest case
	  G4double r = 2.214 - 3.4208 * std::log(1.0 - 0.981 * inuclRndm());
	  G4ThreeVector coord1 = generateWithRandomAngles(r).vect();
	  coordinates.push_back(coord1);
	  coordinates.push_back(-coord1);

	  G4double p = 0.0;
	  G4bool bad = true;
	  G4int itry = 0;

	  while (bad && itry < itry_max) {
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
	    while (badco && itry < itry_max) {
	      if (itry > 0) coordinates.clear();
	      itry++;	
	      G4int i(0);    

	      for (i = 0; i < 2; i++) {
		G4int itry1 = 0;
		G4double s, u, rho; 
		G4double fmax = std::exp(-0.5) / std::sqrt(0.5);

		while (itry1 < itry_max) {
		  itry1++;
		  s = -std::log(inuclRndm());
		  u = fmax * inuclRndm();
		  rho = std::sqrt(s) * std::exp(-s);

		  if (std::sqrt(s) * std::exp(-s) > u && s < s3max) {
		    s = r0forAeq3 * std::sqrt(s);
		    coord1 = generateWithRandomAngles(s).vect();
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
	    b = 1.0 / b;
	    G4double fmax = (1.0 + u * b) * u * std::exp(-u);
	  
	    while (badco && itry < itry_max) {

	      if (itry > 0) coordinates.clear();
	      itry++;
	      G4int i(0);
	    
	      for (i = 0; i < ab-1; i++) {
		G4int itry1 = 0;
		G4double s, u; 

		while (itry1 < itry_max) {
		  itry1++;
		  s = -std::log(inuclRndm());
		  u = fmax * inuclRndm();

		  if (std::sqrt(s) * std::exp(-s) * (1.0 + b * s) > u && s < s4max) {
		    s = r0forAeq4 * std::sqrt(s);
		    coord1 = generateWithRandomAngles(s).vect();
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
	      G4int itry = 0;

	      while(itry < itry_max) {
		itry++;
		u = -std::log(0.879853 - 0.8798502 * inuclRndm());
		x = u * std::exp(-u);

		if(x > inuclRndm()) {
		  p = std::sqrt(0.01953 * u);
		  mom = generateWithRandomAngles(p, massb);
		  momentums.push_back(mom);

		  break;
		}
	      }

	      if(itry == itry_max) {
		G4cout << " can not generate proper momentum for a "
		       << ab << G4endl;

		casparticles.clear();	// Return empty buffer on error
		particles.clear();
		return;
	      } 

	    }
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
    for(ip = 0; ip < G4int(casparticles.size()); ip++) casparticles[ip].print();

    G4cout << " outgoing particles: " << particles.size() << G4endl;
    for(ip = 0; ip < G4int(particles.size()); ip++)
      particles[ip].printParticle();
  }

  return;	// Buffer has been filled
}


G4double G4NucleiModel::absorptionCrossSection(G4double ke, G4int type) const {
  if (type != pionPlus && type != pionMinus && type != pionZero) {
    G4cerr << " absorptionCrossSection only valid for incident pions" << G4endl;
    return 0.;
  }

  // was 0.2 since the beginning, then changed to 1.0 
  // now 0.1 to convert from mb to fm**2
  const G4double corr_fac = 1.0;
  G4double csec = 0.0;
  
  if (ke < 0.3) {
    csec = 0.1106 / std::sqrt(ke) - 0.8 + 0.08 / ((ke-0.123)*(ke-0.123) + 0.0056);

  } else if (ke < 1.0) {
    csec = 3.6735 * (1.0-ke)*(1.0-ke);     
  }

  if (csec < 0.0) csec = 0.0;

  if (verboseLevel > 2) {
    G4cout << " ekin " << ke << " abs. csec " << corr_fac * csec << G4endl;   
  }

  csec *= corr_fac;

  return csec;
}

G4double G4NucleiModel::totalCrossSection(G4double ke, G4int rtype) const
{
  static const G4double keScale[] = {
    0.0,  0.01, 0.013, 0.018, 0.024, 0.032, 0.042, 0.056, 0.075, 0.1,
    0.13, 0.18, 0.24,  0.32,  0.42,  0.56,  0.75,  1.0,   1.3,   1.8,
    2.4,  3.2,  4.2,   5.6,   7.5,  10.0,  13.0,  18.0,  24.0,  32.0};
  static const G4int NBINS = sizeof(keScale)/sizeof(G4double);

  static G4CascadeInterpolator<NBINS> interp(keScale);

  // Pion and nucleon scattering cross-sections are available elsewhere
  switch (rtype) {
  case pro*pro: return G4CascadePPChannel::getCrossSection(ke); break;
  case pro*neu: return G4CascadeNPChannel::getCrossSection(ke); break;
  case pip*pro: return G4CascadePiPlusPChannel::getCrossSection(ke); break;
  case neu*neu: return G4CascadeNNChannel::getCrossSection(ke); break;
  case pim*pro: return G4CascadePiMinusPChannel::getCrossSection(ke); break;
  case pip*neu: return G4CascadePiPlusNChannel::getCrossSection(ke); break;
  case pi0*pro: return G4CascadePiZeroPChannel::getCrossSection(ke); break;
  case pim*neu: return G4CascadePiMinusNChannel::getCrossSection(ke); break;
  case pi0*neu: return G4CascadePiZeroNChannel::getCrossSection(ke); break;
    // Remaining channels are handled locally until arrays are moved
  case kpl*pro:			     
  case k0*neu:  return interp.interpolate(ke, kpPtot); break;
  case kmi*pro:			     
  case k0b*neu: return interp.interpolate(ke, kmPtot); break;
  case kpl*neu:			     
  case k0*pro:  return interp.interpolate(ke, kpNtot); break;
  case kmi*neu:			     
  case k0b*pro: return interp.interpolate(ke, kmNtot); break;
  case lam*pro:			     
  case lam*neu:			     
  case s0*pro:			     
  case s0*neu:  return interp.interpolate(ke, lPtot); break;
  case sp*pro:			     
  case sm*neu:  return interp.interpolate(ke, spPtot); break;
  case sm*pro:			     
  case sp*neu:  return interp.interpolate(ke, smPtot); break;
  case xi0*pro:			     
  case xim*neu: return interp.interpolate(ke, xi0Ptot); break;
  case xim*pro:			     
  case xi0*neu: return interp.interpolate(ke, ximPtot); break;
  default:
    G4cout << " unknown collison type = " << rtype << G4endl; 
  }

  return 0.;		// Failure
}

// Initialize cross-section interpolation tables

const G4double G4NucleiModel::kpPtot[30] = {
   10.0,  10.34, 10.44, 10.61, 10.82, 11.09, 11.43, 11.71, 11.75, 11.8,
   11.98, 12.28, 12.56, 12.48, 12.67, 14.48, 15.92, 17.83, 17.93, 17.88,
   17.46, 17.3,  17.3,  17.4,  17.4,  17.4,  17.4,  17.5,  17.7,  17.8};

const G4double G4NucleiModel::kpNtot[30] = {
    6.64,  6.99,  7.09,  7.27,  7.48,  7.75,  8.1,  8.49,  8.84, 9.31,
    9.8,  10.62, 11.64, 13.08, 14.88, 16.60, 17.5, 18.68, 18.68, 18.29,
   17.81, 17.6,  17.6,  17.6,  17.6,  17.6,  17.7, 17.8,  17.9,  18.0};

const G4double G4NucleiModel::kmPtot[30] = {
 1997.0, 1681.41, 1586.74, 1428.95, 1239.59, 987.12, 671.54, 377.85, 247.30, 75.54,
    71.08, 54.74,   44.08,   44.38,   45.45,  45.07,  41.04,  35.75,  33.22, 30.08,
    27.61, 26.5,    25.2,    24.0,    23.4,   22.8,   22.0,   21.3,   21.0,  20.9};

const G4double G4NucleiModel::kmNtot[30] = {
    6.15,  6.93,  7.16,  7.55,  8.02,  8.65,  9.43, 10.36, 11.34, 12.64,
   14.01, 16.45, 19.32, 23.0,  27.6,  30.92, 29.78, 28.28, 25.62, 23.1,
   22.31, 21.9,  21.73, 21.94, 21.23, 20.5,  20.4,  20.2,  20.1,  20.0};

const G4double G4NucleiModel::lPtot[30] = {
  300.0, 249.07, 233.8, 208.33, 177.78, 137.04, 86.11, 41.41, 28.86, 12.35,
   13.82, 16.76, 20.68,  25.9,   30.37,  31.56, 32.83, 34.5,  34.91, 35.11,
   35.03, 36.06, 35.13,  35.01,  35.0,   35.0,  35.0,  35.0,  35.0,  35.0};

const G4double G4NucleiModel::spPtot[30] = {
  150.0, 146.0, 144.8, 142.8, 140.4, 137.2, 133.2, 127.6, 120.0, 110.0,
   98.06, 84.16, 72.28, 56.58, 43.22, 40.44, 36.14, 30.48, 31.53, 31.92,
   29.25, 28.37, 29.81, 33.15, 33.95, 34.0,  34.0,  34.0,  34.0,  34.0};

const G4double G4NucleiModel::smPtot[30] = {
  937.0, 788.14, 743.48, 669.05, 579.74, 460.65, 311.79, 183.33, 153.65, 114.6,
  105.18, 89.54,  70.58,  45.5,   32.17,  32.54,  32.95,  33.49,  33.55,  33.87,
   34.02, 34.29,  33.93,  33.88,  34.0,   34.0,   34.0,   34.0,   34.0,   34.0};

const G4double G4NucleiModel::xi0Ptot[30] = {
  16.0,  14.72, 14.34, 13.7,  12.93, 11.9,  10.62, 9.29, 8.3,   7.0,
   7.96,  9.56, 11.48, 14.04, 19.22, 25.29, 29.4, 34.8, 34.32, 33.33,
  31.89, 29.55, 27.89, 21.43, 17.0,  16.0,  16.0, 16.0, 16.0,  16.0};

const G4double G4NucleiModel::ximPtot[30] = {
  33.0,  32.5,  32.35, 32.1,  31.8,  31.4,  30.9, 30.2, 29.25, 28.0,
  26.5,  24.6,  22.8,  20.78, 18.22, 19.95, 21.7, 24.0, 24.74, 25.95,
  27.59, 27.54, 23.16, 17.43, 12.94, 12.0,  12.0, 12.0, 12.0,  12.0};
