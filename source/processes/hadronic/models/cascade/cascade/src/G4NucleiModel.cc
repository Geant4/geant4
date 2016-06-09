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
//#define CHC_CHECK

#include "G4NucleiModel.hh"
#include "G4LorentzConvertor.hh"
#include "G4CollisionOutput.hh"


typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4NucleiModel::G4NucleiModel()
  : verboseLevel(2) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::G4NucleiModel" << G4endl;
  }
}

void 
G4NucleiModel::generateModel(G4double a, G4double z) {

  verboseLevel = 2;
  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::generateModel" << G4endl;
  }

  initTotalCrossSections();

  const G4double AU = 1.7234;
  const G4double cuu = 3.3836; 
  const G4double one_third = 1.0 / 3.0;
  const G4double oneBypiTimes4 = 0.0795775; // 1 / 4 Pi
  const G4double pf_coeff = 1.932;
  const G4double pion_vp = 0.007; // in GeV
  const G4double pion_vp_small = 0.007; 
  const G4double radForSmall = 8.0; // fermi
  const G4double piTimes4thirds = 4.189; // 4 Pi/3
  const G4double mproton = 0.93827;
  const G4double mneutron = 0.93957; 
  const G4double alfa3[3] = { 0.7, 0.3, 0.01 }; // listing zone radius
  //  const G4double alfa6[6] = { 0.9, 0.6, 0.4, 0.2, 0.1, 0.05 };

  A = a;
  Z = z;
  neutronNumber = a - z;
  protonNumber = z;
  neutronNumberCurrent = neutronNumber;
  protonNumberCurrent = protonNumber;

// Set binding energies
  G4double dm = bindingEnergy(a, z);

  binding_energies.push_back(0.001 * std::fabs(bindingEnergy(a - 1, z - 1) - dm)); // for P
  binding_energies.push_back(0.001 * std::fabs(bindingEnergy(a - 1, z    ) - dm)); // for N

  G4double CU = cuu * std::pow(a, one_third);
  G4double D1 = CU / AU;
  G4double D = std::exp(-D1);
  G4double CU2 = 0.0; 

  if (a > 3.5) { // a > 3
    std::vector<G4double> ur;

    G4int icase = 0;

    if (a > 11.5) { // a > 11
      // number_of_zones = 6;
      number_of_zones = 3;
      ur.push_back(-D1);

      for (G4int i = 0; i < number_of_zones; i++) {
	// G4double y = std::log((1.0 + D) / alfa6[i] - 1.0);
	G4double y = std::log((1.0 + D)/alfa3[i] - 1.0);
	zone_radii.push_back(CU + AU * y);
	ur.push_back(y);
      };

    } else {
      number_of_zones = 3;
      icase = 1;
      ur.push_back(0.0);
 
      G4double CU1 = CU * CU;
      CU2 = std::sqrt(CU1 * (1.0 - 1.0 / a) + 6.4);

      for (G4int i = 0; i < number_of_zones; i++) {
	G4double y = std::sqrt(-std::log(alfa3[i]));
	zone_radii.push_back(CU2 * y);
	ur.push_back(y);
      };
    }; 

    G4double tot_vol = 0.0;
    std::vector<G4double> v;
    std::vector<G4double> v1;

    G4int i(0);
    for (i = 0; i < number_of_zones; i++) {
      G4double v0;

      if (icase == 0) {
	v0 = volNumInt(ur[i], ur[i + 1], CU, D1);

      } else {
	v0 = volNumInt1(ur[i], ur[i + 1], CU2);
      };
 
      v.push_back(v0);
      tot_vol += v0;
      v0 = (i == 0 ? std::pow(zone_radii[i], G4double(3)) : std::pow(zone_radii[i], G4double(3)) -
	    std::pow(zone_radii[i - 1], G4double(3)));
      v1.push_back(v0);
    }

    // proton
    G4double dd0 = 3.0 * z * oneBypiTimes4 / tot_vol;

    std::vector<G4double> rod;
    std::vector<G4double> pf;
    std::vector<G4double> vz;

    for (i = 0; i < number_of_zones; i++) {
      G4double rd = dd0 * v[i] / v1[i];
      rod.push_back(rd);
      G4double pff = pf_coeff * std::pow(rd, one_third);
      pf.push_back(pff);
      vz.push_back(0.5 * pff * pff / mproton + binding_energies[0]);
    }

    nucleon_densities.push_back(rod);
    zone_potentials.push_back(vz);
    fermi_momenta.push_back(pf);
    //  neutron stuff
    dd0 = 3.0 * (a - z) * oneBypiTimes4 / tot_vol;
    rod.resize(0);
    pf.resize(0);
    vz.resize(0);

    for (i = 0; i < number_of_zones; i++) {
      G4double rd = dd0 * v[i] / v1[i];
      rod.push_back(rd);
      G4double pff = pf_coeff * std::pow(rd, one_third);
      pf.push_back(pff);
      vz.push_back(0.5 * pff * pff / mneutron + binding_energies[1]);
    };

    nucleon_densities.push_back(rod);
    zone_potentials.push_back(vz);
    fermi_momenta.push_back(pf);

    // pion stuff (primitive)
    std::vector<G4double> vp(number_of_zones, pion_vp);
    zone_potentials.push_back(vp);

    // kaon potential (primitive)
    std::vector<G4double> kp(number_of_zones, -0.015);
    zone_potentials.push_back(kp);

    // hyperon potential (primitive)
    std::vector<G4double> hp(number_of_zones, 0.03);
    zone_potentials.push_back(hp);

  } else { // a < 4
    number_of_zones = 1;
    zone_radii.push_back(radForSmall);
    G4double vol = 1.0 / piTimes4thirds / std::pow(zone_radii[0], G4double(3));
    std::vector<G4double> rod;
    std::vector<G4double> pf;
    std::vector<G4double> vz;

    G4int i(0);

    for (i = 0; i < number_of_zones; i++) {
      G4double rd = vol;
      rod.push_back(rd);
      G4double pff = pf_coeff * std::pow(rd, one_third);
      pf.push_back(pff);
      vz.push_back(0.5 * pff * pff / mproton + binding_energies[0]);
    };

    nucleon_densities.push_back(rod);
    zone_potentials.push_back(vz);
    fermi_momenta.push_back(pf);

    // neutron 
    rod.resize(0);
    pf.resize(0);
    vz.resize(0);

    for (i = 0; i < number_of_zones; i++) {
      G4double rd = vol;
      rod.push_back(rd);
      G4double pff = pf_coeff * std::pow(rd, one_third);
      pf.push_back(pff);
      vz.push_back(0.5 * pff * pff / mneutron + binding_energies[1]);
    };

    nucleon_densities.push_back(rod);
    zone_potentials.push_back(vz);
    fermi_momenta.push_back(pf);

    // pion (primitive)
    std::vector<G4double> vp(number_of_zones, pion_vp_small);
    zone_potentials.push_back(vp);
  
    // kaon potential (primitive)
    std::vector<G4double> kp(number_of_zones, -0.015);
    zone_potentials.push_back(kp);

    // hyperon potential (primitive)
    std::vector<G4double> hp(number_of_zones, 0.03);
    zone_potentials.push_back(hp);

  }; 
  nuclei_radius = zone_radii[zone_radii.size() - 1];

}


G4double
G4NucleiModel::volNumInt(G4double r1, G4double r2, 
			 G4double, G4double d1) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::volNumInt" << G4endl;
  }

  const G4double au3 = 5.11864;
  const G4double epsilon = 1.0e-3;
  const G4int itry_max = 1000;
  G4double d2 = 2.0 * d1;
  G4double dr = r2 - r1;
  G4double fi = 0.5 * (r1 * (r1 + d2) / (1.0 + std::exp(r1)) + r2 * (r2 + d2) / (1.0 + std::exp(r2)));
  G4double fun1 = fi * dr;
  G4double fun;
  G4double jc = 1;
  G4double dr1 = dr;
  G4int itry = 0;

  while (itry < itry_max) {
    dr *= 0.5;
    itry++;

    G4double r = r1 - dr;
    fi = 0.0;
    G4int jc1 = G4int(std::pow(2.0, jc - 1) + 0.1);

    for (G4int i = 0; i < jc1; i++) { 
      r += dr1; 
      fi += r * (r + d2) / (1.0 + std::exp(r));
    };

    fun = 0.5 * fun1 + fi * dr;

    if (std::fabs((fun - fun1) / fun) > epsilon) {
      jc++;
      dr1 = dr;
      fun1 = fun;
    } else {
      break;
    }

  }

  if (verboseLevel > 2){
    if(itry == itry_max) G4cout << " volNumInt-> n iter " << itry_max << G4endl;
  }

  return au3 * (fun + d1 * d1 * std::log((1.0 + std::exp(-r1)) / (1.0 + std::exp(-r2))));
}


G4double
G4NucleiModel::volNumInt1(G4double r1, G4double r2, 
			  G4double cu2) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::volNumInt1" << G4endl;
  }

  const G4double epsilon = 1.0e-3;
  const G4int itry_max = 1000;

  G4double dr = r2 - r1;
  G4double fi = 0.5 * (r1 * r1 * std::exp(-r1 * r1) + r2 * r2 * std::exp(-r2 * r2));
  G4double fun1 = fi * dr;
  G4double fun;
  G4double jc = 1;
  G4double dr1 = dr;
  G4int itry = 0;

  while (itry < itry_max) {
    dr *= 0.5;
    itry++;
    G4double r = r1 - dr;
    fi = 0.0;
    G4int jc1 = int(std::pow(2.0, jc - 1) + 0.1);

    for (G4int i = 0; i < jc1; i++) { 
      r += dr1; 
      fi += r * r * std::exp(-r * r);
    }

    fun = 0.5 * fun1 + fi * dr;  

    if (std::fabs((fun - fun1) / fun) > epsilon) {
      jc++;
      dr1 = dr;
      fun1 = fun;

    } else {
      break;
    }

  }

  if (verboseLevel > 2){
    if (itry == itry_max) G4cout << " volNumInt1-> n iter " << itry_max << G4endl;
  }

  return std::pow(cu2, G4double(3)) * fun;
}


void G4NucleiModel::printModel() const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::printModel" << G4endl;
  }

  G4cout << " nuclei model for A " << A << " Z " << Z << G4endl
	 << " proton binding energy " << binding_energies[0] << 
    " neutron binding energy " << binding_energies[1] << G4endl
	 << " Nculei radius " << nuclei_radius << " number of zones " <<
    number_of_zones << G4endl;

  for (G4int i = 0; i < number_of_zones; i++)

    G4cout << " zone " << i+1 << " radius " << zone_radii[i] << G4endl
	   << " protons: density " << getDensity(1,i) << " PF " << 
      getFermiMomentum(1,i) << " VP " << getPotential(1,i) << G4endl
	   << " neutrons: density " << getDensity(2,i) << " PF " << 
      getFermiMomentum(2,i) << " VP " << getPotential(2,i) << G4endl
	   << " pions: VP " << getPotential(3,i) << G4endl;
}


G4InuclElementaryParticle 
G4NucleiModel::generateNucleon(G4int type, G4int zone) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::generateNucleon" << G4endl;
  }

  const G4double one_third = 1.0 / 3.0;

//G4double pmod = getFermiMomentum(type, zone) * std::pow(inuclRndm(), one_third);

  G4double pmod = fermi_momenta[type - 1][zone] * std::pow(inuclRndm(), one_third);

  G4CascadeMomentum mom;
  std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
  G4double FI = randomPHI();
  G4double pt = pmod * COS_SIN.second;

  mom[1] = pt * std::cos(FI);
  mom[2] = pt * std::sin(FI);
  mom[3] = pmod * COS_SIN.first;

  return G4InuclElementaryParticle(mom, type);
}


G4InuclElementaryParticle
G4NucleiModel::generateQuasiDeutron(G4int type1, G4int type2,
				    G4int zone) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::generateQuasiDeutron" << G4endl;
  }

  G4CascadeMomentum mom = generateNucleon(type1, zone).getMomentum(); 
  G4CascadeMomentum mom1 = generateNucleon(type2, zone).getMomentum();
  G4CascadeMomentum dmom;

  for (G4int i = 1; i < 4; i++) dmom[i] = mom[i] + mom1[i]; 

  G4int dtype = 0;

  if (type1 * type2 == 1) {
    dtype = 111;

  } else if (type1 * type2 == 2) { 
    dtype = 112;

  } else if (type1 * type2 == 4) {
    dtype = 122;
  }; 

  return G4InuclElementaryParticle(dmom, dtype);
}


partners 
G4NucleiModel::generateInteractionPartners(G4CascadParticle& cparticle) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::generateInteractionPartners" << G4endl;
  }

  const G4double pi4by3 = 4.1887903; // 4 Pi / 3
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

  partners thePartners;

  G4int ptype = cparticle.getParticle().type();
  G4int zone = cparticle.getCurrentZone();
  G4double pmass = cparticle.getParticle().getMass();
  const G4CascadeMomentum& pmom = cparticle.getParticle().getMomentum();
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
  };  

  G4double path = cparticle.getPathToTheNextZone(r_in, r_out);

  if (verboseLevel > 2){
    G4cout << " r_in " << r_in << " r_out " << r_out << " path " << path << G4endl;
  }

  if (path < -small) { // something wrong
    return thePartners;

  } else if (std::fabs(path) < small) { // just on the boundary
    path = 0.0; 

    G4InuclElementaryParticle particle;

    thePartners.push_back(partner(particle, path));

  } else { // normal case  
    std::vector<G4InuclElementaryParticle> particles;
    G4LorentzConvertor dummy_convertor;

    dummy_convertor.setBullet(pmom, pmass);
  
    for (G4int ip = 1; ip < 3; ip++) { 
      G4InuclElementaryParticle particle = generateNucleon(ip, zone);
      dummy_convertor.setTarget(particle.getMomentum(), particle.getMass());
      G4double ekin = dummy_convertor.getKinEnergyInTheTRS();
      G4double csec = totalCrossSection(ekin, ptype * ip);

      if(verboseLevel > 2){
	G4cout << " ip " << ip << " ekin " << ekin << " csec " << csec << G4endl;
      }

      G4double dens = nucleon_densities[ip - 1][zone];
      G4double rat = getRatio(ip);
      G4double pw = -path * dens * csec * rat;

      if (pw < -huge_num) pw = -huge_num;
      pw = 1.0 - std::exp(pw);

      if (verboseLevel > 2){
	G4cout << " pw " << pw << " rat " << rat << G4endl;
      }

      G4double spath = path;

      if (inuclRndm() < pw) {
	spath = -1.0 / dens / csec / rat * std::log(1.0 - pw * inuclRndm());
	if (cparticle.young(young_cut, spath)) spath = path;

	if (verboseLevel > 2){
	  G4cout << " ip " << ip << " spath " << spath << G4endl;
	}

      };
      if (spath < path) thePartners.push_back(partner(particle, spath));
    };  

    if (verboseLevel > 2){
      G4cout << " after nucleons " << thePartners.size() << " path " << path << G4endl;
    }

    if (cparticle.getParticle().pion()) { // absorption possible

      std::vector<G4InuclElementaryParticle> qdeutrons;
      std::vector<G4double> acsecs;

      G4double tot_abs_csec = 0.0;
      G4double abs_sec;
      G4double vol = std::pow(zone_radii[zone], G4double(3));

      if (zone > 0) vol -= std::pow(zone_radii[zone - 1], G4double(3));
      vol *= pi4by3; 

      G4double rat  = getRatio(1); 
      G4double rat1 = getRatio(2); 

      G4InuclElementaryParticle ppd = generateQuasiDeutron(1, 1, zone);

      if (ptype == 7 || ptype == 5) {
	dummy_convertor.setTarget(ppd.getMomentum(), ppd.getMass());

	G4double ekin = dummy_convertor.getKinEnergyInTheTRS();

	abs_sec = absorptionCrosSection(ekin, ptype);
	abs_sec *= nucleon_densities[0][zone] * nucleon_densities[0][zone]*
	  rat * rat * vol; 

      } else {
	abs_sec = 0.0;
      }; 

      // abs_sec = 0.0;
      tot_abs_csec += abs_sec;
      acsecs.push_back(abs_sec);
      qdeutrons.push_back(ppd);

      G4InuclElementaryParticle npd = generateQuasiDeutron(1, 2, zone);

      dummy_convertor.setTarget(npd.getMomentum(), npd.getMass());

      G4double ekin = dummy_convertor.getKinEnergyInTheTRS();

      abs_sec = absorptionCrosSection(ekin, ptype); 
      abs_sec *= pn_spec * nucleon_densities[0][zone] * nucleon_densities[1][zone] *
	rat * rat1 * vol; 
      tot_abs_csec += abs_sec;
      acsecs.push_back(abs_sec);
      qdeutrons.push_back(npd);

      G4InuclElementaryParticle nnd = generateQuasiDeutron(2, 2, zone);

      if (ptype == 7 || ptype == 3) {
	dummy_convertor.setTarget(nnd.getMomentum(), nnd.getMass());

	G4double ekin = dummy_convertor.getKinEnergyInTheTRS();

	abs_sec = absorptionCrosSection(ekin, ptype); 
	abs_sec *= nucleon_densities[1][zone] * nucleon_densities[1][zone] *
	  rat1 * rat1 * vol; 

      } else {
	abs_sec = 0.0;
      }; 

      // abs_sec = 0.0;
      tot_abs_csec += abs_sec;
      acsecs.push_back(abs_sec);
      qdeutrons.push_back(nnd);

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

	if (apath < path) { // chose the qdeutron

	  G4double sl = inuclRndm() * tot_abs_csec;
	  G4double as = 0.0;

	  for (G4int i = 0; i < 3; i++) {
	    as += acsecs[i];

	    if (sl < as) { 

	      if (verboseLevel > 2){
		G4cout << " deut type " << i << G4endl; 
	      }

	      thePartners.push_back(partner(qdeutrons[i], apath));

	      break;
	    };
	  };
	};    
      };
    };  

    if(verboseLevel > 2){
      G4cout << " after deutrons " << thePartners.size() << G4endl;
    }
  
    if (thePartners.size() > 1) { // sort partners
 
      for (G4int i = 0; i < G4int(thePartners.size()) - 1; i++) {

	for (G4int j = i + 1; j < G4int(thePartners.size()); j++) {

	  if (thePartners[i].second > thePartners[j].second) {

	    G4InuclElementaryParticle particle = thePartners[i].first;
	    G4double pathi = thePartners[i].second;
	    thePartners[i] = partner(thePartners[j].first, thePartners[j].second);
	    thePartners[j] = partner(particle, pathi);
	  };
	};
      };
    };

    G4InuclElementaryParticle particle;

    thePartners.push_back(partner(particle, path));
  }; 
 
  return thePartners;
}


std::vector<G4CascadParticle> 
G4NucleiModel::generateParticleFate(G4CascadParticle& cparticle,
                                    G4ElementaryParticleCollider* theElementaryParticleCollider) {

  if (verboseLevel > 3) G4cout << " >>> G4NucleiModel::generateParticleFate" << G4endl;

  std::vector<G4CascadParticle> outgouing_cparticles;

  partners thePartners = generateInteractionPartners(cparticle);

  if(thePartners.empty()) { // smth. is wrong -> needs special treatment

    G4cout << " generateParticleFate-> can not be here " << G4endl;

  } else {
    G4int npart = thePartners.size();

    if (npart == 1) { // cparticle is on the next zone entry
      // need to go here if particle outside nucleus ?
      //
      cparticle.propagateAlongThePath(thePartners[0].second);
      cparticle.incrementCurrentPath(thePartners[0].second);
      boundaryTransition(cparticle);
      outgouing_cparticles.push_back(cparticle);

      if (verboseLevel > 2){
	G4cout << " next zone " << G4endl;
	cparticle.print();
      }

    } else { // there are possible interactions
 
      std::vector<G4double> old_position = cparticle.getPosition();

      G4InuclElementaryParticle bullet = cparticle.getParticle();

      G4bool no_interaction = true;

      G4int zone = cparticle.getCurrentZone();

      for (G4int i = 0; i < npart - 1; i++) {
	if (i > 0) cparticle.updatePosition(old_position); 

	G4InuclElementaryParticle target = thePartners[i].first; 

	if (verboseLevel > 2){
	  if (target.quasi_deutron()) 
	    G4cout << " try absorption: target " << target.type() << " bullet " <<
	              bullet.type() << G4endl;
	}

	G4CollisionOutput output = theElementaryParticleCollider->collide(&bullet, &target);

	if (verboseLevel > 2) output.printCollisionOutput();

	std::vector<G4InuclElementaryParticle> outgoing_particles = 

	  output.getOutgoingParticles();

        if (passFermi(outgoing_particles, zone)) { // interaction
	  cparticle.propagateAlongThePath(thePartners[i].second);
          std::vector<G4double> new_position = cparticle.getPosition();

	  /*
	  // find jet axis for new particles
          G4double incidentE = cparticle.getParticle().getEnergy();
          G4CascadeMomentum jetAxis;
          for (G4int i = 0; i < G4int(outgoing_particles.size()); i++) {
            for (G4int j = 1; j < 4; j++) jetAxis[j] += (outgoing_particles[i].getMomentum())[j];
          }

          // Find pT wrt jet axis for each secondary
	  */

          for (G4int ip = 0; ip < G4int(outgoing_particles.size()); ip++) { 
            G4CascadParticle temp(outgoing_particles[ip], new_position, zone, 0.0, 0);
	    /*
            G4double pathLength = temp.getPathToTheNextZone(0, nuclei_radius);

	    // Get jet axis
            G4CascadeMomentum pmom = temp.getMomentum();
            G4double secMass = temp.getParticle().getMass();
            G4double dot = 0.0;
            G4double pmod = 0.0;
            G4double jmod = 0.0;
            for (G4int i = 1; i < 4; i++) {
              dot += pmom[i]*jetAxis[i];
              pmod += pmom[i]*pmom[i];
              jmod += jetAxis[i]*jetAxis[i];
	    }

	    //            G4double sinTheta = std::sqrt(1.0 - dot*dot/pmod/jmod);
            G4double pT2 = pmod - dot*dot/jmod;
	    // G4cout << " mass = " << secMass << " Energy = " << incidentE << " pT = " << pT << G4endl;
            G4double formationLength = 1.0*0.1973*incidentE/(pT2 + secMass*secMass);
            if(formationLength > pathLength) {
	      //              G4cout << " formation length = " << formationLength 
	      //                     << " path length = " << pathLength << G4endl;
              temp.propagateAlongThePath(pathLength);
              temp.incrementCurrentPath(pathLength);
              temp.updateZone(number_of_zones-1);
	    }
	    */
            outgouing_cparticles.push_back(temp);
          }

          no_interaction = false;
	  current_nucl1 = 0;
	  current_nucl2 = 0;
#ifdef CHC_CHECK
	  G4double out_charge = 0.0;

	  for (G4int ip = 0; ip < outgoing_particles.size(); ip++) 
	    out_charge += outgoing_particles[ip].getCharge();

	  G4cout << " multiplicity " << outgoing_particles.size() <<
	    " bul type " << bullet.type() << " targ type " << target.type() << 
	    G4endl << " initial charge " << bullet.getCharge() + target.getCharge() 
		 << " out charge " << out_charge << G4endl;  
#endif

	  if (verboseLevel > 2){
	    G4cout << " partner type " << target.type() << G4endl;
	  }

	  if (target.nucleon()) {
	    current_nucl1 = target.type();

	  } else {
	    if (verboseLevel > 2) G4cout << " good absorption " << G4endl;

	    current_nucl1 = (target.type() - 100) / 10;
	    current_nucl2 = target.type() - 100 - 10 * current_nucl1;
          }   
	  
	  if (current_nucl1 == 1) {
	    protonNumberCurrent -= 1.0;

	  } else {
	    neutronNumberCurrent -= 1.0;
	  }; 

	  if (current_nucl2 == 1) {
	    protonNumberCurrent -= 1.0;

	  } else if(current_nucl2 == 2) {
	    neutronNumberCurrent -= 1.0;
	  };
 
	  break;
        }; 
      }  // loop over partners

      if (no_interaction) { // still no interactions
	cparticle.updatePosition(old_position); 
	cparticle.propagateAlongThePath(thePartners[npart - 1].second);
	cparticle.incrementCurrentPath(thePartners[npart - 1].second);
	boundaryTransition(cparticle);
	outgouing_cparticles.push_back(cparticle);
      };
    }; 
  }; 

  return outgouing_cparticles;
}

G4bool G4NucleiModel::passFermi(const std::vector<G4InuclElementaryParticle>& particles, 
				G4int zone) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::passFermi" << G4endl;
  }

  for (G4int i = 0; i < G4int(particles.size()); i++) {

    if (particles[i].nucleon()) {

      if (verboseLevel > 2){
	G4cout << " type " << particles[i].type() << " p " << particles[i].getMomModule()
	       << " pf " << fermi_momenta[particles[i].type() - 1][zone] << G4endl;
      }

      if (particles[i].getMomModule() < fermi_momenta[particles[i].type() - 1][zone]) {

	if (verboseLevel > 2) {
	  G4cout << " rejected by fermi: type " << particles[i].type() << 
	    " p " << particles[i].getMomModule() << G4endl;
	}

	return false;
      };
    };
  };
  return true; 
}

void G4NucleiModel::boundaryTransition(G4CascadParticle& cparticle) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::boundaryTransition" << G4endl;
  }

  G4int zone = cparticle.getCurrentZone();

  if (cparticle.movingInsideNuclei() && zone == 0) {
    G4cout << " boundaryTransition-> in zone 0 " << G4endl;

  } else {
    G4CascadeMomentum mom = cparticle.getMomentum();
    std::vector<G4double> pos = cparticle.getPosition();

    G4int type = cparticle.getParticle().type();

    G4double pr = 0.0;

    G4double r = 0.0;

    G4int i(0);

    for (i = 0; i < 3; i++) {
      pr += pos[i] * mom[i + 1];
      r += pos[i] * pos[i];
    };

    r = std::sqrt(r);
    pr /= r;

    G4int next_zone = cparticle.movingInsideNuclei() ? zone - 1 : zone + 1;

    G4double dv = getPotential(type,zone) - getPotential(type, next_zone);
    //    G4cout << "Potentials for type " << type << " = " 
    //           << getPotential(type,zone) << " , "
    //	   << getPotential(type,next_zone) << G4endl;

    G4double qv = dv * dv - 2.0 * dv * mom[0] + pr * pr;

    G4double p1r;

    if (verboseLevel > 2){
      G4cout << " type " << type << " zone " << zone
             << " next " << next_zone
             << " qv " << qv << " dv " << dv << G4endl;
    }

    if(qv <= 0.0) { // reflection 
      p1r = -pr;
      cparticle.incrementReflectionCounter();

    } else { // transition
      p1r = std::sqrt(qv);
      if(pr < 0.0) p1r = -p1r;
      cparticle.updateZone(next_zone);
      cparticle.resetReflection();
    };
 
    G4double prr = (p1r - pr) / r;  

    for (i = 0; i < 3; i++) mom[i + 1] += pos[i] * prr;

    cparticle.updateParticleMomentum(mom);
  }; 
}

G4bool G4NucleiModel::worthToPropagate(const G4CascadParticle& cparticle) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::worthToPropagate" << G4endl;
  }

  const G4double cut_coeff = 2.0;

  G4bool worth = true;

  if (cparticle.reflectedNow()) {
    G4int zone = cparticle.getCurrentZone();

    G4int ip = cparticle.getParticle().type();

    if (cparticle.getParticle().getKineticEnergy() < cut_coeff *    
       getFermiKinetic(ip, zone)) worth = false; 

  };

  return worth;
}

G4double G4NucleiModel::getRatio(G4int ip) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::getRatio" << G4endl;
  }

  G4double rat;
  //  G4double ratm;

  // Calculate number of protons and neutrons in local region
  //  G4double Athird = std::pow(A, 0.3333);
  //  G4double Nneut = Athird*(A-Z)/A;
  //  G4double Nprot = Athird*Z/A;

  // Reduce number of 
  if (ip == 1) {
    if (verboseLevel > 2){
      G4cout << " current " << protonNumberCurrent << " inp " << protonNumber << G4endl;
    }

    rat = protonNumberCurrent/protonNumber;

    // Calculate ratio modified for local region
    //    G4double deltaP = protonNumber - protonNumberCurrent;
    //    G4cout << " deltaP = " << deltaP << G4endl;
    //    ratm = std::max(0.0, (Nprot - deltaP)/Nprot);

  } else {
    if (verboseLevel > 2){
      G4cout << " current " << neutronNumberCurrent << " inp " << neutronNumber << G4endl;
    }

    rat = neutronNumberCurrent/neutronNumber;

    // Calculate ratio modified for local region
    //    G4double deltaN = neutronNumber - neutronNumberCurrent;
    //   G4cout << " deltaN = " << deltaN << G4endl;
    //    ratm = std::max(0.0, (Nneut - deltaN)/Nneut);
  }

  //  G4cout << " get ratio: ratm =  " << ratm << G4endl;
  return rat;
  //  return ratm;
}

G4CascadParticle G4NucleiModel::initializeCascad(G4InuclElementaryParticle* particle) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::initializeCascad(G4InuclElementaryParticle* particle)" << G4endl;
  }

  const G4double large = 1000.0;

  G4double s1 = std::sqrt(inuclRndm()); 
  G4double phi = randomPHI();
  G4double rz = nuclei_radius * s1;

  std::vector<G4double> pos(3);

  pos[0] = rz * std::cos(phi);
  pos[1] = rz * std::sin(phi);
  pos[2] = -nuclei_radius * std::sqrt(1.0 - s1 * s1);
 
  G4CascadParticle cpart(*particle, pos, number_of_zones, large, 0);

  if (verboseLevel > 2){
    cpart.print();
  }

  return cpart;
}

std::pair<std::vector<G4CascadParticle>, std::vector<G4InuclElementaryParticle> >
G4NucleiModel::initializeCascad(G4InuclNuclei* bullet, 
				G4InuclNuclei* target) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::initializeCascad(G4InuclNuclei* bullet, G4InuclNuclei* target)" << G4endl;
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

  std::vector<G4CascadParticle> casparticles;
  std::vector<G4InuclElementaryParticle> particles;

// first decide whether it will be cascad or compound final nuclei

  G4double ab = bullet->getA();
  G4double zb = bullet->getZ();
  G4double at = target->getA();
  G4double zt = target->getZ();

  if (ab < max_a_for_cascad) {

    G4double benb = 0.001 * bindingEnergy(ab, zb) / ab;
    G4double bent = 0.001 * bindingEnergy(at, zt) / at;
    G4double ben = benb < bent ? bent : benb;

    if (bullet->getKineticEnergy()/ab > ekin_cut*ben) {
      G4int itryg = 0;

      while (casparticles.size() == 0 && itryg < itry_max) {      
	itryg++;

	if(itryg > 0) particles.resize(0);
      
	//    nucleons coordinates and momenta in nuclei rest frame
	std::vector<std::vector<G4double> > coordinates;
	std::vector<G4CascadeMomentum> momentums;
     
	if (ab < 3.0) { // deutron, simplest case
	  G4double r = 2.214 - 3.4208 * std::log(1.0 - 0.981 * inuclRndm());
	  G4double s = 2.0 * inuclRndm() - 1.0;
	  G4double r1 = r * std::sqrt(1.0 - s * s);
	  std::vector<G4double> coord1(3);
	  G4double phi = randomPHI();
	  coord1[0] = r1 * std::cos(phi);
	  coord1[1] = r1 * std::sin(phi);
	  coord1[2] = r * s;   
	  coordinates.push_back(coord1);
	  G4int i(0);

	  for (i = 0; i < 3; i++) coord1[i] *= -1;
	  coordinates.push_back(coord1);
	  G4double p = 0.0;
	  G4bool bad = true;
	  G4int itry = 0;

	  while (bad && itry < itry_max) {
	    itry++;
	    p = 456.0 * inuclRndm();

	    if (p * p / (p * p + 2079.36) / (p * p + 2079.36) > 1.2023e-4 * inuclRndm() &&
	       p * r > 312.0) bad = false;
	  };

	  if (itry == itry_max)
	    if (verboseLevel > 2){ 
	      G4cout << " deutron bullet generation-> itry = " << itry_max << G4endl;	
	    }

	  p = 0.0005 * p;

	  if (verboseLevel > 2){ 
	    G4cout << " p nuc " << p << G4endl;
	  }
	  G4CascadeMomentum mom;
	  std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
	  G4double FI = randomPHI();
	  G4double P1 = p * COS_SIN.second;
	  mom[1] = P1 * std::cos(FI);
	  mom[2] = P1 * std::sin(FI);
	  mom[3] = p * COS_SIN.first;
	  momentums.push_back(mom);

	  for (i = 1; i < 4; i++) mom[i] *= -1;
	  momentums.push_back(mom);

	} else {
	  G4int ia = int(ab + 0.5);

	  std::vector<G4double> coord1(3);

	  G4bool badco = true;

	  G4int itry = 0;
        
	  if (ab < 4.0) { // a == 3
	    while (badco && itry < itry_max) {
	      if (itry > 0) coordinates.resize(0);
	      itry++;	
	      G4int i(0);    

	      for (i = 0; i < 2; i++) {
		G4int itry1 = 0;
		G4double s; 
		G4double u;
		G4double rho;
		G4double fmax = std::exp(-0.5) / std::sqrt(0.5);

		while (itry1 < itry_max) {
		  itry1++;
		  s = -std::log(inuclRndm());
		  u = fmax * inuclRndm();
		  rho = std::sqrt(s) * std::exp(-s);

		  if (std::sqrt(s) * std::exp(-s) > u && s < s3max) {
		    s = r0forAeq3 * std::sqrt(s);
		    std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
		    u = s * COS_SIN.second;  
		    G4double phi = randomPHI();
		    coord1[0] = u * std::cos(phi);
		    coord1[1] = u * std::sin(phi);
		    coord1[2] = s * COS_SIN.first;   
		    coordinates.push_back(coord1);

		    if (verboseLevel > 2){
		      G4cout << " i " << i << " r " << std::sqrt(coord1[0] * coord1[0] +
							    coord1[1] * coord1[1] + 
							    coord1[2] * coord1[2]) << G4endl;
		    }
		    break;
		  };
		};

		if (itry1 == itry_max) { // bad case
		  coord1[0] = coord1[1] = coord1[2] = 10000.;
		  coordinates.push_back(coord1);

		  break;
		};
	      };

	      for (i = 0; i < 3; i++) coord1[i] = - coordinates[0][i] -
				       coordinates[1][i]; 
	      if (verboseLevel > 2) {
		G4cout << " 3  r " << std::sqrt(coord1[0] * coord1[0] +
					   coord1[1] * coord1[1] + 
					   coord1[2] * coord1[2]) << G4endl;
	      }

	      coordinates.push_back(coord1);  	    
	    
	      G4bool large_dist = false;

	      for (i = 0; i < 2; i++) {
		for (G4int j = i+1; j < 3; j++) {
		  G4double r2 = std::pow(coordinates[i][0] - coordinates[j][0], G4double(2)) +
		    std::pow(coordinates[i][1] - coordinates[j][1], G4double(2)) +
		    std::pow(coordinates[i][2] - coordinates[j][2], G4double(2));

		  if (verboseLevel > 2) {
		    G4cout << " i " << i << " j " << j << " r2 " << r2 << G4endl;
		  }

		  if (r2 > r_large2for3) {
		    large_dist = true;

		    break; 
		  };      
		};

		if (large_dist) break;
	      }; 

	      if(!large_dist) badco = false;

	    };

	  } else { // a >= 4
	    G4double b = 3./(ab - 2.0);
	    G4double b1 = 1.0 - b / 2.0;
	    G4double u = b1 + std::sqrt(b1 * b1 + b);
	    b = 1.0 / b;
	    G4double fmax = (1.0 + u * b) * u * std::exp(-u);
	  
	    while (badco && itry < itry_max) {

	      if (itry > 0) coordinates.resize(0);
	      itry++;
	      G4int i(0);
	    
	      for (i = 0; i < ia-1; i++) {
		G4int itry1 = 0;
		G4double s; 
		G4double u;

		while (itry1 < itry_max) {
		  itry1++;
		  s = -std::log(inuclRndm());
		  u = fmax * inuclRndm();

		  if (std::sqrt(s) * std::exp(-s) * (1.0 + b * s) > u && s < s4max) {
		    s = r0forAeq4 * std::sqrt(s);
		    std::pair<double, double> COS_SIN = randomCOS_SIN();
		    u = s * COS_SIN.second;  
		    G4double phi = randomPHI();
		    coord1[0] = u*std::cos(phi);
		    coord1[1] = u*std::sin(phi);
		    coord1[2] = s*COS_SIN.first;   
		    coordinates.push_back(coord1);

		    if (verboseLevel > 2) {
		      G4cout << " i " << i << " r " << std::sqrt(coord1[0]  * coord1[0] +
							    coord1[1] * coord1[1] + 
							    coord1[2] * coord1[2]) << G4endl;
		    }

		    break;
		  };
		};

		if (itry1 == itry_max) { // bad case
		  coord1[0] = coord1[1] = coord1[2] = 10000.0;
		  coordinates.push_back(coord1);

		  break;
		};
	      };

	      for(i = 0; i < 3; i++) {
		coord1[i] = 0.0;

		for(G4int j = 0; j < ia -1; j++) coord1[i] -= coordinates[j][i];
	      };

	      coordinates.push_back(coord1);   

	      if (verboseLevel > 2){
		G4cout << " last r " << std::sqrt(coord1[0] * coord1[0] +
					     coord1[1] * coord1[1] + 
					     coord1[2] * coord1[2]) << G4endl;
	      }
	    
	      G4bool large_dist = false;

	      for (i = 0; i < ia-1; i++) {
		for (G4int j = i+1; j < ia; j++) {
	     
		  G4double r2 = std::pow(coordinates[i][0] - coordinates[j][0], G4double(2)) +
		   
		    std::pow(coordinates[i][1]-coordinates[j][1], G4double(2)) +
		    std::pow(coordinates[i][2] - coordinates[j][2], G4double(2));

		  if (verboseLevel > 2){
		    G4cout << " i " << i << " j " << j << " r2 " << r2 << G4endl;
		  }

		  if (r2 > r_large2for4) {
		    large_dist = true;

		    break; 
		  };      
		};

		if (large_dist) break;
	      }; 

	      if (!large_dist) badco = false;
	    };
	  }; 

	  if(badco) {
	    G4cout << " can not generate the nucleons coordinates for a " << ab <<
	      G4endl;	

	    return std::pair<std::vector<G4CascadParticle>, std::vector<G4InuclElementaryParticle> >
	      (casparticles, particles);

	  } else { // momentums
	    G4double p;
	    G4double u;
	    G4double x;
	    G4CascadeMomentum mom;
	    //G4bool badp = True;
	    G4int i(0);

	    for (i = 0; i < ia - 1; i++) {
	      G4int itry = 0;

	      while(itry < itry_max) {
		itry++;
		u = -std::log(0.879853 - 0.8798502 * inuclRndm());
		x = u * std::exp(-u);

		if(x > inuclRndm()) {
		  p = std::sqrt(0.01953 * u);
		  std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
		  G4double pt = p * COS_SIN.second;  
		  G4double phi = randomPHI();
		  mom[1] = pt * std::cos(phi);
		  mom[2] = pt * std::sin(phi);
		  mom[3] = p * COS_SIN.first;   
		  momentums.push_back(mom);

		  break;
		};
	      };

	      if(itry == itry_max) {
		G4cout << " can not generate proper momentum for a " << ab << G4endl;

		return std::pair<std::vector<G4CascadParticle>, std::vector<G4InuclElementaryParticle> >
		  (casparticles, particles);
	      }; 

	    };
	    // last momentum

	    for(i = 1; i < 4; i++) {
	      mom[i] = 0.;

	      for(G4int j = 0; j < ia -1; j++) mom[i] -= momentums[j][i]; 
	    };

	    momentums.push_back(mom);
	  }; 
	}
 
	// Coordinates and momenta at rest are generated, now back to the lab
	G4double rb = 0.0;
	G4int i(0);

	for(i = 0; i < G4int(coordinates.size()); i++) {      
	  G4double rp = std::sqrt(coordinates[i][0] * coordinates[i][0] +
			     coordinates[i][1] * coordinates[i][1] +
			     coordinates[i][2] * coordinates[i][2]);

	  if(rp > rb) rb = rp;
	};

	// nuclei i.p. as a whole
	G4double s1 = std::sqrt(inuclRndm()); 
	G4double phi = randomPHI();
	G4double rz = (nuclei_radius + rb) * s1;
	std::vector<double> global_pos(3);
	global_pos[0] = rz * std::cos(phi);
	global_pos[1] = rz * std::sin(phi);
	global_pos[2] = -(nuclei_radius + rb) * std::sqrt(1.0 - s1 * s1);

	for (i = 0; i < G4int(coordinates.size()); i++) {
	  coordinates[i][0] += global_pos[0];
	  coordinates[i][1] += global_pos[1];
	  coordinates[i][2] += global_pos[2];
	};  

	// all nucleons at rest
	std::vector<G4InuclElementaryParticle> raw_particles;
	G4int ia = G4int(ab + 0.5);
	G4int iz = G4int(zb + 0.5);

	for (G4int ipa = 0; ipa < ia; ipa++) {
	  G4int knd = ipa < iz ? 1 : 2;
	  raw_particles.push_back(G4InuclElementaryParticle(momentums[ipa], knd));
	}; 
      
	G4InuclElementaryParticle dummy(small_ekin, 1);
	G4LorentzConvertor toTheBulletRestFrame;
	toTheBulletRestFrame.setBullet(dummy.getMomentum(), dummy.getMass());
	toTheBulletRestFrame.setTarget(bullet->getMomentum(),bullet->getMass());
	toTheBulletRestFrame.toTheTargetRestFrame();

	particleIterator ipart;

	for (ipart = raw_particles.begin(); ipart != raw_particles.end(); ipart++) {
	  G4CascadeMomentum mom = 
	    toTheBulletRestFrame.backToTheLab(ipart->getMomentum());
	  ipart->setMomentum(mom); 
	};

	// fill cascad particles and outgoing particles

	for(G4int ip = 0; ip < G4int(raw_particles.size()); ip++) {
	  const G4CascadeMomentum& mom = raw_particles[ip].getMomentum();
	  G4double pmod = std::sqrt(mom[1] * mom[1] + mom[2] * mom[2] + mom[3] * mom[3]);
	  G4double t0 = -(mom[1] * coordinates[ip][0] + mom[2] * coordinates[ip][1] +
			  mom[3] * coordinates[ip][2]) / pmod;
	  G4double det = t0 * t0 + nuclei_radius * nuclei_radius - 
	    coordinates[ip][0] * coordinates[ip][0] - 
	    coordinates[ip][1] * coordinates[ip][1] - 
	    coordinates[ip][2] * coordinates[ip][2];
	  G4double tr = -1.0;

	  if(det > 0.0) {
	    G4double t1 = t0 + std::sqrt(det);
	    G4double t2 = t0 - std::sqrt(det);

	    if(std::fabs(t1) <= std::fabs(t2)) {	 
	      if(t1 > 0.0) {
		if(coordinates[ip][2] + mom[3] * t1 / pmod <= 0.0) tr = t1;
	      };

	      if(tr < 0.0 && t2 > 0.0) {

		if(coordinates[ip][2] + mom[3] * t2 / pmod <= 0.0) tr = t2;
	      };

	    } else {
	      if(t2 > 0.0) {

		if(coordinates[ip][2] + mom[3] * t2 / pmod <= 0.0) tr = t2;
	      };

	      if(tr < 0.0 && t1 > 0.0) {
		if(coordinates[ip][2] + mom[3] * t1 / pmod <= 0.0) tr = t1;
	      };
	    }; 

	  };

	  if(tr >= 0.0) { // cascad particle
	    coordinates[ip][0] += mom[1] * tr / pmod;
	    coordinates[ip][1] += mom[2] * tr / pmod;
	    coordinates[ip][2] += mom[3] * tr / pmod;
	    casparticles.push_back(
				   G4CascadParticle(raw_particles[ip], coordinates[ip], 
						    number_of_zones, large, 0));

	  } else {
	    particles.push_back(raw_particles[ip]); 
	  }; 
	};
      };    

      if(casparticles.size() == 0) {
	particles.resize(0);

	G4cout << " can not generate proper distribution for " << itry_max << " steps " << G4endl;

      };    
    };
  };

  if(verboseLevel > 2){
    G4cout << " cascad particles: " << casparticles.size() << G4endl;
    G4int ip(0);

    for(ip = 0; ip < G4int(casparticles.size()); ip++) casparticles[ip].print();
    G4cout << " outgoing particles: " << particles.size() << G4endl;

    for(ip = 0; ip < G4int(particles.size()); ip++) particles[ip].printParticle();
  }

  return std::pair<std::vector<G4CascadParticle>, std::vector<G4InuclElementaryParticle> >
    (casparticles, particles);
}


G4double G4NucleiModel::totalCrossSection(G4double ke, G4int rtype) const
{
  const G4double keScale[30] = {
    0.0,  0.01, 0.013, 0.018, 0.024, 0.032, 0.042, 0.056, 0.075, 0.1,
    0.13, 0.18, 0.24,  0.32,  0.42,  0.56,  0.75,  1.0,   1.3,   1.8,
    2.4,  3.2,  4.2,   5.6,   7.5,  10.0,  13.0,  18.0,  24.0,  32.0};

  G4int ik = 29;
  G4double sk = 1.0;
  for (G4int i = 1; i < 30; i++) {
    if (ke <= keScale[i]) {
      ik = i;
      sk = (ke - keScale[ik - 1]) / (keScale[ik] - keScale[ik - 1]);
      break;
    }
  }

  G4double csec = 0.0;

  // pp, nn
  if (rtype == 1 || rtype == 4) {
    csec = PPtot[ik - 1] + sk * (PPtot[ik] - PPtot[ik - 1]);

  // np
  } else if (rtype == 2) {
    csec = NPtot[ik - 1] + sk * (NPtot[ik] - NPtot[ik - 1]);

  // pi+p, pi-n  
  } else if (rtype == 3 || rtype == 10) { 
    csec = pipPtot[ik - 1] + sk * (pipPtot[ik] - pipPtot[ik - 1]);

  // pi-p, pi+n 
  } else if (rtype == 5 || rtype == 6) {
    csec = pimPtot[ik - 1] + sk * (pimPtot[ik] - pimPtot[ik - 1]);

  // pi0p, pi0n
  } else if (rtype == 7 || rtype == 14) {
    csec = pizPtot[ik - 1] + sk * (pizPtot[ik] - pizPtot[ik - 1]);

    // k+ p, k0 n 
  } else if (rtype == 11 || rtype == 30) {
    csec = kpPtot[ik - 1] + sk * (kpPtot[ik] - kpPtot[ik - 1]);

  // k- p, k0b n
  } else if (rtype == 13 || rtype == 34) {
    csec = kmPtot[ik - 1] + sk * (kmPtot[ik] - kmPtot[ik - 1]);

  // k+ n, k0 p
  } else if (rtype == 22 || rtype == 15) {
    csec = kpNtot[ik - 1] + sk * (kpNtot[ik] - kpNtot[ik - 1]);

  // k- n, k0b p
  } else if (rtype == 26 || rtype == 17) {
    csec = kmNtot[ik - 1] + sk * (kmNtot[ik] - kmNtot[ik - 1]);

  // L p, L n, S0 p, S0 n
  } else if (rtype == 21 || rtype == 25 || rtype == 42 || rtype == 50) {
    csec = lPtot[ik - 1] + sk * (lPtot[ik] - lPtot[ik - 1]);

  // Sp p, Sm n
  } else if (rtype == 23 || rtype == 54) {
    csec = spPtot[ik - 1] + sk * (spPtot[ik] - spPtot[ik - 1]);

  // Sm p, Sp n
  } else if (rtype == 27 || rtype == 46) {
    csec = smPtot[ik - 1] + sk * (smPtot[ik] - smPtot[ik - 1]);

  // Xi0 p, Xi- n
  } else if (rtype == 29 || rtype == 62) {
    csec = xi0Ptot[ik - 1] + sk * (xi0Ptot[ik] - xi0Ptot[ik - 1]);

  // Xi- p, Xi0 n
  } else if (rtype == 31 || rtype == 58) {
    csec = ximPtot[ik - 1] + sk * (ximPtot[ik] - ximPtot[ik - 1]);

  } else {
    G4cout << " unknown collison type = " << rtype << G4endl; 
  }

  return csec;
}


void G4NucleiModel::initTotalCrossSections()
{
  const G4double PPtotData[30] = {
  17613.0, 302.9, 257.1, 180.6, 128.4,  90.5,  66.1,  49.4,  36.9, 29.6,
     26.0,  23.1,  22.6,  23.0,  27.0,  32.0,  44.0,  47.04, 44.86, 46.03,
     44.09, 41.81, 41.17, 40.65, 40.15, 40.18, 39.26, 38.36, 38.39, 38.41};

  const G4double NPtotData[30] = {
  20357.0, 912.6, 788.6, 582.1, 415.0, 272.0, 198.8, 145.0, 100.4,  71.1,
     58.8,  45.7,  38.9,  34.4,  34.0,  35.0,  37.5,  39.02, 40.29, 40.72,
     42.36, 41.19, 42.04, 41.67, 40.96, 39.48, 39.79, 39.39, 39.36, 39.34};

  const G4double pipPtotData[30] = {
    0.0,   1.2,   2.5,   3.8,   5.0,  7.0,   9.0,  15.0, 30.0,  64.0,
  130.0, 190.0, 130.0,  56.0,  28.0, 17.14, 19.28, 27.4, 40.05, 32.52,
   30.46, 29.0,  27.26, 25.84, 25.5, 24.5,  24.0,  23.5, 23.0,  23.0};

  const G4double pimPtotData[30] = {
    6.13,  6.4,   6.67,  6.94,  7.22,  7.5,  8.3,  12.0,  14.4,  24.0,
   46.0,  72.04, 43.02, 27.19, 27.32, 43.8, 37.08, 51.37, 34.21, 34.79,
   32.08, 31.19, 30.32, 28.5,  27.0,  25.9, 25.5,  25.2,  25.0,  24.8};

  //  const G4double pizPtotData[30] = {
  //    0.0,   3.55,  4.65,  5.9,   7.75, 10.1,  11.8,  18.0,  27.7, 52.5,
  //  102.0, 150.0, 102.64, 51.03, 34.94, 34.52, 32.45, 44.05, 40.2, 34.93,
  //   32.0,  30.0,  28.29, 26.91, 26.25, 25.25, 24.75, 24.35, 24.0, 23.9};

  // New test
  const G4double pizPtotData[30] = {
    6.43,  7.18,  7.54,  8.01,  8.52,  9.13, 10.22, 14.37, 20.96, 34.73,
   61.07, 98.23, 61.97, 32.62, 28.07, 31.37, 35.15, 40.17, 37.27, 33.49,
   31.06, 29.52, 28.29, 26.91, 26.25, 25.25, 24.75, 24.35, 24.0,  23.9};

  const G4double kpPtotData[30] = {
   10.0,  10.34, 10.44, 10.61, 10.82, 11.09, 11.43, 11.71, 11.75, 11.8,
   11.98, 12.28, 12.56, 12.48, 12.67, 14.48, 15.92, 17.83, 17.93, 17.88,
   17.46, 17.3,  17.3,  17.4,  17.4,  17.4,  17.4,  17.5,  17.7,  17.8};

  const G4double kpNtotData[30] = {
    6.64,  6.99,  7.09,  7.27,  7.48,  7.75,  8.1,  8.49,  8.84, 9.31,
    9.8,  10.62, 11.64, 13.08, 14.88, 16.60, 17.5, 18.68, 18.68, 18.29,
   17.81, 17.6,  17.6,  17.6,  17.6,  17.6,  17.7, 17.8,  17.9,  18.0};

  const G4double kmPtotData[30] = {
 1997.0, 1681.41, 1586.74, 1428.95, 1239.59, 987.12, 671.54, 377.85, 247.30, 75.54,
    71.08, 54.74,   44.08,   44.38,   45.45,  45.07,  41.04,  35.75,  33.22, 30.08,
    27.61, 26.5,    25.2,    24.0,    23.4,   22.8,   22.0,   21.3,   21.0,  20.9};

  const G4double kmNtotData[30] = {
    6.15,  6.93,  7.16,  7.55,  8.02,  8.65,  9.43, 10.36, 11.34, 12.64,
   14.01, 16.45, 19.32, 23.0,  27.6,  30.92, 29.78, 28.28, 25.62, 23.1,
   22.31, 21.9,  21.73, 21.94, 21.23, 20.5,  20.4,  20.2,  20.1,  20.0};

  const G4double lPtotData[30] = {
  300.0, 249.07, 233.8, 208.33, 177.78, 137.04, 86.11, 41.41, 28.86, 12.35,
   13.82, 16.76, 20.68,  25.9,   30.37,  31.56, 32.83, 34.5,  34.91, 35.11,
   35.03, 36.06, 35.13,  35.01,  35.0,   35.0,  35.0,  35.0,  35.0,  35.0};

  const G4double spPtotData[30] = {
  150.0, 146.0, 144.8, 142.8, 140.4, 137.2, 133.2, 127.6, 120.0, 110.0,
   98.06, 84.16, 72.28, 56.58, 43.22, 40.44, 36.14, 30.48, 31.53, 31.92,
   29.25, 28.37, 29.81, 33.15, 33.95, 34.0,  34.0,  34.0,  34.0,  34.0};

  const G4double smPtotData[30] = {
  937.0, 788.14, 743.48, 669.05, 579.74, 460.65, 311.79, 183.33, 153.65, 114.6,
  105.18, 89.54,  70.58,  45.5,   32.17,  32.54,  32.95,  33.49,  33.55,  33.87,
   34.02, 34.29,  33.93,  33.88,  34.0,   34.0,   34.0,   34.0,   34.0,   34.0};

  const G4double xi0PtotData[30] = {
  16.0,  14.72, 14.34, 13.7,  12.93, 11.9,  10.62, 9.29, 8.3,   7.0,
   7.96,  9.56, 11.48, 14.04, 19.22, 25.29, 29.4, 34.8, 34.32, 33.33,
  31.89, 29.55, 27.89, 21.43, 17.0,  16.0,  16.0, 16.0, 16.0,  16.0};

  const G4double ximPtotData[30] = {
  33.0,  32.5,  32.35, 32.1,  31.8,  31.4,  30.9, 30.2, 29.25, 28.0,
  26.5,  24.6,  22.8,  20.78, 18.22, 19.95, 21.7, 24.0, 24.74, 25.95,
  27.59, 27.54, 23.16, 17.43, 12.94, 12.0,  12.0, 12.0, 12.0,  12.0};

  for (G4int i = 0; i < 30; i++) {
    PPtot[i] = PPtotData[i];
    NPtot[i] = NPtotData[i];
    pipPtot[i] = pipPtotData[i];
    pimPtot[i] = pimPtotData[i];
    pizPtot[i] = pizPtotData[i];
    kpPtot[i] = kpPtotData[i];
    kpNtot[i] = kpNtotData[i];
    kmPtot[i] = kmPtotData[i];
    kmNtot[i] = kmNtotData[i];
    lPtot[i] = lPtotData[i];
    spPtot[i] = spPtotData[i];
    smPtot[i] = smPtotData[i];
    xi0Ptot[i] = xi0PtotData[i];
    ximPtot[i] = ximPtotData[i];
  }

}
