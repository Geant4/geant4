//#define CHC_CHECK

#include "G4NucleiModel.hh"
#include "G4LorentzConvertor.hh"
#include "G4CollisionOutput.hh"

typedef G4std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4NucleiModel::G4NucleiModel()
  : verboseLevel(2) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::G4NucleiModel" << G4endl;
  }
};

void G4NucleiModel::generateModel(G4double a, 
				  G4double z) {

  verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::generateModel" << G4endl;
  }

  const G4double AU = 1.7234;

  const G4double cuu = 3.3836; 

  const G4double one_third = 1.0 / 3.0;

  const G4double oneBypiTimes4 = 0.0795775; // 1 / 4 Pi

  const G4double pf_coeff = 1.932;

  const G4double pion_vp = 0.007; // ::: make enetgy units explisit

  const G4double pion_vp_small = 0.007; 

  const G4double radForSmall = 8.0; // fermi

  const G4double piTimes4thirds = 4.189; // 4 Pi/3

  const G4double mproton = 0.93827;

  const G4double mneutron = 0.93957; 

  const G4double alfa3[3] = { 0.7, 0.3, 0.01 };

  //  const G4double alfa6[6] = { 0.9, 0.6, 0.4, 0.2, 0.1, 0.05 };

  A = a;
  Z = z;
  neutronNumber = a - z;
  protonNumber = z;
  neutronNumberCurrent = neutronNumber;
  protonNumberCurrent = protonNumber;

// set binding energies

  G4double dm = bindingEnergy(a, z);

  binding_energies.push_back(0.001 * fabs(bindingEnergy(a - 1, z - 1) - dm)); // for P
  binding_energies.push_back(0.001 * fabs(bindingEnergy(a - 1, z) - dm)); // for N

  G4double CU = cuu * pow(a, one_third);

  G4double D1 = CU / AU;

  G4double D = exp(-D1);

  G4double CU2 = 0.0; 

  if(a > 3.5) { // a > 3

    G4std::vector<G4double> ur;

    G4int icase = 0;

    if(a > 11.5) { // a > 11
      //    number_of_zones = 6;
      number_of_zones = 3;
      ur.push_back(-D1);
      for(G4int i = 0; i < number_of_zones; i++) {
	//      G4double y = log((1.0 + D) / alfa6[i] - 1.0);
	G4double y = log((1.0 + D)/alfa3[i] - 1.0);

	zone_radii.push_back(CU + AU * y);
	ur.push_back(y);
      };
    } else {
      number_of_zones = 3;
      icase = 1;
      ur.push_back(0.0);
 
      G4double CU1 = CU * CU;

      CU2 = sqrt(CU1 * (1.0 - 1.0 / a) + 6.4);
      for(G4int i = 0; i < number_of_zones; i++) {

	G4double y = sqrt(-log(alfa3[i]));

	zone_radii.push_back(CU2 * y);
	ur.push_back(y);
      };
    }; 

    G4double tot_vol = 0.0;

    G4std::vector<G4double> v;

    G4std::vector<G4double> v1;

    G4int i(0);
    for(i = 0; i < number_of_zones; i++) {

      G4double v0;

      if(icase == 0) {
	v0 = volNumInt(ur[i], ur[i + 1], CU, D1);
      } else {
	v0 = volNumInt1(ur[i], ur[i + 1], CU2);
      }; 
      v.push_back(v0);
      tot_vol += v0;
      v0 = (i == 0 ? pow(zone_radii[i], 3) : pow(zone_radii[i], 3) -
	    pow(zone_radii[i - 1], 3));
      v1.push_back(v0);
    };
    //  proton stuff

    G4double dd0 = 3.0 * z * oneBypiTimes4 / tot_vol;

    G4std::vector<G4double> rod;

    G4std::vector<G4double> pf;

    G4std::vector<G4double> vz;

    for(i = 0; i < number_of_zones; i++) {

      G4double rd = dd0 * v[i] / v1[i];

      rod.push_back(rd);

      G4double pff = pf_coeff * pow(rd, one_third);

      pf.push_back(pff);
      vz.push_back(0.5 * pff * pff / mproton + binding_energies[0]);
    };
    nucleon_densities.push_back(rod);
    zone_potentials.push_back(vz);
    fermi_momenta.push_back(pf);
    //  neutron stuff
    dd0 = 3.0 * (a - z) * oneBypiTimes4 / tot_vol;
    rod.resize(0);
    pf.resize(0);
    vz.resize(0);
    for(i = 0; i < number_of_zones; i++) {

      G4double rd = dd0 * v[i] / v1[i];

      rod.push_back(rd);

      G4double pff = pf_coeff * pow(rd, one_third);

      pf.push_back(pff);
      vz.push_back(0.5 * pff * pff / mneutron + binding_energies[1]);
    };
    nucleon_densities.push_back(rod);
    zone_potentials.push_back(vz);
    fermi_momenta.push_back(pf);
    //  pion stuff (primitive)
    G4std::vector<G4double> vp(number_of_zones, pion_vp);
    zone_potentials.push_back(vp);
  } else { // a < 4
    number_of_zones = 1;
    zone_radii.push_back(radForSmall);

    G4double vol = 1.0 / piTimes4thirds / pow(zone_radii[0], 3);

    G4std::vector<G4double> rod;

    G4std::vector<G4double> pf;

    G4std::vector<G4double> vz;

    G4int i(0);
    for(i = 0; i < number_of_zones; i++) {

      G4double rd = vol;

      rod.push_back(rd);

      G4double pff = pf_coeff * pow(rd, one_third);

      pf.push_back(pff);
      vz.push_back(0.5 * pff * pff / mproton + binding_energies[0]);
    };
    nucleon_densities.push_back(rod);
    zone_potentials.push_back(vz);
    fermi_momenta.push_back(pf);
    //  neutron stuff
    rod.resize(0);
    pf.resize(0);
    vz.resize(0);
    for(i = 0; i < number_of_zones; i++) {

      G4double rd = vol;

      rod.push_back(rd);

      G4double pff = pf_coeff * pow(rd, one_third);

      pf.push_back(pff);
      vz.push_back(0.5 * pff * pff / mneutron + binding_energies[1]);
    };
    nucleon_densities.push_back(rod);
    zone_potentials.push_back(vz);
    fermi_momenta.push_back(pf);
    //  pion stuff (primitive)

    G4std::vector<G4double> vp(number_of_zones, pion_vp_small);

    zone_potentials.push_back(vp);  
  }; 
  nuclei_radius = zone_radii[zone_radii.size() - 1];
}

G4double G4NucleiModel::volNumInt(G4double r1, 
				  G4double r2, 
				  G4double cu,
				  G4double d1) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::volNumInt" << G4endl;
  }

  const G4double au3 = 5.11864;

  const G4double epsilon = 1.0e-3;

  const G4int itry_max = 1000;

  G4double d2 = 2.0 * d1;

  G4double dr = r2 - r1;

  G4double fi = 0.5 * (r1 * (r1 + d2) / (1.0 + exp(r1)) + r2 * (r2 + d2) / (1.0 + exp(r2)));

  G4double fun1 = fi * dr;

  G4double fun;

  G4double jc = 1;

  G4double dr1 = dr;

  G4int itry = 0;

  while(itry < itry_max) {
    dr *= 0.5;
    itry++;

    G4double r = r1 - dr;

    fi = 0.0;

    G4int jc1 = int(pow(2.0, jc - 1) + 0.1);

    for(G4int i = 0; i < jc1; i++) { 
      r += dr1; 
      fi += r * (r + d2) / (1.0 + exp(r));
    };
    fun = 0.5 * fun1 + fi * dr;
    if(fabs((fun - fun1) / fun) > epsilon) {
      jc++;
      dr1 = dr;
      fun1 = fun;
    } else {

      break;

    }; 
  }; 

  if(verboseLevel > 2){
    if(itry == itry_max) G4cout << " volNumInt-> n iter " << itry_max << G4endl;
  }

  return au3 * (fun + d1 * d1 * log((1.0 + exp(-r1)) / (1.0 + exp(-r2))));
}

G4double G4NucleiModel::volNumInt1(G4double r1, 
				   G4double r2, 
				   G4double cu2) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::volNumInt1" << G4endl;
  }

  const G4double epsilon = 1.0e-3;

  const G4int itry_max = 1000;

  G4double dr = r2 - r1;
  G4double fi = 0.5 * (r1 * r1 * exp(-r1 * r1) + r2 * r2 * exp(-r2 * r2));
  G4double fun1 = fi * dr;
  G4double fun;
  G4double jc = 1;
  G4double dr1 = dr;
  G4int itry = 0;

  while(itry < itry_max) {
    dr *= 0.5;
    itry++;
    G4double r = r1 - dr;
    fi = 0.0;
    G4int jc1 = int(pow(2.0, jc - 1) + 0.1);
    for(G4int i = 0; i < jc1; i++) { 
      r += dr1; 
      fi += r * r * exp(-r * r);
    };
    fun = 0.5 * fun1 + fi * dr;  
    if(fabs((fun - fun1) / fun) > epsilon) {
      jc++;
      dr1 = dr;
      fun1 = fun;
    } else {
      break;
    }; 
  }; 
  if(verboseLevel > 2){
    if(itry == itry_max) G4cout << " volNumInt1-> n iter " << itry_max << G4endl;
  }

  return pow(cu2, 3) * fun;
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

  for(G4int i = 0; i < number_of_zones; i++)

    G4cout << " zone " << i+1 << " radius " << zone_radii[i] << G4endl
	   << " protons: density " << getDensity(1,i) << " PF " << 
      getFermiMomentum(1,i) << " VP " << getPotential(1,i) << G4endl
	   << " neutrons: density " << getDensity(2,i) << " PF " << 
      getFermiMomentum(2,i) << " VP " << getPotential(2,i) << G4endl
	   << " pions: VP " << getPotential(3,i) << G4endl;

}; 

G4InuclElementaryParticle G4NucleiModel::generateNucleon(G4int type, 
							 G4int zone) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::generateNucleon" << G4endl;
  }

  const G4double one_third = 1.0 / 3.0;

//G4double pmod = getFermiMomentum(type, zone) * pow(inuclRndm(), one_third);

  G4double pmod = fermi_momenta[type - 1][zone] * pow(inuclRndm(), one_third);

  G4std::vector<G4double> mom(4);

  G4std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();

  G4double FI = randomPHI();

  G4double pt = pmod * COS_SIN.second;

  mom[1] = pt * cos(FI);
  mom[2] = pt * sin(FI);
  mom[3] = pmod * COS_SIN.first;

  return G4InuclElementaryParticle(mom, type);
}

G4InuclElementaryParticle G4NucleiModel::generateQuasiDeutron(G4int type1, 
							      G4int type2,
							      G4int zone) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::generateQuasiDeutron" << G4endl;
  }

  G4std::vector<G4double> mom = generateNucleon(type1, zone).getMomentum(); 

  G4std::vector<G4double> mom1 = generateNucleon(type2, zone).getMomentum();

  G4std::vector<G4double> dmom(4);

  for(G4int i = 1; i < 4; i++) dmom[i] = mom[i] + mom1[i]; 

  G4int dtype = 0;

  if(type1 * type2 == 1) {
    dtype = 111;
  } else if(type1 * type2 == 2) { 
    dtype = 112;
  } else if(type1 * type2 == 4) {
    dtype = 122;
  }; 

  return G4InuclElementaryParticle(dmom, dtype);
}

partners G4NucleiModel::generateInteractionPartners(G4CascadParticle& cparticle) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::generateInteractionPartners" << G4endl;
  }

  const G4double pi4by3 = 4.1887903; // 4 Pi / 3

  const G4double small = 1.0e-10;

  const G4double huge = 50.0;

  const G4double pn_spec = 1.0;

  //const G4double pn_spec = 0.5;

  //const G4double young_cut = sqrt(10.0) * 0.1;
  //const G4double young_cut = sqrt(10.0) * 0.5;
  //const G4double young_cut = sqrt(10.0) * 0.45;

  const G4double young_cut = sqrt(10.0) * 0.25;

  //const G4double young_cut = sqrt(10.0) * 0.2;
  //const G4double young_cut = 0.0;

  partners thePartners;

  G4int ptype = cparticle.getParticle().type();

  G4int zone = cparticle.getCurrentZone();

  G4double pmass = cparticle.getParticle().getMass();

  G4std::vector<G4double> pmom = cparticle.getParticle().getMomentum();

  G4double r_in;

  G4double r_out;

  if(zone == number_of_zones) { // particle is outside 
    r_in = nuclei_radius;
    r_out = 0.0;
  }
  else if(zone == 0) { // particle is outside core
    r_in = 0.0;
    r_out = zone_radii[0];
  } else {
    r_in = zone_radii[zone - 1];
    r_out = zone_radii[zone];
  };  

  G4double path = cparticle.getPathToTheNextZone(r_in, r_out);

  if(verboseLevel > 2){
    G4cout << " r_in " << r_in << " r_out " << r_out << " path " << path << G4endl;
  }

  if(path < -small) { // something wrong
    return thePartners;
  }
  else if(fabs(path) < small) { // just on the bounday
    path = 0.0; 

    G4InuclElementaryParticle particle;

    thePartners.push_back(partner(particle, path));
  }
  else { // normal case  
  
    G4std::vector<G4InuclElementaryParticle> particles;

    G4LorentzConvertor dummy_convertor;

    dummy_convertor.setBullet(pmom, pmass);
  
    for(G4int ip = 1; ip < 3; ip++) { 

      G4InuclElementaryParticle particle = generateNucleon(ip, zone);

      dummy_convertor.setTarget(particle.getMomentum(), particle.getMass());

      G4double ekin = dummy_convertor.getKinEnergyInTheTRS();

      G4double csec = crossSection(ekin, ptype * ip);

      if(verboseLevel > 2){
	G4cout << " ip " << ip << " ekin " << ekin << " csec " << csec << G4endl;
      }

      G4double dens = nucleon_densities[ip - 1][zone];

      G4double rat = getRatio(ip);

      //    double rat = 1.0;

      G4double pw = -path * dens * csec * rat;

      if(pw < -huge) pw = -huge;
      pw = 1.0 - exp(pw);

      if(verboseLevel > 2){
	G4cout << " pw " << pw << " rat " << rat << G4endl;
      }

      G4double spath = path;

      if(inuclRndm() < pw) {
	spath = -1.0 / dens / csec / rat * log(1.0 - pw * inuclRndm());
	if(cparticle.young(young_cut, spath)) spath = path;

	if(verboseLevel > 2){
	  G4cout << " ip " << ip << " spath " << spath << G4endl;
	}

      };
      if(spath < path) thePartners.push_back(partner(particle, spath));
    };  

    if(verboseLevel > 2){
      G4cout << " after nucleons " << thePartners.size() << " path " << path << G4endl;
    }

    if(cparticle.getParticle().pion()) { // absorption possible

      G4std::vector<G4InuclElementaryParticle> qdeutrons;

      G4std::vector<G4double> acsecs;

      G4double tot_abs_csec = 0.0;

      G4double abs_sec;

      G4double vol = pow(zone_radii[zone], 3);

      if(zone > 0) vol -= pow(zone_radii[zone - 1], 3);
      vol *= pi4by3; 

      G4double rat = getRatio(1); 

      G4double rat1 = getRatio(2); 

      G4InuclElementaryParticle ppd = generateQuasiDeutron(1, 1, zone);
      if(ptype == 7 || ptype == 5) {
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

      if(ptype == 7 || ptype == 3) {
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

      if(verboseLevel > 2){
	G4cout << " rod1 " << acsecs[0] << " rod2 " << acsecs[1]  
	       << " rod3 " << acsecs[2] << G4endl;
      }

      if(tot_abs_csec > small) {
     
	G4double pw = -path * tot_abs_csec;

	if(pw < -huge) pw = -huge;
	pw = 1.0 - exp(pw);

	if(verboseLevel > 2){
	  G4cout << " pw " << pw << G4endl;
	}

	G4double apath = path;

	if(inuclRndm() < pw) 
	  apath = -1.0 / tot_abs_csec * log(1.0 - pw * inuclRndm());
	if(cparticle.young(young_cut, apath)) apath = path;  

	if(verboseLevel > 2){
	  G4cout << " apath " << apath << " path " << path << G4endl;
	}

	if(apath < path) { // chose the qdeutron

	  G4double sl = inuclRndm() * tot_abs_csec;

	  G4double as = 0.0;

	  for(G4int i = 0; i < 3; i++) {
	    as += acsecs[i];
	    if(sl < as) { 
	      if(verboseLevel > 2){
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
  
    if(thePartners.size() > 1) { // sort partners
      for(G4int i = 0; i < G4int(thePartners.size()) - 1; i++) {
	for(G4int j = i + 1; j < G4int(thePartners.size()); j++) {
	  if(thePartners[i].second > thePartners[j].second) {

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

G4std::vector<G4CascadParticle> G4NucleiModel::generateParticleFate(G4CascadParticle& cparticle,
							     G4ElementaryParticleCollider* theElementaryParticleCollider) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::generateParticleFate" << G4endl;
  }

  G4std::vector<G4CascadParticle> outgouing_cparticles;

  partners thePartners = generateInteractionPartners(cparticle);

  if(thePartners.empty()) { // smth. is wrong -> needs special treatment

    G4cout << " generateParticleFate-> can not be here " << G4endl;

  } else {

    G4int npart = thePartners.size();

    if(npart == 1) { // cparticle is on the next zone entry
      cparticle.propagateAlongThePath(thePartners[0].second);
      cparticle.incrementCurrentPath(thePartners[0].second);
      boundaryTransition(cparticle);
      outgouing_cparticles.push_back(cparticle);

      if(verboseLevel > 2){
	G4cout << " next zone " << G4endl;
	cparticle.print();
      }

    } else { // there are possible interactions
 
      G4std::vector<G4double> old_position = cparticle.getPosition();

      G4InuclElementaryParticle bullet = cparticle.getParticle();

      G4bool no_interaction = true;

      G4int zone = cparticle.getCurrentZone();

      for(G4int i = 0; i < npart - 1; i++) {
	if(i > 0) cparticle.updatePosition(old_position); 

	G4InuclElementaryParticle target = thePartners[i].first; 

	if(verboseLevel > 2){
	  if(target.quasi_deutron()) 

	    G4cout << " try absorption: target " << target.type() << " bullet " <<

	      bullet.type() << G4endl;
	}

	G4CollisionOutput output = theElementaryParticleCollider->collide(&bullet, &target);

	if(verboseLevel > 2){
	  output.printCollisionOutput();
	}

	G4std::vector<G4InuclElementaryParticle> outgoing_particles = 

	  output.getOutgoingParticles();
	if(passFermi(outgoing_particles, zone)) { // interaction
	  cparticle.propagateAlongThePath(thePartners[i].second);

	  G4std::vector<G4double> new_position = cparticle.getPosition();

	  for(G4int ip = 0; ip < G4int(outgoing_particles.size()); ip++) 
	    outgouing_cparticles.push_back(G4CascadParticle(outgoing_particles[ip],
							    new_position, zone, 0.0));
	  no_interaction = false;
	  current_nucl1 = 0;
	  current_nucl2 = 0;
#ifdef CHC_CHECK
	  G4double out_charge = 0.0;

	  for(G4int ip = 0; ip < outgoing_particles.size(); ip++) 
	    out_charge += outgoing_particles[ip].getCharge();

	  G4cout << " multiplicity " << outgoing_particles.size() <<
	    " bul type " << bullet.type() << " targ type " << target.type() << 
	    G4endl << " initial charge " << bullet.getCharge() + target.getCharge() 
		 << " out charge " << out_charge << G4endl;  
#endif

	  if(verboseLevel > 2){
	    G4cout << " partner type " << target.type() << G4endl;
	  }

	  if(target.nucleon()) {
	    current_nucl1 = target.type();
	  } else {

	    if(verboseLevel > 2){
	      G4cout << " good absorption " << G4endl;
	    }

	    current_nucl1 = (target.type() - 100) / 10;
	    current_nucl2 = target.type() - 100 - 10 * current_nucl1;
	  };   
	  
	  if(current_nucl1 == 1) {
	    protonNumberCurrent -= 1.0;
	  } else {
	    neutronNumberCurrent -= 1.0;
	  }; 
	  if(current_nucl2 == 1) {
	    protonNumberCurrent -= 1.0;
	  } else if(current_nucl2 == 2) {
	    neutronNumberCurrent -= 1.0;
	  };
 
	  break;
	}; 
      };
      if(no_interaction) { // still now interactions
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

G4bool G4NucleiModel::passFermi(const G4std::vector<G4InuclElementaryParticle>& particles, 
				G4int zone) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::passFermi" << G4endl;
  }

  for(G4int i = 0; i < G4int(particles.size()); i++) {
    if(particles[i].nucleon()) {

      if(verboseLevel > 2){
	G4cout << " type " << particles[i].type() << " p " << particles[i].getMomModule()
	       << " pf " << fermi_momenta[particles[i].type() - 1][zone] << G4endl;
      }

      if(particles[i].getMomModule() < fermi_momenta[particles[i].type() - 1][zone]) {

	if(verboseLevel > 2) {
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

  if(cparticle.movingInsideNuclei() && zone == 0) {
    G4cout << " boundaryTransition-> in zone 0 " << G4endl;
  } else {
   
    G4std::vector<G4double> mom = cparticle.getMomentum();

    G4std::vector<G4double> pos = cparticle.getPosition();

    G4int type = cparticle.getParticle().type();

    G4double pr = 0.0;

    G4double r = 0.0;

    G4int i(0);
    for(i = 0; i < 3; i++) {
      pr += pos[i] * mom[i + 1];
      r += pos[i] * pos[i];
    };
    r = sqrt(r);
    pr /= r;

    G4int next_zone = cparticle.movingInsideNuclei() ? zone - 1 : zone + 1;

    G4double dv = getPotential(type,zone) - getPotential(type, next_zone);

    G4double qv = dv * dv - 2.0 * dv * mom[0] + pr * pr;

    G4double p1r;

    if(verboseLevel > 2){
      cout << " type " << type << " zone " << zone << " next " << next_zone <<
	" qv " << qv << " dv " << dv << endl;
    }

    if(qv <= 0.0) { // reflection 
      p1r = -pr;
      cparticle.incrementReflectionCounter();
    } else { // transition
      p1r = sqrt(qv);
      if(pr < 0.0) p1r = -p1r;
      cparticle.updateZone(next_zone);
      cparticle.resetReflection();
    };
 
    G4double prr = (p1r - pr) / r;  

    for(i = 0; i < 3; i++) mom[i + 1] += pos[i] * prr;
    cparticle.updateParticleMomentum(mom);
  }; 
}

G4bool G4NucleiModel::worthToPropagate(const G4CascadParticle& cparticle) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::worthToPropagate" << G4endl;
  }

  const G4double cut_coeff = 2.0;

  G4bool worth = true;

  if(cparticle.reflectedNow()) {

    G4int zone = cparticle.getCurrentZone();

    G4int ip = cparticle.getParticle().type();

    if(cparticle.getParticle().getKineticEnergy() < cut_coeff *    
       getFermiKinetic(ip, zone)) worth = false; 
  };

  return worth;
}

G4double G4NucleiModel::getRatio(G4int ip) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::getRatio" << G4endl;
  }

  G4double rat;

  if(ip == 1) {

    if(verboseLevel > 2){
      G4cout << " current " << protonNumberCurrent << " inp " << protonNumber << G4endl;
    }

    rat = protonNumberCurrent / protonNumber;
  } else {

    if(verboseLevel > 2){
      G4cout << " current " << neutronNumberCurrent << " inp " << neutronNumber << G4endl;
    }
    rat = neutronNumberCurrent / neutronNumber;
  }; 

  return rat;
}

G4CascadParticle G4NucleiModel::initializeCascad(G4InuclElementaryParticle* particle) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4NucleiModel::initializeCascad(G4InuclElementaryParticle* particle)" << G4endl;
  }

  const G4double large = 1000.0;

  G4double s1 = sqrt(inuclRndm()); 

  G4double phi = randomPHI();

  G4double rz = nuclei_radius * s1;

  G4std::vector<G4double> pos(3);

  pos[0] = rz * cos(phi);
  pos[1] = rz * sin(phi);
  pos[2] = -nuclei_radius * sqrt(1.0 - s1 * s1);
 
  G4CascadParticle cpart(*particle, pos, number_of_zones, large);

  if(verboseLevel > 2){
    cpart.print();
  }

  return cpart;
}

G4std::pair<G4std::vector<G4CascadParticle>, G4std::vector<G4InuclElementaryParticle> >
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

  G4std::vector<G4CascadParticle> casparticles;
  G4std::vector<G4InuclElementaryParticle> particles;

// first decide whether it will be cascad or compound final nuclei

  G4double ab = bullet->getA();

  G4double zb = bullet->getZ();

  G4double at = target->getA();

  G4double zt = target->getZ();

  if(ab < max_a_for_cascad) {

    G4double benb = 0.001 * bindingEnergy(ab, zb) / ab;

    G4double bent = 0.001 * bindingEnergy(at, zt) / at;

    G4double ben = benb < bent ? bent : benb;

    if(bullet->getKineticEnergy()/ab > ekin_cut*ben) {

      int itryg = 0;
      while(casparticles.size() == 0 && itryg < itry_max) {
      
	itryg++;
	if(itryg > 0) particles.resize(0);
      
	//    nucleons coordinates and momenta in nuclei rest frame

	G4std::vector<G4std::vector<G4double> > coordinates;

	G4std::vector<G4std::vector<G4double> > momentums;
     
	if(ab < 3.0) { // deutron, simplest case

	  G4double r = 2.214 - 3.4208 * log(1.0 - 0.981 * inuclRndm());

	  G4double s = 2.0 * inuclRndm() - 1.0;

	  G4double r1 = r * sqrt(1.0 - s * s);

	  G4std::vector<G4double> coord1(3);

	  G4double phi = randomPHI();

	  coord1[0] = r1 * cos(phi);
	  coord1[1] = r1 * sin(phi);
	  coord1[2] = r * s;   
	  coordinates.push_back(coord1);
	  G4int i(0);
	  for(i = 0; i < 3; i++) coord1[i] *= -1;
	  coordinates.push_back(coord1);
        
	  G4double p = 0.0;

	  G4bool bad = true;

	  G4int itry = 0;

	  while(bad && itry < itry_max) {
	    itry++;
	    p = 456.0 * inuclRndm();
	    if(p * p / (p * p + 2079.36) / (p * p + 2079.36) > 1.2023e-4 * inuclRndm() &&
	       p * r > 312.0) bad = false;
	  };
	  if(itry == itry_max)
	    
	    if(verboseLevel > 2){ 
	      cout << " deutron bullet generation-> itry = " << itry_max << endl;	
	    }
	  p = 0.0005 * p;

	  if(verboseLevel > 2){ 
	    cout << " p nuc " << p << endl;
	  }

	  G4std::vector<G4double> mom(4);

	  G4std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();

	  G4double FI = randomPHI();

	  G4double P1 = p * COS_SIN.second;

	  mom[1] = P1 * cos(FI);
	  mom[2] = P1 * sin(FI);
	  mom[3] = p * COS_SIN.first;
	  momentums.push_back(mom);
	  for(i = 1; i < 4; i++) mom[i] *= -1;
	  momentums.push_back(mom);
	} else {

	  G4int ia = int(ab + 0.5);

	  G4std::vector<G4double> coord1(3);

	  G4bool badco = true;

	  G4int itry = 0;
        
	  if(ab < 4.0) { // a == 3
	    while(badco && itry < itry_max) {
	      if(itry > 0) coordinates.resize(0);
	      itry++;	
	      G4int i(0);    
	      for(i = 0; i < 2; i++) {

		G4int itry1 = 0;

		G4double s; 

		G4double u;

		G4double rho;

		G4double fmax = exp(-0.5) / sqrt(0.5);

		while(itry1 < itry_max) {
		  itry1++;
		  s = -log(inuclRndm());
		  u = fmax * inuclRndm();
		  rho = sqrt(s) * exp(-s);
		  if(sqrt(s) * exp(-s) > u && s < s3max) {
		    s = r0forAeq3 * sqrt(s);

		    G4std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();

		    u = s * COS_SIN.second;  

		    G4double phi = randomPHI();

		    coord1[0] = u * cos(phi);
		    coord1[1] = u * sin(phi);
		    coord1[2] = s * COS_SIN.first;   
		    coordinates.push_back(coord1);

		    if (verboseLevel > 2){
		      G4cout << " i " << i << " r " << sqrt(coord1[0] * coord1[0] +
							    coord1[1] * coord1[1] + 
							    coord1[2] * coord1[2]) << G4endl;
		    }

		    break;
		  };
		};
		if(itry1 == itry_max) { // bad case
		  coord1[0] = coord1[1] = coord1[2] = 10000.;
		  coordinates.push_back(coord1);

		  break;
		};
	      };
	      for(i = 0; i < 3; i++) coord1[i] = - coordinates[0][i] -
				       coordinates[1][i]; 

	      if (verboseLevel > 2) {
		G4cout << " 3  r " << sqrt(coord1[0] * coord1[0] +
					   coord1[1] * coord1[1] + 
					   coord1[2] * coord1[2]) << G4endl;
	      }

	      coordinates.push_back(coord1);  	    
	    
	      G4bool large_dist = false;

	      for(i = 0; i < 2; i++) {
		for(G4int j = i+1; j < 3; j++) {

		  G4double r2 = pow(coordinates[i][0] - coordinates[j][0], 2) +
		    pow(coordinates[i][1] - coordinates[j][1], 2) +
		    pow(coordinates[i][2] - coordinates[j][2], 2);

		  if (verboseLevel > 2) {
		    G4cout << " i " << i << " j " << j << " r2 " << r2 << G4endl;
		  }

		  if(r2 > r_large2for3) {
		    large_dist = true;

		    break; 
		  };      
		};
		if(large_dist) break;
	      }; 
	      if(!large_dist) badco = false;
	    };
	  } else { // a >= 4
	
	    G4double b = 3./(ab - 2.0);

	    G4double b1 = 1.0 - b / 2.0;

	    G4double u = b1 + sqrt(b1 * b1 + b);

	    b = 1.0 / b;

	    G4double fmax = (1.0 + u * b) * u * exp(-u);
	  
	    while(badco && itry < itry_max) {
	      if(itry > 0) coordinates.resize(0);
	      itry++;
	      G4int i(0);	    
	      for(i = 0; i < ia-1; i++) {

		G4int itry1 = 0;

		G4double s; 

		G4double u;

		while(itry1 < itry_max) {
		  itry1++;
		  s = -log(inuclRndm());
		  u = fmax * inuclRndm();
		  if(sqrt(s) * exp(-s) * (1.0 + b * s) > u && s < s4max) {
		    s = r0forAeq4 * sqrt(s);

		    G4std::pair<double, double> COS_SIN = randomCOS_SIN();

		    u = s * COS_SIN.second;  
		 
		    G4double phi = randomPHI();

		    coord1[0] = u*cos(phi);
		    coord1[1] = u*sin(phi);
		    coord1[2] = s*COS_SIN.first;   
		    coordinates.push_back(coord1);

		    if (verboseLevel > 2) {
		      G4cout << " i " << i << " r " << sqrt(coord1[0]  * coord1[0] +
							    coord1[1] * coord1[1] + 
							    coord1[2] * coord1[2]) << G4endl;
		    }

		    break;
		  };
		};
		if(itry1 == itry_max) { // bad case
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
		G4cout << " last r " << sqrt(coord1[0] * coord1[0] +
					     coord1[1] * coord1[1] + 
					     coord1[2] * coord1[2]) << G4endl;
	      }
	    
	      G4bool large_dist = false;

	      for(i = 0; i < ia-1; i++) {
		for(G4int j = i+1; j < ia; j++) {
	     
		  G4double r2 = pow(coordinates[i][0] - coordinates[j][0], 2) +
		   
		    pow(coordinates[i][1]-coordinates[j][1], 2) +
		    pow(coordinates[i][2] - coordinates[j][2], 2);

		  if (verboseLevel > 2){
		    G4cout << " i " << i << " j " << j << " r2 " << r2 << G4endl;
		  }

		  if(r2 > r_large2for4) {
		    large_dist = true;

		    break; 
		  };      
		};
		if(large_dist) break;
	      }; 
	      if(!large_dist) badco = false;
	    };
	  }; 
	  if(badco) {

	    G4cout << " can not generate the nucleons coordinates for a " << ab <<
	      G4endl;	

	    return G4std::pair<G4std::vector<G4CascadParticle>, G4std::vector<G4InuclElementaryParticle> >
	      (casparticles, particles);

	  } else { // momentums

	    G4double p;

	    G4double u;

	    G4double x;

	    G4std::vector<G4double> mom(4);

	    //	    G4bool badp = True;

	    G4int i(0);
	    for(i = 0; i < ia - 1; i++) {

	      G4int itry = 0;

	      while(itry < itry_max) {
		itry++;
		u = -log(0.879853 - 0.8798502 * inuclRndm());
		x = u * exp(-u);
		if(x > inuclRndm()) {
		  p = sqrt(0.01953 * u);

		  G4std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();

		  G4double pt = p * COS_SIN.second;  

		  G4double phi = randomPHI();

		  mom[1] = pt * cos(phi);
		  mom[2] = pt * sin(phi);
		  mom[3] = p * COS_SIN.first;   
		  momentums.push_back(mom);

		  break;
		};
	      };
	      if(itry == itry_max) {

		G4cout << " can not generate proper momentum for a " << ab << G4endl;

		return G4std::pair<G4std::vector<G4CascadParticle>, G4std::vector<G4InuclElementaryParticle> >
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
	}; 
	// coordinates and momentums at rest are generated, now back to the lab;


	G4double rb = 0.0;

	G4int i(0);
	for(i = 0; i < G4int(coordinates.size()); i++) {
      
	  G4double rp = sqrt(coordinates[i][0] * coordinates[i][0] +
			     coordinates[i][1] * coordinates[i][1] +
			     coordinates[i][2] * coordinates[i][2]);
	  if(rp > rb) rb = rp;
	};
	//    nuclei i.p. as a whole

	G4double s1 = sqrt(inuclRndm()); 

	G4double phi = randomPHI();
    
	G4double rz = (nuclei_radius + rb) * s1;

	G4std::vector<double> global_pos(3);

	global_pos[0] = rz * cos(phi);
	global_pos[1] = rz * sin(phi);
	global_pos[2] = -(nuclei_radius + rb) * sqrt(1.0 - s1 * s1);
	cout <<"+++++> 2" << endl;
	for(i = 0; i < G4int(coordinates.size()); i++) {
	  coordinates[i][0] += global_pos[0];
	  coordinates[i][1] += global_pos[1];
	  coordinates[i][2] += global_pos[2];
	};  
	//   all nucleons at rest

	G4std::vector<G4InuclElementaryParticle> raw_particles;

	G4int ia = int(ab + 0.5);

	G4int iz = int(zb + 0.5);

	for(G4int ipa = 0; ipa < ia; ipa++) {

	  G4int knd = ipa < iz ? 1 : 2;

	  raw_particles.push_back(G4InuclElementaryParticle(momentums[ipa], knd));
	}; 
      
	G4InuclElementaryParticle dummy(small_ekin, 1);
  
	G4LorentzConvertor toTheBulletRestFrame;
	cout <<"+++++> 3" << endl;
	toTheBulletRestFrame.setBullet(dummy.getMomentum(), dummy.getMass());
	toTheBulletRestFrame.setTarget(bullet->getMomentum(),bullet->getMass());
	toTheBulletRestFrame.toTheTargetRestFrame();

	particleIterator ipart;

	for(ipart = raw_particles.begin(); ipart != raw_particles.end(); ipart++) {

	  G4std::vector<G4double> mom = 
	    toTheBulletRestFrame.backToTheLab(ipart->getMomentum());

	  ipart->setMomentum(mom); 
	};
	//  fill cascad particles and outgoing particles
	for(G4int ip = 0; ip < G4int(raw_particles.size()); ip++) {

	  G4std::vector<G4double> mom = raw_particles[ip].getMomentum();

	  G4double pmod = sqrt(mom[1] * mom[1] + mom[2] * mom[2] + mom[3] * mom[3]);

	  G4double t0 = -(mom[1] * coordinates[ip][0] + mom[2] * coordinates[ip][1] +
			  mom[3] * coordinates[ip][2]) / pmod;

	  G4double det = t0 * t0 + nuclei_radius * nuclei_radius - 
	    coordinates[ip][0] * coordinates[ip][0] - 
	    coordinates[ip][1] * coordinates[ip][1] - 
	    coordinates[ip][2] * coordinates[ip][2];
	 
	  G4double tr = -1.0;

	  if(det > 0.0) {

	    G4double t1 = t0 + sqrt(det);

	    G4double t2 = t0 - sqrt(det);

	    if(fabs(t1) <= fabs(t2)) {	  
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
						    number_of_zones, large));
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
    for(ip = 0; ip < G4int(casparticles.size()); ip++)
      casparticles[ip].print();
    G4cout << " outgoing particles: " << particles.size() << G4endl;
    for(ip = 0; ip < G4int(particles.size()); ip++)
      particles[ip].printParticle();
  }

  return G4std::pair<G4std::vector<G4CascadParticle>, G4std::vector<G4InuclElementaryParticle> >
    (casparticles, particles);
}
