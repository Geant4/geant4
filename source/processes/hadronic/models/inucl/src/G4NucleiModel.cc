
//#define DEBUG
//#define CHC_CHECK

#include "G4NucleiModel.hh"
#include "G4LorentzConvertor.hh"
#include "G4CollisionOutput.hh"

typedef vector<G4InuclElementaryParticle>::iterator particleIterator;

void G4NucleiModel::generateModel(G4double a, G4double z) {

const G4double AU = 1.7234;
const G4double cuu = 3.3836; 
const G4double one_third = 1.0/3.0;
const G4double oneBypiTimes4 = 0.0795775; // 1 / 4 Pi
const G4double pf_coeff = 1.932;
const G4double pion_vp = 0.007 * GeV; 
const G4double pion_vp_small = 0.007 * GeV; 
const G4double radForSmall = 8.0 * fermi;
const G4double piTimes4thirds = 4.189; // 4 Pi/3
const G4double mproton = 0.93827 * GeV;
const G4double mneutron = 0.93957 * GeV; 

const G4double alfa3[3] = { 0.7, 0.3, 0.01 };
const G4double alfa6[6] = { 0.9, 0.6, 0.4, 0.2, 0.1, 0.05 };

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
G4double CU2; 

if(a > 3.5) { // a > 3
  vector<G4double> ur;
  int icase = 0;
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
  }
   else {
    number_of_zones = 3;
    icase = 1;
    ur.push_back(0.0);
    double CU1 = CU * CU;
    CU2 = sqrt(CU1 * (1.0 - 1.0 / a) + 6.4);
    for(int i = 0; i < number_of_zones; i++) {
      G4double y = sqrt(-log(alfa3[i]));
      zone_radii.push_back(CU2 * y);
      ur.push_back(y);
    };
  }; 
  G4double tot_vol = 0.;
  vector<G4double> v;
  vector<G4double> v1;
  for(G4int i = 0; i < number_of_zones; i++) {
    G4double v0;
    if(icase == 0) {
      v0 = volNumInt(ur[i], ur[i + 1], CU, D1);
    }
     else {
      v0 = volNumInt1(ur[i], ur[i + 1], CU2);
    }; 
    v.push_back(v0);
    tot_vol += v0;
    v0 = (i == 0 ? pow(zone_radii[i], 3) : pow(zone_radii[i], 3) -
         pow(zone_radii[i - 1], 3));
    v1.push_back(v0);
  };
//  proton stuff
  double dd0 = 3.0 * z * oneBypiTimes4 / tot_vol;
  vector<G4double> rod;
  vector<G4double> pf;
  vector<G4double> vz;
  for(G4int i = 0; i < number_of_zones; i++) {
    G4double rd = dd0 * v[i] / v1[i];
    rod.push_back(rd);
    double pff = pf_coeff * pow(rd, one_third);
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
  for(G4int i = 0; i < number_of_zones; i++) {
    G4double rd = dd0 * v[i] / v1[i];
    rod.push_back(rd);
    double pff = pf_coeff * pow(rd, one_third);
    pf.push_back(pff);
    vz.push_back(0.5 * pff * pff / mneutron + binding_energies[1]);
  };
  nucleon_densities.push_back(rod);
  zone_potentials.push_back(vz);
  fermi_momenta.push_back(pf);
//  pion stuff (primitive)
  vector<G4double> vp(number_of_zones, pion_vp);
  zone_potentials.push_back(vp);
}
 else { // a < 4
  number_of_zones = 1;
  zone_radii.push_back(radForSmall);
  G4double vol = 1.0 / piTimes4thirds / pow(zone_radii[0], 3);
  vector<G4double> rod;
  vector<G4double> pf;
  vector<G4double> vz;
  for(G4int i = 0; i < number_of_zones; i++) {
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
  for(G4int i = 0; i < number_of_zones; i++) {
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
  vector<G4double> vp(number_of_zones, pion_vp_small);
  zone_potentials.push_back(vp);  
}; 
nuclei_radius = zone_radii[zone_radii.size() - 1];

}

G4double G4NucleiModel::volNumInt(G4double r1, G4double r2, G4double cu,
                             G4double d1) const {
const G4double au3 = 5.11864;
const G4double epsilon = 1.e-3;
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
  G4int jc1 = G4int(pow(2.0, jc - 1) + 0.1);
  for(G4int i = 0; i < jc1; i++) { 
    r += dr1; 
    fi += r * (r + d2) / (1.0 + exp(r));
  };
  fun = 0.5 * fun1 + fi * dr;
  if(fabs((fun - fun1) / fun) > epsilon) {
    jc++;
    dr1 = dr;
    fun1 = fun;
  }
   else {
    break;
  }; 
}; 
#ifdef DEBUG
if(itry == itry_max) G4cout << " volNumInt-> n iter " << itry_max << G4endl;
#endif
return au3 * (fun + d1 * d1 * log((1.0 + exp(-r1)) / (1.0 + exp(-r2))));

}

G4double G4NucleiModel::volNumInt1(G4double r1, G4double r2, G4double cu2) const {

const G4double epsilon = 1.e-3;
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
  }
   else {
    break;
  }; 
}; 
#ifdef DEBUG
if(itry == itry_max) G4cout << " volNumInt1-> n iter " << itry_max << G4endl;
#endif
return pow(cu2, 3) * fun;

}

void G4NucleiModel::printModel() const {
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

G4InuclElementaryParticle NucleiModel::generateNucleon(G4int type, G4int zone) const {
const G4double one_third = 1.0 / 3.0;
//G4double pmod = getFermiMomentum(type, zone) * pow(inuclRndm(), one_third);
G4double pmod = fermi_momenta[type - 1][zone] * pow(inuclRndm(), one_third);
vector<G4double> mom(4);
pair<G4double, G4double> COS_SIN = randomCOS_SIN();
G4double FI = randomPHI();
G4double pt = pmod * COS_SIN.second;
mom[1] = pt * cos(FI);
mom[2] = pt * sin(FI);
mom[3] = pmod * COS_SIN.first;

return InuclElementaryParticle(mom, type);

}

G4InuclElementaryParticle G4NucleiModel::generateQuasiDeutron(G4int type1, G4int type2,
                   G4int zone) const {
vector<G4double> mom = generateNucleon(type1, zone).getMomentum(); 
vector<G4double> mom1 = generateNucleon(type2, zone).getMomentum();
vector<G4double> dmom(4);
for(G4int i = 1; i < 4; i++) dmom[i] = mom[i] + mom1[i]; 
G4int dtype;
if(type1 * type2 == 1) {
  dtype = 111;
}
 else if(type1 * type2 == 2) { 
  dtype = 112;
}
 else if(type1 * type2 == 4) {
  dtype = 122;
}; 

return G4InuclElementaryParticle(dmom, dtype);

}

partners G4NucleiModel::generateInteractionPartners(
                                      CascadParticle& cparticle) const {

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
vector<G4double> pmom = cparticle.getParticle().getMomentum();

G4double r_in;
G4double r_out;

if(zone == number_of_zones) { // particle is outside 
  r_in = nuclei_radius;
  r_out = 0.0;
}
 else if(zone == 0) { // particle is outside core
  r_in = 0.0;
  r_out = zone_radii[0];
}
 else {
  r_in = zone_radii[zone - 1];
  r_out = zone_radii[zone];
};  

G4double path = cparticle.getPathToTheNextZone(r_in, r_out);
#ifdef DEBUG
G4cout << " r_in " << r_in << " r_out " << r_out << " path " << path << G4endl;
#endif

if(path < -small) { // something wrong
  return thePartners;
}
 else if(fabs(path) < small) { // just on the bounday
  path = 0.0; 
  G4InuclElementaryParticle particle;
  thePartners.push_back(partner(particle, path));
}
 else { // normal case  
  
  vector<G4InuclElementaryParticle> particles;
  G4LorentzConvertor dummy_convertor;
  dummy_convertor.setBullet(pmom, pmass);
  
  for(G4int ip = 1; ip < 3; ip++) { 
    G4InuclElementaryParticle particle = generateNucleon(ip, zone);
    dummy_convertor.setTarget(particle.getMomentum(), particle.getMass());
    G4double ekin = dummy_convertor.getKinEnergyInTheTRS();
    G4double csec = crossSection(ekin, ptype * ip);
#ifdef DEBUG
    G4cout << " ip " << ip << " ekin " << ekin << " csec " << csec << G4endl;
#endif    
    G4double dens = nucleon_densities[ip - 1][zone];
    G4double rat = getRatio(ip);
//    double rat = 1.0;
    G4double pw = -path * dens * csec * rat;
    if(pw < -huge) pw = -huge;
    pw = 1.0 - exp(pw);
#ifdef DEBUG
    G4cout << " pw " << pw << " rat " << rat << G4endl;
#endif    
    G4double spath = path;
    if(inuclRndm() < pw) {
      spath = -1.0 / dens / csec / rat * log(1.0 - pw * inuclRndm());
      if(cparticle.young(young_cut, spath)) spath = path;
#ifdef DEBUG
      G4cout << " ip " << ip << " spath " << spath << G4endl;
#endif    
    };
    if(spath < path) thePartners.push_back(partner(particle, spath));
  };  
#ifdef DEBUG
G4cout << " after nucleons " << thePartners.size() << " path " << path << G4endl;
#endif

  if(cparticle.getParticle().pion()) { // absorption possible
    vector<G4InuclElementaryParticle> qdeutrons;
    vector<G4double> acsecs;
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
    }
     else {
      abs_sec = 0.0;
    }; 
//      abs_sec = 0.0;
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
    }
     else {
      abs_sec = 0.0;
    }; 
//      abs_sec = 0.0;
    tot_abs_csec += abs_sec;
    acsecs.push_back(abs_sec);
    qdeutrons.push_back(nnd);
#ifdef DEBUG
G4cout << " rod1 " << acsecs[0] << " rod2 " << acsecs[1]  
       << " rod3 " << acsecs[2] << G4endl;
#endif

    if(tot_abs_csec > small) {
      double pw = -path * tot_abs_csec;
      if(pw < -huge) pw = -huge;
      pw = 1.0 - exp(pw);
#ifdef DEBUG
    G4cout << " pw " << pw << G4endl;
#endif
      G4double apath = path;
      if(inuclRndm() < pw) 
        apath = -1.0 / tot_abs_csec * log(1.0 - pw * inuclRndm());
      if(cparticle.young(young_cut, apath)) apath = path;  
#ifdef DEBUG
      G4cout << " apath " << apath << " path " << path << G4endl;
#endif
      if(apath < path) { // chose the qdeutron
        G4double sl = inuclRndm() * tot_abs_csec;
        G4double as = 0.0;
        for(G4int i = 0; i < 3; i++) {
          as += acsecs[i];
	  if(sl < as) { 
#ifdef DEBUG
G4cout << " deut type " << i << G4endl; 
#endif	  
            thePartners.push_back(partner(qdeutrons[i], apath));
            break;
	  };
        };
      };    
    };
  };  
#ifdef DEBUG
G4cout << " after deutrons " << thePartners.size() << G4endl;
#endif
  
  if(thePartners.size() > 1) { // sort partners
    for(G4int i = 0; i < thePartners.size() - 1; i++) {
      for(G4int j = i + 1; j < thePartners.size(); j++) {
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

vector<G4CascadParticle> G4NucleiModel::generateParticleFate(G4CascadParticle& cparticle,
     G4ElementaryParticleCollider* theElementaryParticleCollider) {

vector<G4CascadParticle> outgouing_cparticles;

partners thePartners = generateInteractionPartners(cparticle);

if(thePartners.empty()) { // smth. is wrong -> needs special treatment
  G4cout << " generateParticleFate-> can not be here " << G4endl;
}
 else {
  int npart = thePartners.size();
  if(npart == 1) { // cparticle is on the next zone entry
    cparticle.propagateAlongThePath(thePartners[0].second);
    cparticle.incrementCurrentPath(thePartners[0].second);
    boundaryTransition(cparticle);
    outgouing_cparticles.push_back(cparticle);
#ifdef DEBUG
G4cout << " next zone " << G4endl;
    cparticle.print();
#endif
  }
   else { // there are possible interactions
    vector<G4double> old_position = cparticle.getPosition();
    G4InuclElementaryParticle bullet = cparticle.getParticle();
    G4bool no_interaction = true;
    G4int zone = cparticle.getCurrentZone();
    for(G4int i = 0; i < npart - 1; i++) {
      if(i > 0) cparticle.updatePosition(old_position); 
      G4InuclElementaryParticle target = thePartners[i].first; 
#ifdef DEBUG
      if(target.quasi_deutron()) 
        G4cout << " try absorption: target " << target.type() << " bullet " <<
	  bullet.type() << G4endl;
#endif
      G4CollisionOutput output = theElementaryParticleCollider->collide(
                                                        &bullet, &target);
#ifdef DEBUG
      output.printCollisionOutput();
#endif
      vector<G4InuclElementaryParticle> outgoing_particles = 
        output.getOutgoingParticles();
      if(passFermi(outgoing_particles, zone)) { // interaction
	cparticle.propagateAlongThePath(thePartners[i].second);
	vector<G4double> new_position = cparticle.getPosition();
        for(G4int ip = 0; ip < outgoing_particles.size(); ip++) 
	   outgouing_cparticles.push_back(CascadParticle(outgoing_particles[ip],
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
#ifdef DEBUG
G4cout << " partner type " << target.type() << G4endl;
#endif
	
	if(target.nucleon()) {
	  current_nucl1 = target.type();
	}
	 else {
#ifdef DEBUG
	  G4cout << " good absorption " << G4endl;
#endif
	  current_nucl1 = (target.type() - 100) / 10;
	  current_nucl2 = target.type() - 100 - 10 * current_nucl1;
	};   
	  
        if(current_nucl1 == 1) {
	  protonNumberCurrent -= 1.0;
	}
	 else {
	  neutronNumberCurrent -= 1.0;
	}; 
        if(current_nucl2 == 1) {
	  protonNumberCurrent -= 1.0;
	}
	 else if(current_nucl2 == 2) {
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

G4bool G4NucleiModel::passFermi(const vector<G4InuclElementaryParticle>& particles, G4int zone)
                                       {
for(G4int i = 0; i < particles.size(); i++) {
  if(particles[i].nucleon()) {
#ifdef DEBUG
    G4cout << " type " << particles[i].type() << " p " << particles[i].getMomModule()
      << " pf " << fermi_momenta[particles[i].type() - 1][zone] << G4endl;
#endif
#ifdef DEBUG
    if(particles[i].getMomModule() < fermi_momenta[particles[i].type()-1][zone])
     {
        cout << " rejected by fermi: type " << particles[i].type() << 
	" p " << particles[i].getMomModule() << endl;
        return false;
     };
#else
    if(particles[i].getMomModule() < fermi_momenta[particles[i].type()-1][zone])
        return false;
#endif
  };
};

return true; 

}

void NucleiModel::boundaryTransition(CascadParticle& cparticle) {
  int zone = cparticle.getCurrentZone();
  if(cparticle.movingInsideNuclei() && zone == 0) {
    cout << " boundaryTransition-> in zone 0 " << endl;
  }
   else {
    vector<double> mom = cparticle.getMomentum();
    vector<double> pos = cparticle.getPosition();
    int type = cparticle.getParticle().type();
    double pr = 0.;
    double r = 0.;
    for(int i = 0; i < 3; i++) {
      pr += pos[i]*mom[i+1];
      r += pos[i]*pos[i];
    };
    r = sqrt(r);
    pr /= r;
    int next_zone = cparticle.movingInsideNuclei() ? zone - 1 : zone + 1;
    double dv = getPotential(type,zone) - getPotential(type,next_zone);
    double qv = dv*dv - 2.*dv*mom[0] + pr*pr;
    double p1r;
#ifdef DEBUG
cout << " type " << type << " zone " << zone << " next " << next_zone <<
                  " qv " << qv << " dv " << dv << endl;
#endif
    if(qv <= 0.) { // reflection 
      p1r = -pr;
      cparticle.incrementReflectionCounter();
    }
     else { // transition
      p1r = sqrt(qv);
      if(pr < 0.) p1r = -p1r;
      cparticle.updateZone(next_zone);
      cparticle.resetReflection();
    }; 
    double prr = (p1r-pr)/r;  
    for(int i = 0; i < 3; i++) mom[i+1] += pos[i]*prr;
    cparticle.updateParticleMomentum(mom);
  }; 
}

bool NucleiModel::worthToPropagate(const CascadParticle& cparticle) const {
const double cut_coeff = 2.;

  bool worth = true;
  if(cparticle.reflectedNow()) {
    int zone = cparticle.getCurrentZone();
    int ip = cparticle.getParticle().type();
    if(cparticle.getParticle().getKineticEnergy() < cut_coeff*    
      getFermiKinetic(ip,zone)) worth = false; 
  };
  return worth;
}

double NucleiModel::getRatio(int ip) const {
  double rat;
  if(ip == 1) {
#ifdef DEBUG
cout << " current " << protonNumberCurrent << " inp " << protonNumber << endl;
#endif
    rat = protonNumberCurrent/protonNumber;
  }
   else {
#ifdef DEBUG
cout << " current " << neutronNumberCurrent << " inp " << neutronNumber << endl;
#endif
    rat = neutronNumberCurrent/neutronNumber;
  }; 
  return rat;
}

CascadParticle NucleiModel::initializeCascad(
           InuclElementaryParticle* particle) {

const double large = 1000.;

double s1 = sqrt(inuclRndm()); 
double phi = randomPHI();
double rz = nuclei_radius*s1;
vector<double> pos(3);
pos[0] = rz*cos(phi);
pos[1] = rz*sin(phi);
pos[2] = -nuclei_radius*sqrt(1. - s1*s1);

#ifdef DEBUG
CascadParticle cpart(*particle, pos, number_of_zones, large);
cpart.print();
return cpart;
#else
return CascadParticle(*particle, pos, number_of_zones, large);
#endif

}

pair<vector<CascadParticle>, vector<InuclElementaryParticle> >
   NucleiModel::initializeCascad(InuclNuclei* bullet, InuclNuclei* target) {

const double large = 1000.;
const double max_a_for_cascad = 5.;
const double ekin_cut = 2.;
const double small_ekin = 1.e-6;
const double r_large2for3 = 62.;
const double r0forAeq3 = 3.92;
const double s3max = 6.5;
const double r_large2for4 = 69.14;
const double r0forAeq4 = 4.16;
const double s4max = 7.;
const int itry_max = 100;

vector<CascadParticle> casparticles;
vector<InuclElementaryParticle> particles;

// first decide whether it will be cascad or compound final nuclei
double ab = bullet->getA();
double zb = bullet->getZ();
double at = target->getA();
double zt = target->getZ();

if(ab < max_a_for_cascad) {

  double benb = 0.001*bindingEnergy(ab,zb)/ab;
  double bent = 0.001*bindingEnergy(at,zt)/at;
  double ben = benb < bent ? bent : benb;

  if(bullet->getKineticEnergy()/ab > ekin_cut*ben) {

    int itryg = 0;
    while(casparticles.size() == 0 && itryg < itry_max) {
      
      itryg++;
      if(itryg > 0) particles.resize(0);
      
//    nucleons coordinates and momenta in nuclei rest frame
      vector<vector<double> > coordinates;
      vector<vector<double> > momentums;
     
      if(ab < 3.) { // deutron, simplest case
        double r = 2.214 - 3.4208*log(1. - 0.981*inuclRndm());
        double s = 2.*inuclRndm()- 1.;
        double r1 = r*sqrt(1. - s*s);
        vector<double> coord1(3);
        double phi = randomPHI();
        coord1[0] = r1*cos(phi);
        coord1[1] = r1*sin(phi);
        coord1[2] = r*s;   
        coordinates.push_back(coord1);
        for(int i = 0; i < 3; i++) coord1[i] *= -1;
        coordinates.push_back(coord1);
        
	double p;
        bool bad = true;
        int itry = 0;
        while(bad && itry < itry_max) {
          itry++;
	  p = 456.*inuclRndm();
	  if(p*p/(p*p + 2079.36)/(p*p + 2079.36) > 1.2023e-4*inuclRndm() &&
	     p*r > 312.) bad = false;
        };
        if(itry == itry_max) 
          cout << " deutron bullet generation-> itry = " << itry_max << endl;	
        p = 0.0005*p;
#ifdef DEBUG
        cout << " p nuc " << p << endl;
#endif
	vector<double> mom(4);
        pair<double,double> COS_SIN = randomCOS_SIN();
        double FI = randomPHI();
        double P1 = p*COS_SIN.second;
        mom[1] = P1*cos(FI);
        mom[2] = P1*sin(FI);
        mom[3] = p*COS_SIN.first;
        momentums.push_back(mom);
        for(int i = 1; i < 4; i++) mom[i] *= -1;
        momentums.push_back(mom);
      }
       else {
        
	int ia = int(ab + 0.5);
	vector<double> coord1(3);
        bool badco = true;
	int itry = 0;
        
	if(ab < 4.) { // a == 3
	  while(badco && itry < itry_max) {
	    if(itry > 0) coordinates.resize(0);
	    itry++;	    
	    for(int i = 0; i < 2; i++) {
	      int itry1 = 0;
	      double s; 
	      double u;
	      double rho;
	      double fmax = exp(-0.5)/sqrt(0.5);
	      while(itry1 < itry_max) {
	        itry1++;
		s = -log(inuclRndm());
	        u = fmax*inuclRndm();
	        rho = sqrt(s)*exp(-s);
	        if(sqrt(s)*exp(-s) > u && s < s3max) {
                  s = r0forAeq3*sqrt(s);
                  pair<double,double> COS_SIN = randomCOS_SIN();
		  u = s*COS_SIN.second;  
		  double phi = randomPHI();
                  coord1[0] = u*cos(phi);
                  coord1[1] = u*sin(phi);
                  coord1[2] = s*COS_SIN.first;   
                  coordinates.push_back(coord1);
#ifdef DEBUG
		  cout << " i " << i << " r " << sqrt(coord1[0]*coord1[0] +
		    coord1[1]*coord1[1] + coord1[2]*coord1[2]) << endl;
#endif
		  break;
		};
	      };
	      if(itry1 == itry_max) { // bad case
	        coord1[0] = coord1[1] = coord1[2] = 10000.;
                coordinates.push_back(coord1);
		break;
	      };
	    };
            for(int i = 0; i < 3; i++) coord1[i] = - coordinates[0][i] -
	       coordinates[1][i]; 
#ifdef DEBUG
	    cout << " 3  r " << sqrt(coord1[0]*coord1[0] +
		    coord1[1]*coord1[1] + coord1[2]*coord1[2]) << endl;
#endif
	    coordinates.push_back(coord1);  	    
	    
	    bool large_dist = false;
	    for(int i = 0; i < 2; i++) {
	      for(int j = i+1; j < 3; j++) {
	        double r2 = pow(coordinates[i][0]-coordinates[j][0],2) +
	             pow(coordinates[i][1]-coordinates[j][1],2) +
	             pow(coordinates[i][2]-coordinates[j][2],2);
#ifdef DEBUG
		cout << " i " << i << " j " << j << " r2 " << r2 << endl;
#endif
		if(r2 > r_large2for3) {
		  large_dist = true;
                  break; 
		};      
	      };
	      if(large_dist) break;
	    }; 
            if(!large_dist) badco = false;
	  };
	}
	 else { // a >= 4
	  double b = 3./(ab - 2.);
	  double b1 = 1. - b/2.;
	  double u = b1 + sqrt(b1*b1 + b);
	  b = 1./b;
	  double fmax = (1. + u*b)*u*exp(-u);
	  
	  while(badco && itry < itry_max) {
	    if(itry > 0) coordinates.resize(0);
	    itry++;	    
	    for(int i = 0; i < ia-1; i++) {
	      int itry1 = 0;
	      double s; 
	      double u;
	      while(itry1 < itry_max) {
	        itry1++;
		s = -log(inuclRndm());
	        u = fmax*inuclRndm();
	        if(sqrt(s)*exp(-s)*(1.+b*s) > u && s < s4max) {
                  s = r0forAeq4*sqrt(s);
                  pair<double,double> COS_SIN = randomCOS_SIN();
		  u = s*COS_SIN.second;  
		  double phi = randomPHI();
                  coord1[0] = u*cos(phi);
                  coord1[1] = u*sin(phi);
                  coord1[2] = s*COS_SIN.first;   
                  coordinates.push_back(coord1);
#ifdef DEBUG
		  cout << " i " << i << " r " << sqrt(coord1[0]*coord1[0] +
		    coord1[1]*coord1[1] + coord1[2]*coord1[2]) << endl;
#endif
		  break;
		};
	      };
	      if(itry1 == itry_max) { // bad case
	        coord1[0] = coord1[1] = coord1[2] = 10000.;
                coordinates.push_back(coord1);
		break;
	      };
	    };
	    for(int i = 0; i < 3; i++) {
	      coord1[i] = 0.;
	      for(int j = 0; j < ia -1; j++) coord1[i] -= coordinates[j][i];
	    };
	    coordinates.push_back(coord1);   
#ifdef DEBUG
            cout << " last r " << sqrt(coord1[0]*coord1[0] +
		    coord1[1]*coord1[1] + coord1[2]*coord1[2]) << endl;
#endif
	    
	    bool large_dist = false;
	    for(int i = 0; i < ia-1; i++) {
	      for(int j = i+1; j < ia; j++) {
	        double r2 = pow(coordinates[i][0]-coordinates[j][0],2) +
	             pow(coordinates[i][1]-coordinates[j][1],2) +
	             pow(coordinates[i][2]-coordinates[j][2],2);
#ifdef DEBUG
		cout << " i " << i << " j " << j << " r2 " << r2 << endl;
#endif
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
	  cout << " can not generate the nucleons coordinates for a " << ab <<
	           endl;	
          return pair<vector<CascadParticle>, vector<InuclElementaryParticle> >
               (casparticles,particles);
        }
	 else { // momentums
	  double p;
          double u;
	  double x;
	  vector<double> mom(4);
          bool badp = true;
          for(int i = 0; i < ia - 1; i++) {
	    int itry = 0;
	    while(itry < itry_max) {
              itry++;
	      u = -log(0.879853 - 0.8798502*inuclRndm());
	      x = u*exp(-u);
	      if(x > inuclRndm()) {
	        p = sqrt(0.01953*u);
                pair<double,double> COS_SIN = randomCOS_SIN();
		double pt = p*COS_SIN.second;  
		double phi = randomPHI();
                mom[1] = pt*cos(phi);
                mom[2] = pt*sin(phi);
                mom[3] = p*COS_SIN.first;   
                momentums.push_back(mom);
	        break;
	      };
            };
            if(itry == itry_max) {
              cout << " can not generate proper momentum for a " << ab << endl;
              return pair<vector<CascadParticle>, vector<InuclElementaryParticle> >
               (casparticles,particles);
	    }; 
          };
//           last momentum
          for(int i = 1; i < 4; i++) {
	    mom[i] = 0.;
	    for(int j = 0; j < ia -1; j++) mom[i] -= momentums[j][i]; 
	  };
	  momentums.push_back(mom);
	}; 
      }; 
//     coordinates and momentums at rest are generated, now back to the lab;
      double rb = 0.;
      for(int i = 0; i < coordinates.size(); i++) {
        double rp = sqrt(coordinates[i][0]*coordinates[i][0] +
              coordinates[i][1]*coordinates[i][1] +
	            coordinates[i][2]*coordinates[i][2]);
        if(rp > rb) rb = rp;
      };
//    nuclei i.p. as a whole
      double s1 = sqrt(inuclRndm()); 
      double phi = randomPHI();
      double rz = (nuclei_radius+rb)*s1;
      vector<double> global_pos(3);
      global_pos[0] = rz*cos(phi);
      global_pos[1] = rz*sin(phi);
      global_pos[2] = -(nuclei_radius+rb)*sqrt(1. - s1*s1);

      for(int i = 0; i < coordinates.size(); i++) {
        coordinates[i][0] += global_pos[0];
        coordinates[i][1] += global_pos[1];
        coordinates[i][2] += global_pos[2];
      };  

//   all nucleons at rest
      vector<InuclElementaryParticle> raw_particles;
      int ia = int(ab + 0.5);
      int iz = int(zb + 0.5);
      for(int ipa = 0; ipa < ia; ipa++) {
        int knd = ipa < iz ? 1 : 2;
        raw_particles.push_back(InuclElementaryParticle(momentums[ipa],knd));
      }; 
      
      InuclElementaryParticle dummy(small_ekin,1);
      LorentzConvertor toTheBulletRestFrame;
      toTheBulletRestFrame.setBullet(dummy.getMomentum(),
                                               dummy.getMass());
      toTheBulletRestFrame.setTarget(bullet->getMomentum(),bullet->getMass());
      toTheBulletRestFrame.toTheTargetRestFrame();
      particleIterator ipart;
      for(ipart = raw_particles.begin(); ipart != raw_particles.end(); ipart++) {
	vector<double> mom = 
	  toTheBulletRestFrame.backToTheLab(ipart->getMomentum());
        ipart->setMomentum(mom); 
      };

//  fill cascad particles and outgoing particles
      for(int ip = 0; ip < raw_particles.size(); ip++) {
        vector<double> mom = raw_particles[ip].getMomentum();
        double pmod = sqrt(mom[1]*mom[1] + mom[2]*mom[2] + mom[3]*mom[3]);
        double t0 = -(mom[1]*coordinates[ip][0] + mom[2]*coordinates[ip][1] +
           mom[3]*coordinates[ip][2])/pmod;
        double det = t0*t0 + nuclei_radius*nuclei_radius - 
	    coordinates[ip][0]*coordinates[ip][0] - 
	 coordinates[ip][1]*coordinates[ip][1] - 
	 coordinates[ip][2]*coordinates[ip][2];
	 
        double tr = -1.;
	if(det > 0.) {
          double t1 = t0 + sqrt(det);
	  double t2 = t0 - sqrt(det);
	  if(fabs(t1) <= fabs(t2)) {	  
	    if(t1 > 0.) {
	      if(coordinates[ip][2] + mom[3]*t1/pmod <= 0.) tr = t1;
	    };
            if(tr < 0. && t2 > 0.) {
	      if(coordinates[ip][2] + mom[3]*t2/pmod <= 0.) tr = t2;
	    };
	  }
	   else {
	    if(t2 > 0.) {
	      if(coordinates[ip][2] + mom[3]*t2/pmod <= 0.) tr = t2;
	    };
            if(tr < 0. && t1 > 0.) {
	      if(coordinates[ip][2] + mom[3]*t1/pmod <= 0.) tr = t1;
	    };
	  }; 
        };
	if(tr >= 0.) { // cascad particle
          coordinates[ip][0] += mom[1]*tr/pmod;
          coordinates[ip][1] += mom[2]*tr/pmod;
          coordinates[ip][2] += mom[3]*tr/pmod;
	  casparticles.push_back(
	     CascadParticle(raw_particles[ip], coordinates[ip], 
	                                    number_of_zones, large));
        }
         else {
          particles.push_back(raw_particles[ip]); 
        }; 
      };
    };    
    if(casparticles.size() == 0) {
      particles.resize(0);
      cout << " can not generate proper distribution for " << itry_max << "
                     steps " << endl;
    };    
  };
};

#ifdef DEBUG
  cout << " cascad particles: " << casparticles.size() << endl;
  for(int ip = 0; ip < casparticles.size(); ip++)
    casparticles[ip].print();
  cout << " outgoing particles: " << particles.size() << endl;
  for(int ip = 0; ip < particles.size(); ip++)
     particles[ip].printParticle();
#endif
return pair<vector<CascadParticle>, vector<InuclElementaryParticle> >
   (casparticles,particles);

}
