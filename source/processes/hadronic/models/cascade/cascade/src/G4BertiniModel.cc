#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4Nucleus.hh"

#include "G4BertiniModel.hh"

G4BertiniModel::G4BertiniModel() {
verboseLevel = 1  ;
}

G4BertiniModel::~G4BertiniModel(){
  ;
}



void G4BertiniModel::azio(G4double &s, 
	  G4double &c) {
  // description :
  // parameters  : s, c
  // uses        :
  // changes     :
  G4double r1;
  G4double r2;
  G4double sum = 2.0;
  while (sum > 1.0) {
    r1 = dflran();
    r2 = dflran();
    sum = sqr(r1) + sqr(r2); // ysq
  }
  sum /= 2.0;                // (xsq + ysq) / 2
  c = (sum - sqr(r1)) / sum;   // (ysq - xsq) / (xsq + ysq)
  s = (r1 * r2) / sum;         // (2 * x * y) / (xsq + ysq)
  r1 = dflran();
  if (r1 >= 0.5) s = -s;
  return;
}

G4double G4BertiniModel::bovera(const G4double v, 
		const G4double ve) {
  // description :
  // parameters  : v, ve
  // uses        :
  // changes     :
  G4double ver = sqrt(1.0 - sqr(ve) / sqr(v * rcpmv + ve));
  return sqrt(ver * (ver * 1.109 + 0.691) + 0.108) / ver;;
}

void G4BertiniModel::prot(G4double p0, 
	  G4double e, 
	  G4double p, 
	  G4double x, 
	  G4double y, 
	  G4double z) {
  // prot selects from the ranft distribution over all momenta up to p
  // over angles from 0 to 45 degrees
  // to be used for (p, p) and (n, p) inelastic events
  // p0 is momentum of incident particle [GeV/c]
  // p is momentum of secondary nucleon
  // e is energy of secondary nucleon
  // x, y, z are the x, y, z-direction cosines respectively
  G4cout << "  Entering prot" << endl;
  G4double dndpm; // is the maximum of dndp = a5 + a6/pow(p0, a7)
  G4double dndp;
  //  G4double b1;
  G4double b2;
  G4double b3;
  G4double b4;
  G4double b5;
  G4double b6;
  G4double b7;
  G4double r1;
  G4double r2;
  G4double r3;
  G4double r4;
  G4double r5;
  G4double r6;
  G4double csp;
  G4double snp;
  const G4double protonMass = 0.940075;
  const G4double protonMass2 = 0.883741005625; // sqr(protonMass)
  const G4double a1 = 0.885; // a1 - a3 are the ranft coefficients
  const G4double a2 = 0.101;
  const G4double a3 = 4.256;
  const G4double a4 = 0.616849233; // a4 is 45-deg. in rads squared // ::: use CLHEP style (pi/8)^2 *toRad
  const G4double a5 = 0.216470;
  const G4double a6 = 0.672377;
  const G4double a7 = 1.1004872;
  //if (e = 0.0) { // if (e != 0.0) goto L4; // ::: not resolved if correction
  if (e > 1e-6) { //::: ks. orginal  
    r1 = p0;
    G4double  b1 = sqr(r1);
    b2 = sqrt(1.0 + b1/protonMass2);
    dndpm = a5 + a6/pow(r1, a7);
    while(1){
      r2 = G4UniformRand();
      r2 = r1 * r2; // r2 = chosen p - rejection follows
      b3 = sqr(r2);
      b4 = sqrt(1.0 + b3/protonMass2);
      b5 = abs(int(sqr(r2) * b2 - r1 * b4));
      b6 = 1.0 + b2 - r1 * r2/(protonMass2 * b4);
      b7 = -a3 * b3 * a4;
      if (b7 < -50.0) {  // ::: verify if else structure
	b7 = 1.0;
      } else {
	b7 = exp(b7);
	b7 = 1.0 - b7;
      }
      dndp = (pi / a3) * (a1 / r1 + (a2 / b1) * b5) * b6 * b7;
      r3 = dndp/dndpm;
      r4 = G4UniformRand();
      if (r4 < r3) break;
    }
    // rejection completed - accept p = r2. Next follows z
  L6:
    r4 = G4UniformRand();
    r5 = (-1.0/(a3 * b3)) * log(1.0 - r4 * b7);
    z = cos(sqrt(r5));
    azio(csp, snp);
    r6 = sqrt(1.0 - sqr(z));
    x = r6 * csp;
    y = r6 * snp;
    p = r2;
    e = sqrt(b3 + protonMass2) - protonMass;
    G4cout << "  Leaving prot" << endl;
    return;  
  } // end if // L4:
  r2 = p;
  b3 = sqr(r2);
  b7 = -a3 * b3 * a4;
  //  if (b7 < -50.0) goto L5;
  //  b7 = exp(b7);
  //  b7 = 1.0 - b7;
  //  goto L6;
  // L5:
  //  b7 = 1.0;
  //  goto L6; 
  
  if (b7 < -50.0) {// ::: verify this if else structure (see original above)
    b7 = 1.0;
  } else {
    b7 = exp(b7);
    b7 = 1.0 - b7;
  }
  goto L6;
}


G4double G4BertiniModel::exprn() { 
  // description :
  // parameters  : -
  // uses        :
  // changes     :
  G4double expa;
  G4double expb; 
  G4double expao; 
  G4double whole;  
  whole = 0.0;
  while(1) {
    expa = dflran();
    expao = expa; 
    while(1) {
      expb = dflran();
      if (expb >= expa) {
	// random2 > random1
	expa = expao + whole;
	return expa;
      }
      expa = dflran();
      if (expa < expb) continue;
    }
    whole += 1.0;
  }
}


void G4BertiniModel::crdet(G4int nodata, 
	   G4double data[], 
	   G4double ener) {
  // answers stored in crdt  
  G4int i, i1, k, n, l;
  G4int ie;
  // ie = (d1 = *ener / 20.0, (G4int) abs(d1)); //::: fix
  // energy interval 
  G4double univ = (ener - G4double(ie) * 20.0) / 20.0;
  // inpt = 0 if whole interval considered 
  // nodata = data per energy interval 
  for (i = 1; i <= 25; ++i) crdt[i - 1] = 0.0;
  k = nodata * ie + 1;
  if (inpt < 0) {
    goto L40;
  } else if (inpt == 0) {
    goto L50;
  } else {
    goto L80;
  }
  //  98 stop 
 L40:
  G4Exception("crdet-1");
 L50:
  n = nodata;
 L60:
  l = k + nodata;
  i1 = n;
  for (i = 1; i <= i1; ++i) {
    crdt[i - 1] = (data[l] - data[k]) * univ + data[k];
    ++k;
    ++l;
  }
  inpt = 0;
  return;
 L80:
  k = inpt - 1 + k;
  n = 2;
  goto L60;
  // not all parts evaluated 
}

G4double G4BertiniModel::dflran() { 
  // description : u(0,1). First try rejects 0 and 1
  // parameters  : -
  // uses        :
  // changes     :
  G4double value = G4UniformRand();
  while (value == 0.0 || value == 1.0) value = G4UniformRand();
  return value;
}


void G4BertiniModel::pol1(G4double &x, G4double &y){};

G4double G4BertiniModel::massPionCharged = 7.08e12;
G4double G4BertiniModel::massPionZero = 6.84e12;   
G4double G4BertiniModel::oneThird = 1.0 / 3.0;
G4double G4BertiniModel::twoThirds = 1.0 / 3.0;
G4double G4BertiniModel::fourThirds = 4.0 / 3.0;
G4double G4BertiniModel::bindingEnergy(7.0);


