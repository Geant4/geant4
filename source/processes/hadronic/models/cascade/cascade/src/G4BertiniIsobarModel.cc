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
#include "G4BertiniIsobarModel.hh"

G4BertiniIsobarModel::G4BertiniIsobarModel() {
  ;
}

G4BertiniIsobarModel::~G4BertiniIsobarModel(){
  ;
}


void G4BertiniIsobarModel::pi(G4double po, 
			  G4int    particle, 
			  G4double ener, 
			  G4double pmom, 
			  G4double a, 
			  G4double b, 
			  G4double g) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4cout << "  Entering pi" << G4endl;
  G4double ratt[12];
  G4double r2;
  G4double csp;
  G4double snp;
  G4double a1 = 3.386;
  G4double a2 = 4.146;
  G4double a3 = 4.556;
  G4double ta3 = a3 + a3;
  G4double a5 = 9.600;
  G4double a6 = 4.823;
  G4double topi = 6.28318;
  G4double xp   = 0.7853982;
  G4double xp2a3 = sqr(xp) * a3;
  //  G4double sxp  = sqrt(xp); //::: not used in original
  G4double spo  = sqrt(po);
  G4double pos  = sqr(po);
  G4double dmax;
  G4double a4;
  G4double rm;
 
  switch (particle)
    {
    case 1:  // p
      G4cout << "    pi test1" << G4endl;
    case 2:  // n
      a4 = 7.141;
      rm = 0.940;
      break;
    case 3:  // pi+
    case 4:  // pi0
      a4 = 7.141;
      rm = 0.139;
      break;
    case 5:  // pi-
      a4 = 1.853;
      rm = 0.139;
    }
  G4double ps = po * G4UniformRand(); // sample uniformly from po. Find dndp(max) for ps  // ::: dnd ? p
  G4int inx;
  if (po < pmxxx[6]) {
    if (po < pmxxx[4]) {
      if (po > pmxxx[2]) {
	inx = 3;
	if (po > pmxxx[3]) inx = 4;
      } else {
	inx = 1;
	if (po > pmxxx[1]) inx = 2;
      }
      goto L50;
    }
    inx = 5;
    if (po > pmxxx[5]) inx = 6;
  } else if (po < pmxxx[9]) {
    inx = 7;
    if (po >= pmxxx[7]) {
      inx = 8;
      if (po > pmxxx[8]) inx = 9;
    }
  } else {
    inx = 10;
    if (po >= pmxxx[10]) {
      inx = 11;
      if (po > pmxxx[11]) inx = 12;
    }
  }
 L50: 
  if (particle == 5) dmax = dndpim[1][inx] * sqr(po) + dndpim[2][inx] * po + dndpim[3][inx]; // pi-
  else dmax = dndpip[1][inx] * sqr(po) + dndpip[2][inx] * po + dndpip[3][inx];
 L53:
  G4cout << "    pi test2" << G4endl;
  if (ener != 0.0) ps = pmom;
  G4double pss = sqr(ps);
  G4double psa6 = ps * a6;
  G4double f1 = a1 * pss * exp(-a2 * ps / spo);
  G4double f2 = (a4 * pss / po) * exp(-a5 * pss / pos);
  G4double e1 = xp2a3 * ps * spo;
  if (e1 > 50.0) e1 = 50.0;
  G4double e2 = xp * psa6;
  if (e2 > 50.0) e2 = 50.0;
  G4double dpsmx = topi * f1 / (ta3 * ps * spo) * (1.0 - exp(-e1)) + topi * f2 / sqr(psa6) * (1.0 - exp(-e2) * 100.0);
  G4double ratio = dpsmx / dmax;
  r2 = G4UniformRand();
  if (ener == 0.0 && r2 > ratio) {
    ps =  po * G4UniformRand();
    goto L53;
  }
  // for current ps and r2, find angles. Set angle intervals - coarse - test for ratio and r2
  G4double f3 = topi * f1;
  G4double f4 = topi * f2;
  G4double a7 = 1.0 / (ta3 * ps * spo);
  G4double a8 = 1.0 / sqr(psa6);
  G4double f3a7 = f3 * a7;
  G4double f4a8 = f4 * a8;
  a7 *= 2;
  r2 = G4UniformRand();
  G4double r3 = r2;
  G4double exl = -1.0;
  G4double f6l = -1.0;
  ratt[1] = 0.0;
  G4double anglefrs = 0.0;
  G4double anglelas = angle[9];
 L58:
  G4double tu = (anglefrs + anglelas) * 0.5;
  G4double txu = tu * tu / a7;
  if (txu > 50.0) txu = 50.0;
  G4double exu = -exp(-txu);
  G4double f6u = 0.0;
  if ((psa6 * tu) <= 50.0) f6u = exp(-psa6 * tu) * (-psa6 * tu - 1.0);
  G4double f5 = exu - exl;
  G4double f6 = f6u - f6l;
  G4double rat = f3a7 * f5 + f4a8 * f6;
  ratt[2] = rat / dpsmx;
  if (fabs(anglelas - tu) / anglelas < 1.0e-04) goto L80;
  if (r3 <= ratt[2]) goto L68;
  anglefrs = tu;
  goto L58;
 L68:
  anglelas = tu;
  goto L58;
 L80:
  G4double tht = tu; // find angles
  g = cos(tht);
  azio(csp, snp);
  G4double z = sqrt(1.0 - sqr(g));
  a = z * csp;
  b = z * snp;
  ener = sqrt(pss + sqr(rm)) - rm;
  if (ener < 0.0) ener = 1.0e-06;
  pmom = ps;
  return;
}

void G4BertiniIsobarModel::nucnuc(G4double p0, 
			      G4int nofas, 
			      G4int itype) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4cout << "  Entering nucnuc" << G4endl;
  // no nucleon conservation
  G4double val; 
  G4double csa;
  G4double sna;
  G4double csp;
  G4double snp;
  G4double sp[5];
  G4double count(0.0);
  G4double witty(0.03);
  G4double ddd(0.15);
  G4double c45d(0.707);
  G4int ntimes(0);
  G4double erem = einc;
  st0r = einc - 0.139;
  nofas = 0;
  prob(itype, p0);
  sp[1] = pb[1] / (pb[1] + pb[2]);
  sp[2] = 1.0;
  sp[3] = pb[3] / (pb[3] + pb[4] + pb[5]);
  sp[4] = pb[4] / (pb[3] + pb[4] + pb[5]) + sp[3];
  sp[5] = 1.0;
  numberOfNucleus = 0;
 L8:
  val = G4UniformRand();
  G4bool test = true;
  G4int j;
  for (j = 1; j < 6; ++j) {
    if (val < pb[j]) {
      G4cout << "    nucnuc test1" << G4endl;
      test = false;
      break;
    }
  }
  if (test) G4cerr << "nucnuc2" << G4endl;
  j = 1;
  G4double e;
  G4double a;
  G4double b;
  G4double g;
  G4double w;
  G4double p;
  G4double mass;
  G4int i;
  stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
  G4cout << "    nucnuc test2" << G4endl;
  if (w != 0.0) {
    erem -= e + mass;
    if (erem < 0.0) goto L10;
  L23:
    if (++nofas > 60) {  // ::: verify
      if (verboseLevel > 1) {
	G4cout << "nofas greater than 60 in nucnuc" << G4endl;
      }
      if (++ntimes >= 10) G4cerr << "nucnuc1" << G4endl;
      erem = einc;
      st0r = einc - 0.139;
      nofas = 0;
      prob(itype, p0);
      sp[1] = pb[1] / (pb[1] + pb[2]);
      sp[2] = 1.0;
      sp[3] = pb[3] / (pb[3] + pb[4] + pb[5]);
      sp[4] = pb[4] / (pb[3] + pb[4] + pb[5]) + sp[3];
      sp[5] = 1.0;
      numberOfNucleus = 0;
    } else {
      if (i > 2 || ++numberOfNucleus < 3) {
	efas[nofas]   = e * MeV; 
	itxxx[nofas]  = i - 1;
	alpfas[nofas] = a;
	betfas[nofas] = b;
	gamfas[nofas] = g;
	wtfas[nofas]  = fabs(w);
	if (w < 0.0) return;
      } else {
	--nofas;
	j = 2;
	stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
	erem += e + mass;
      }
    }
    goto L8;
  }
  w = 1.0;
  test = false;
  do
    {
      e = 0.0;
      if (i <= 2) {
	prot(p0, e, p, a, b, g);
	if (e > st0r) test = true;
      } else {
	pi(p0, i, e, p, a, b, g);
	if (e > st0r) test = true;
      }
    }
  while (test);
  if (i > 2) {
    j = 4;
    stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
    if (i == 0) goto L8;
    mass = 0.139;
    goto L5;
  }
  j = 4;
  stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
  if (i != 0) {
    mass = 0.0;
  L5:
    erem -= e + mass;
    if (erem < 0.0) goto L10;
    if (++nofas > 60) {
      G4cout << "nofas greater than 60 in nucnuc" << G4endl;
      if (++ntimes >= 10) G4cerr << "nucnuc1" << G4endl;
      erem = einc;
      st0r = einc - 0.139;
      nofas = 0;
      prob(itype, p0);
      sp[1] = pb[1] / (pb[1] + pb[2]);
      sp[2] = 1.0;
      sp[3] = pb[3] / (pb[3] + pb[4] + pb[5]);
      sp[4] = pb[4] / (pb[3] + pb[4] + pb[5]) + sp[3];
      sp[5] = 1.0;
      numberOfNucleus = 0;
    } else {
      if (i > 2 || ++numberOfNucleus < 3) {
	efas[nofas] = e * MeV;
	itxxx[nofas] = i - 1;
	alpfas[nofas] = a;
	betfas[nofas] = b;
	gamfas[nofas] = g;
	wtfas[nofas] = fabs(w);
	if (w < 0.0) return;
      } else {
	--nofas;
	j = 2;
	stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
	erem += e + mass;
      }
    }
  }
  goto L8;
 L10:
  erem += e + mass;
  j = 2;
  stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
  if (numberOfNucleus > 1) {
    val  = G4UniformRand();
    i = 3;
    if (val > sp[3]) i = 4;
    if (val > sp[4]) i = 5;
    if (erem <= einc * ddd) {
      val = G4UniformRand();
      if (val < witty) {
	if ((i > 2) && (erem < 0.139)) {
	  if (erem < 0.0) G4cerr << "nucnuc4" << G4endl;
	  ++count;
	  //	for (G4int i = 0; i < nofas; ++i) efas[i] += erem / nofas * MeV;
	  return;
	}
	azio(csa, sna);
	do {
	  pol1(csp, snp);
	}
	while (csp > c45d);
	a = snp * csa;
	b = snp * sna;
	g = csp;
	i > 2.0 ? mass  = 0.139 : mass = 0.0;
	e = erem - mass;
	w = -1.0;
	goto L23;
      }
    }
    if (i > 2) {
      e = erem - 0.139;
      if (e < 0.0) {
	if (erem < 0.0) G4cerr << "nucnuc4" << G4endl;
	++count;
	for (G4int i = 0; i < nofas; ++i) efas[i] += erem / nofas * MeV;
	return;
      }
      p = sqrt(e * e + 0.278 * e);
      w = 1.0;
      pi(p0, i, e, p, a, b, g);
      j = 3;
      stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
      w = -1.0;
      goto L23;
    }
  } else if (numberOfNucleus > 0) {
    val = G4UniformRand();
    i = 1;
    if (val > sp[1]) i = 2;
  } else {
    numberOfNucleus++;
    val = G4UniformRand();
    val > sp[1] ? i = 2 : i = 1;
    val = G4UniformRand();
    e = erem * val;
    erem -= e;
    p = sqrt(e * (e + 1.88));
    w = 1.0;
    prot(p0, e, p, a, b, g);
    j = 3;
    stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
    if (++nofas > 60) {
      G4cout << "nofas greater than 60 in nucnuc" << G4endl;
      if (++ntimes >= 10) G4cerr << "nucnuc1" << G4endl;
      erem = einc;
      st0r = einc - 0.139;
      nofas = 0;
      prob(itype, p0);
      sp[1] = pb[1] / (pb[1] + pb[2]);
      sp[2] = 1.0;
      sp[3] = pb[3] / (pb[3] + pb[4] + pb[5]);
      sp[4] = pb[4] / (pb[3] + pb[4] + pb[5]) + sp[3];
      sp[5] = 1.0;
      numberOfNucleus = 0;
      goto L8;
    }
    efas[nofas] = e * MeV;
    itxxx[nofas] = i - 1;
    alpfas[nofas] = a;
    betfas[nofas] = b;
    gamfas[nofas] = g;
    wtfas[nofas] = fabs(w);
    val = G4UniformRand();
    i = 1;
    if (val > sp[1]) i = 2;
  }
  e = erem;
  if (e < 0.0) G4cerr << "nucnuc3" << G4endl;
  p = sqrt(sqr(e) + 1.88 * e);
  w = 1.0;
  prot(p0, e, p, a, b, g);
  j = 3;
  stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
  w = -1.0;
  goto L23;
}

void G4BertiniIsobarModel::pinuc(G4double p0, 
			     G4int nofas, 
			     G4int itype) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4cout << "  Entering pinuc" << G4endl;
  G4double val; 
  G4double csa; 
  G4double sna; 
  G4double csp; 
  G4double snp;
  G4double sp[5];
  G4double a;
  G4double p;
  G4double e;
  G4double b;
  G4double g;
  G4double w;
  G4double erem;
  G4double mass;
  G4double temp;
  G4double count(0.0); 
  G4double witty(0.03);
  G4double ddd(0.15);
  G4double c45d(0.707);
  G4int ntimes(0);
  G4int i(1);
  G4int j(1);
  goto L8889;
  while(1) {
    G4cout << "nofas greater than 60 in pinuc" << G4endl;
    if (++ntimes >= 10) G4cerr << "pinuc1" << G4endl;  
  L8889:
    G4double erem = einc + 0.139;
    st0r = einc - 0.139;
    nofas = 0;
    prob(itype, p0);
    G4double cons = pb[1] + pb[2];
    sp[1] = pb[1] / cons;
    sp[2] = 1.0;
    cons = pb[3] + pb[4] + pb[5];
    sp[3] = pb[3] / cons;
    sp[4] = pb[4] / cons + sp[3];
    sp[5] = 1.0;
    numberOfNucleus = 0;
  L8:
    G4cout << "    pinuc test1" << G4endl;
    val = G4UniformRand();
    G4int j;
    for (j = 0; j < 6; j++) { //::: verify
      if (val < pb[j]) goto L2;
    }
    G4cerr << "pinuc2" << G4endl;
  L2:
    stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
    if (w != 0.0) goto L5;
    w = 1.0; 
    mass = 0.0;
    if (i > 2) mass = 0.139;
  L101:
    e = 0.0;
    pi(p0, i, e, p, a, b, g);
    if (e > st0r) goto L101;
    j = 4;
    stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
    if (i == 0) goto L8;
  L5:
    erem =+ - e - mass;
    if (erem < 0.0) goto L10;
  L23:
    nofas++;
    if (nofas <= 60) break;
  } // end while
  if (i > 2) goto L24;
  numberOfNucleus++;
  if (numberOfNucleus > 2) {
    nofas--;
    j = 2;
    stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
    erem =+ e + mass;
    goto L8;
  }
 L24:
  efas[nofas] = e * MeV;
  itxxx[nofas] = i - 1; // ::: it renamed temporarily to itxxx
  alpfas[nofas] = a;
  betfas[nofas] = b;
  gamfas[nofas] = g;
  wtfas[nofas] = fabs(w);
  if (w < 0.0) return;
  goto L8;
 L10:
  erem =+ e + mass;
  j = 2;
  stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
  if (numberOfNucleus > 0) goto L330;
  val = G4UniformRand();
  i = 1;
  if (val > sp[1]) i = 2;
  e = erem;
  p = sqrt(e * (e + 1.88));
  goto L663;
 L330:
  val = G4UniformRand();
  i = 3;
  if (val > sp[3]) i = 4;
  if (val > sp[4]) i = 5;
  if (erem <= (einc * ddd)) {
    val = G4UniformRand();
    if (val < witty) goto L200;
  }
  if (erem < 0.139 && i > 2) goto L25;
  mass = 0.0;
  if (i > 2) mass = 0.139;
  e = erem - mass;
  if (e < 0.0) G4cerr << "mass" << G4endl;
  if (i > 2) {
    p = sqrt(sqr(e) + 1.88 * e);
  } else p = sqrt(sqr(e) + 0.278 * e);
 L663:
  w = 1.0;
  pi(p0, i, e, p, a, b, g);
  j = 3;
  stor(p0, i, e, a, b, g, w, j, erem, p, mass, itype);
  w =  -1.0;
  goto L23;
 L25:
  if (erem < 0.0) G4cerr << "pinuc3" << G4endl;
  count =+ 1.0;
  temp = (erem / nofas) * MeV;
  for (i = 0; i < nofas; i++) efas[i] += temp;
  G4cout << "  Leaving pinuc" << G4endl;
  return;
 L200:
  if (i > 2 && erem < 0.139) goto L25;
  azio(csa, sna);
  do {
    pol1(csp, snp);
  } while (csp > c45d);
  a = snp * csa;
  b = snp * sna;
  g = csp;
  mass = 0.0;
  if (i > 2) mass = 0.139;
  e = erem - mass;
  w = -1.0;
  goto L23;
}

void G4BertiniIsobarModel::prob(G4int itype, G4double pinc) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double cnu[5];
  G4double et[5];
  G4double hold;
  G4double fectr;
  G4double factr;
  G4int inx;
  G4double dene; 
  G4double fune;
  G4double et1;
  G4double et2;
  G4double dnu1;
  G4double dnu2;
  G4double etot;
  const G4double sumas(1.079965);
  const G4double tpmas(1.880150);
  const G4double cpn(1.75);
  const G4double cpp(0.875);
  const G4double witty(2.0); // witty used in (p, p) to get nu, n
  const G4double cnunn(2.0); // cnunp and cnunn used in (n, p) to get nu, p and nu, n
  const G4double cnunp(1.0);
  const G4double abc(0.5);
  const G4double witt2(1.0); // witt2, abc, and ab2 used in (pi+, p) to get nu, p, nu, n, and nu, pi+
  const G4double ab2(1.0);
  const G4double aaa(0.5);
  const G4double wit(1.0); // wit, aaa, abk used in (pi-, p) to get nu, p,  nu, n,  and nu, pi-
  const G4double abk(1.0);
  const G4double cc3(0.0); // cc3 used in (pi+, p) to get nu, pi-
  const G4double cc4(0.0); // cc4 used in (pi-, p) to get nu, pi+
  G4double pinc2 = sqr(pinc);
  for (G4int i = 0; i < 5; i++) pb[i] = 0.0; // initialize pb array
  if (pinc < pz[5]) goto L35; // find pinc position in pz array to get all coefficients
  if (pinc < pz[8]) goto L30;
  inx = 10;
  if (pinc < pz[9]) goto L50;
  inx = 11;
  if (pinc > pz[10]) inx = 12;
  goto L50;
 L30:
  inx = 7;
  if (pinc < pz[6]) goto L50;
  inx = 8;
  if (pinc > pz[7]) inx = 9;
  goto L50;
 L35:
  if (pinc < pz[3]) {
    inx = 5;
    if (pinc > pz[4]) inx = 6;
    goto L50;
  }
  if (pinc > pz[1]) {
    inx = 1;
    if (pinc > pz[0]) inx = 2;
  } else {
    inx = 3;
    if (pinc > pz[2]) inx = 4;
  }
 L50:
  switch (itype) {
  case 1:  // (p, p) collisions
    etot = einc + tpmas;
    cnu[0] = pinc2 *   cpnu[0][inx] + pinc *   cpnu[1][inx]   + cpnu[2][inx];
    cnu[2] = pinc2 * cpipnu[0][inx] + pinc * cpipnu[1][inx] + cpipnu[2][inx];
    cnu[4] = pinc2 * cpimnu[0][inx] + pinc * cpimnu[1][inx] + cpimnu[2][inx];
    et[0] =  pinc2 *    cpk[0][inx] + pinc *    cpk[1][inx]    + cpk[2][inx];
    et[2] =  pinc2 *  cpipk[0][inx] + pinc *  cpipk[1][inx]  + cpipk[2][inx];
    et[4] =  pinc2 *  cpimk[0][inx] + pinc *  cpimk[1][inx]  + cpimk[2][inx];
    dnu1 = cnu[0];
    dnu2 = witty - cnu[0];
    et2 = et[0] * dnu2/cnu[0];
    et1 = et[0];
    if (pinc < 18.0) goto L220;
    cnu[1] = dnu2;
    et[1] = et2;
    fune = etot - et[0] - et[1] - et[2] - et[4];
    dene = et[2]/cnu[2];
  case 2:  // (n, p) collisions
    etot = einc + tpmas;
    cnu[0] = pinc2 *   cpnu[0][inx] + pinc *   cpnu[1][inx]   + cpnu[2][inx];
    cnu[2] = pinc2 * cpipnu[0][inx] + pinc * cpipnu[1][inx] + cpipnu[2][inx];
    cnu[4] = pinc2 * cpimnu[0][inx] + pinc * cpimnu[1][inx] + cpimnu[2][inx];
    et[0] =  pinc2 *    cpk[0][inx] + pinc *    cpk[1][inx]    + cpk[2][inx];
    et[2] =  pinc2 *  cpipk[0][inx]  + pinc * cpipk[1][inx]  + cpipk[2][inx];
    et[4] =  pinc2 *  cpimk[0][inx]  + pinc * cpimk[1][inx]  + cpimk[2][inx];
    dnu2 = cnunn - cnunp;
    et2 = et[0] * dnu2 / cnu[0];
    et1 = et[0] * cnunp / cnu[0];
    dnu1 = cnunp;
    if (pinc < 18.0) goto L220;
    cnu[1] = dnu2;
    et[1]  = et2;
    et[0]  = et1;
    cnu[0] = dnu1;
    dene = et[2] / cnu[2];
    fune = etot - et[0] - et[1] - et[2] - et[4];
  case 3:  // (pi+, p) collisions
    etot   = einc + sumas;
    cnu[0] = abc;
    cnu[1] = witt2 - abc;
    cnu[2] = pinc2 * cpipnu[0][inx] + pinc * cpipnu[1][inx] + cpipnu[2][inx];
    cnu[4] = pinc2 * cpimnu[0][inx] + pinc * cpimnu[1][inx] + cpimnu[2][inx];
    et[2]  = pinc2 *  cpipk[0][inx] + pinc *  cpipk[1][inx] +  cpipk[2][inx];
    et[4]  = pinc2 *  cpimk[0][inx] + pinc *  cpimk[1][inx] +  cpimk[2][inx];
    dene   = et[2] / cnu[2];
    cnu[2] = cnu[2] + ab2;
    et[0]  = dene * cnu[0];
    et[1]  = dene * cnu[1];
    et[2]  = dene * cnu[2];
    et[4]  = et[4] * (1.0 + cc3 / cnu[4]);
    cnu[4] = cnu[4] + cc3;
    fune   = etot - et[0] - et[1] - et[2] - et[4];
  case 4:  // (pi0, p) collisions
    // you cant get here from there 
    G4cerr << "prob" << G4endl;
  case 5:  // (pi-, p) collisions
    etot = einc + sumas;
    cnu[0] = aaa;
    cnu[1] = wit - aaa;
    cnu[2] = pinc2 * cpipnu[0][inx] + pinc * cpipnu[1][inx] + cpipnu[2][inx];
    cnu[4] = pinc2 * cpimnu[0][inx] + pinc * cpimnu[1][inx] + cpimnu[2][inx];
    et[2]  = pinc2 *  cpipk[0][inx] + pinc *  cpipk[1][inx] +  cpipk[2][inx];
    et[4]  = pinc2 *  cpimk[0][inx] + pinc *  cpimk[1][inx] +  cpimk[2][inx];
    dene   = et[2] / cnu[2];
    et[0]  = dene * cnu[0];
    et[1]  = dene * cnu[1];
    hold   = et[4] / cnu[4];
    cnu[4] += abk;
    et[4]  = hold * cnu[4];
    et[2]  *= (1.0 + cc4 / cnu[2]);
    cnu[2] += cc4;
    fune   = etot - et[0] - et[1] - et[2] - et[4];
  }
  cnu[3] = fune / dene;
  pb[0]  = cnu[0]; 
  for (i = i; i <5; i++) pb[i] = cnu[i] + pb[i - 1]; 
  for (i = 1; i < 6; i++) pb[i] /= pb[5];
  return;
 L220:
  cnu[1] = cpn - cpp;
  et[1] = et[0] * cnu[1] / cnu[0];
  et[0] = et[0] * cpp / cnu[0];
  cnu[1] = cpp;
  fectr = (18.0 - pinc) / 15.5;
  factr = 1.00001 - fectr;
  cnu[1] = cnu[1] * fectr + dnu2 * factr;
  cnu[0] = cnu[0] * fectr + dnu1 * factr;
  et[1] =   et[1] * fectr +  et2 * factr;
  et[0] =   et[0] * fectr +  et1 * factr;
  dene = et[2] / cnu[2];
  fune = etot - et[0] - et[1] - et[2] - et[4];
  cnu[3] = fune / dene;
  pb[0]  = cnu[0]; 
  for (i = 1; i < 5; i++) pb[i] = cnu[i] + pb[i - 1]; 
  for (i = 0; i < 5; i++) pb[i] /= pb[4];
}

void G4BertiniIsobarModel::rout1() {
  // decription: set nuclei potential
  // parameters:
  // uses:
  // changes:
  out[10] = zee;
  G4double value2 = pow(zee, twoThirds);
  for (G4int i = 4; i < 7; i++) {
    G4cout << "    rout1 test1" << i << G4endl; 
    space[i + 2] = out[i  ] * out[10];
    space[i + 5] = out[i+3] * value2 + 7.0;
  }
  // scaled p per cc and potential p well depth (MeV) in each region
  out[11] = amasno - out[10]; // number of neutrons n,  stored
  value2 = pow(out[11], twoThirds);
  for (i = 4; i < 7; i++) {
    G4cout << "    rout1 test2 " << i << G4endl; 
    space[i - 4] = out[i    ] * out[11];
    space[i - 1] = out[i + 3] * value2 + 7.0;
  }
  // scaled n. per cc and pot. neut. well depth (MeV) in each region
  for (i = 0; i < 3; i++) {
    hvn[i]  = 0.5 * space[i + 3];
    hvp[i]  = 0.5 * space[i + 9];
    awd[i]  = hvn[i] + hvp[i];
    fvnp[i] = 0.5 * awd[i];
    vnvp[i] = space[i + 3] - space[i + 9];
    pmac[i] =  vnvp[i]  - hvn[i];
    ppan[i] = -vnvp[i]  - hvp[i];
    thpn[i] =   hvp[i] - vnvp[i];
    ffptfn[i] = - vnvp[i] + fvnp[i];
    tffn[i] = space[i + 9] - fvnp[i];
    tffp[i] = vnvp[i] + tffn[i];
  }
  pppda = (2.0 * zee) / (zee + amasno - 1.0);
  ppmda = (2.0 * out[12]) / (amasno + out[12] - 1.0);
  // G4double ppnda = (2.0 * zee * out[12]) / (amasno * amasno - amasno); //::: is it used anywhere
  ppnna = (out[12] * out[12] - out[12]) / (out[12] * out[12] + zee * zee - amasno);
  // pi absorption calculated for each region  1/2 nwd,  (-pnan,  -ppac)
  // 1/2 pwd,  (-pnap,  -pnac) average well depth,  1/4 average well depth, 
  // (n-p) well depth,  (-vpvn),  1/2 nwd -pwd,  (pmap,  -vphn),  1/2 pwd -nwd, 
  // (-vnhp),  3/2 pwd -nwd,  5/4 pwd -3/4 nwd,  3/4 pwd -1/4 nwd, 
  // 3/4 nwd -1/4 pwd,  probility (pi+, d) absorption,  probility (pi-, d) absorption
  // probility  (pi-, d) abs, prob (pi, nn) abs rather than pp
  G4int k = 15;
  for (i = 4; i < 7; i++) {
    out[k] = space[i + 6] + einc;
    out[k + 3] = space[i] + einc;
    --k;
  }
  // total kinetic energy in MeV incident p(n) particle in each region
  out[30] = zee / out[4] * 1.4412 * 10e-13;
  // Coulomb potential at surface in MeV. Conversion factor = MeV-cm/p (bg33)
  if (ctofe == 0) ctofe = out[30]; //if ctofe = 0, then equate it to 1/2 potential energy at surface
  for (i = 0; i < 3; i++) {
    cfepn[i + 3] = space[i + 3] + ctofen;
    cfepn[i    ] = space[i + 9] + ctofe;
  }
  // bg33p-cutoff energies in each region for n(p)
  in = 0;
  G4double value1 = 6.28318531e+10 * pow(0.119366207, oneThird); // Fermi momenta / cm. pf = 2*PI*((3/8*PI)^(1/3))*e10
  for (i = 1; i < 4; i++) {
    fmpn[i  ] = value1 * pow(space[i + 6], oneThird);
    fmpn[i+3] = value1 * pow(space[i    ], oneThird);
  }
  // Fermi momenta / cm. of p(n)
  for (i = 1; i < 5; i++) rands[i] = randi[i];
  G4cout << "  Leaving rout1 " << G4endl;
  return;
}

void G4BertiniIsobarModel::rout2(G4double *t) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4cout << "  Entering rout2" << G4endl;
  i1 = 1;
  value2 = einc + space[12];
  G4double ans;
  ans = bovera(value2, massPionCharged);
  space[14] =   0.20e-24 * ans; // (pi+, p)
  space[15] =  0.023e-24 * ans; // (pi-, p) el
  space[16] = 0.0451e-24 * ans; // (pi-, p) ex;
  G4cout << "    rout2 test1 " << G4endl; 
  if (value1 <= 100.0) {
    fmax[1] = space[14];
    fmax[2] = space[15];
    fmax[3] = space[16];
    s[1] = 0.0;
    s[2] = 0.0;
    G4cout << "    rout2 test2" << G4endl;
  } else goto L130;
 L12500:
  crdet(1, t, einc); // (pi+, p) absorption cross-section
  space[17] = crdt[1];
  if (i1 > 0) fmax[4] = space[17];
  return;
 L130:
  if (value2 > 2600.0) {
    i1 = 0;
    return;
  } else {
    G4cout << "    rout2 test4" << G4endl; 
    space[17] = 0.0;
    if (value2 <= 220.0) {
      s[1] = 0.75e-27 * ans; // (pi+, p) s.p. 400 MeV
      s[2] = 4.7e-27 * ans; // (pi-, p) s.p. 400 MeV
      goto L146;
    }
  }
  if (value2 <= 400.0) {
    space[14] = 0.20e-24;
    space[16] = 0.0451e-24;
    s[1] =  7.8e-27 * ans;
    s[2] = 21.8e-27 * ans;
  } else {
    goto L160;}
  // 660 MeV;
 L146:
  if (einc < 360.0) {
    goto L12500;
  } else return;
 L160:
  if (value2 <= 500.0) {
    space[14] = 0.113e-24; // 250 MeV;
    space[15] = 20.5e-27 * ans; // 620 MeV
    space[16] = 27.7e-27; // 250 MeV
    s[1] = 13.8e-27 * ans;
    s[2] = 24.4e-27 * ans; // 800 MeV
    i1 = -1;
    goto L146;
  }
  s[2] = 30.4e-27 * ans;
  space[15]= 26.3e-27 * ans;
  // 900
  if (value2 <= 600.0) {
    // 940, 325
    space[14] = 53.0e-27;
    // 325
    // value2 <= 600
    // 325
    space[16] = 16.2e-27;
    s[1] = 15.2e-27 * ans;
    // 940
    return;
  }
  if (value2 <= 800.0) { // 1200, 400
    space[14] = 33.0e-27; // 400
    space[16]= 12.0e-27 * ans; // 400
    s[1] = 20.9e-27 * ans; // 1200
    return;
  }
  space[14] = 19.3e-27 * ans; // 1300
  space[16] = 8.2e-27 * ans; // 540
  s[1] = 23.3e-27 * ans; // 1400
  s[2] = 30.4e-27 * ans;
  // sigma(a) 900 + 2 mb for future correction
  return;
  // values of fi for both incident and cascade particles for pi+-, single production.  
  // s[1] = n (s.p.,  s[2] = p) s.p.,  space[14] = n (s,  space[15]= p) d.s.,  
  // space[16]= p)  abs,  space(17) = p) abs  (ppac);
}

void G4BertiniIsobarModel::rout6() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  if (i3 < 0) {
    abz = 1.0;
    medium = static_cast<G4int>(clsm);
    knot = not;
    value1 = G4UniformRand();
  } else goto L345;
  if (isw[11] == 0 ) {
    if (value1 >=  ppmda) {
      i3 = 1;
      return;
    } else goto L330; // probability (pi-, d) abs.
  } else goto L325;
  return;
 L325:
  if (value1 >= pppda) {  // probability (pi+, d) abs.
    i3 = 1;
    return;
  }
 L330:
  if (isw[11] == 0 ) {
    it = 13;
    absec = pmac[medium];
    goto L345;
  } else {
    it = 14;
    absec = - hvn[medium];
  }
 L345:
  strkp = - 1.0;
  i1 = 0;
  i2 = medium;
  // bb(); // ::: 
  i3 = 0;
  return;
}
	  
void G4BertiniIsobarModel::rout7() {
  G4cout << "  Entering rout7" << G4endl;
  i3 = 0;
  if (curr[1] > 3.0) { 
    if (curr[1] > 5.0) {
      goto L135;
    } else {
      if (curr[1] < 5) {
	goto L385; 
      } else goto L375;
    }
    // p, n not permitted
  } else {
    if (curr[1] < 3.0) {
    L135:
      i3 = - 1;
    } else goto L380;
  }
  return;
 L375:
  it    = 7;
  ifca  = 5; // pi - [5]
  absec = pmac[medium]; // (pi-, pp) abs.  tyor = pmapp(20021);
  goto L400;
 L380:
  it    = 10;
  ifca  = 3;
  // tyor = ppan(2000) (pi+, nn) abs. energy correction pi+;
  absec = ppan[medium];
  goto L405;
 L385:
  value1 = G4UniformRand();
  if (value1 < ppnna) {
    it = 8;
    ifca = 4; // pnann(20015) = tyor (pi-, nn) abs  pi
    absec = - hvn[medium];
    goto L405;
  }
  it = 9;
  ifca = 2; // pnapp(20011) = tyor (pi-, pp) abs. pi
  absec = - hvp[medium];
 L400:
  strkp = - 1.0;
  energy[1] = wkrpn[medium] * rcpmv + massParticle[1];
  goto L410;
 L405:
  strkp = - 2.0;
  energy[1] = wkrpn[medium + 3] * rcpmv + massParticle[1];
 L410:
  if (inc == 0) {
    // p1clc(); // :::
  } else {
    // p1cli(); // :::
  }
  return;
}

void G4BertiniIsobarModel::rout7a() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  {
    i1 = - 1;
    spisom();
    if (ifca == 1 || ifca == 3 || ifca == 4) value1 = space[medium + 3] - 7.0;
    if (ifca == 2 || ifca == 5) value1 = space[medium + 9] - 7.0;
    if (value1 < 0.0){
      i3 = 2;
      return;
    }
  } while ((value1 * 2.0 * rcpmv) <= energy[2]);
  massParticle[3] = massNucleon;
  massParticle[2] = 2.0 * massNucleon;
  energy[2] = massParticle[2] + energy[2];
  value1 = ex;
  if (medium < 2) { 
    i3 = 3;
    return;
  } else if (medium == 2) {
    i3 = 4;
    return;
  }
  if (inc != 0) { 
    if (isw[1] != 0) {
      i3 = 5;
      return;
    }
    if (isw[2] != 0) {
      i3 = 6;
      return;
    }
    i3 = 9;
    return;
  }
  if (isw[1] == 0 ) {
    i3 = 7;
    return;
  }
  if (isw[2] == 0) { 
    i3 = 8;
    return;
  }
  i3 = 1;
  return;
}

void G4BertiniIsobarModel::rout8() {
  // decription:
  // parameters:
  // uses:
  // changes:
  i3 = 1;
  if (iv >= 0) {
    if (value1 > value2) goto L600;
  }
  if (isw[3] != 0) goto L595;
  ifc = 7 + ifcc;
  // 7 = bg6e(2461) 8 = bg6ia(4026) ntnt(21626) bg48x(12762)= 19
  if (in != 0) {
    i3 = 2;
    return;
  }
  c[3] = d[2];
  for (;;) {
    i3 = 3;
    return;
  L595:
    ifc = 8 + ifcc;
    if (in != 0) {
      i3 = 4;
      return;
    }
    c[3] = d[2] + d[3] + d[4];
  }
 L600:
  signex();
  return;
}


void G4BertiniIsobarModel::rou13() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4cout << "  Entering rou13" << G4endl;
  i3 = 0;
  if (iv >= 0) goto L705;
  if (abz == 0) goto L705;
  if (ifca == 2) goto L703;
  if (ifca > 2) goto L702;
 L701:
  in = 0;
 L350:
  i3 = 1;
  return;
 L705:
  signex();
  if (ifc <= 12) goto L706;
  else goto L707;
 L702:
  if (ifca == 6) {
    goto L701;
  } else if (ifca > 6) goto L710;
 L703:
  in = 0;
  i3 = - 1;
  G4cout << "  Leaving rou13" << G4endl;
  return;
 L710:
  if (ifca < 8) goto L701;
  in = 1;
  goto L350;
 L706:
  in = 0;
  return;
 L707:
  if (ifc <= 18) {
    in = - 1;
    return;
  }
  in = 1;
  return;
}
	  
void G4BertiniIsobarModel::rou14() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4cout << "  Entering rou14" << G4endl;
  G4int m;
  if (i3 < 0) goto L861;
  if (i3 == 0) goto L850;
  G4cout << "    rou14 test1" << G4endl;
  if (i1 >= 0 ) goto L740;
  i1 = 0;
  value1 = (energy[3] - massParticle[3])/rcpmv;
  if (abz != 0) value1 =+ absec;
  if (pt[2] < 2.0) goto L710;
  if (pt[2] == 2.0) goto L685;
  else goto L740;
  // pt[2] = 1 = p,  pt[2] = 2 = n, pt[2] = 3, 4, 5 = pi
 L685:
  i3 = 2;
 L686:
  goto L1000;
 L710:
  i3 = 1;
  goto L686;
 L740:
  pinst();
  if (i1 != 0 && verboseLevel > 0) {
    G4cout << "rou14:just called pinst: i1 = " << i5 << " i3 = " << i5 << G4endl;
  }
  if (i1 == 0) goto L800;
 L135:
  i3 = 3;
  goto L686;
 L800:
  i1 = 0;
  m = static_cast<G4int>(pt[2]);
  value2 = value1;
  if (m >= 3) goto L820;
 L805:
  if (eco[m] < value2) goto L815;
  if (eco[m] > value2) goto L810;
  if (eco[m] <= 0.0) goto L135;
 L810:
  pt[i1 + 3]= 0.0;
  pnbc[m] = pnbc[m] + 1.0;
  goto L834;
 L815:
  pt[i1 + 3] = value2;
  if (i1 == 0) goto L835;
  else goto L860;
 L820:
  ccofe = clcfe;
  if (m >= 4) goto L822;
  if (strkp < -2.0) goto L826;
  if (strkp == -2.0) goto L824;
  else goto L826;
 L822:
  if (strkp <= -2.0) goto L826;
  ccofe = clcfe - ctofe + ctofen;
  goto L826;
 L824:
  ccofe = clcfe + ctofe - ctofen;
 L826:
  if (value2 <= ccofe) goto L810;
  if (strkp > -2.0) goto L815;
  pt[3] = value1 - space[medium + 3] + space[medium + 9];
 L834:
  if (i1 != 0) goto L845;
 L835:
  m = static_cast<G4int>(pt[14]);
  if (m >= 3) goto L135;
  value2 = com;
  i1 = 12;
  goto L805;
 L845:
  if (pt[3] != 0.0) goto L860;
 L850:
  capunp(); 
  if (i1 != 0 && verboseLevel > 0) G4cout << "rou14: just called punp: i1 = " << i5 << G4endl;
  if (i1 < 0 ) goto L135; // error
  if ( i1 > 0 ) goto L4415; // piscc(6607)
  i3 = 4; // end of record
  goto L686;
 L4415:
  i3 = 5;
  goto L686;
 L860:
  collm(- 1);
 L861:
  if (ke > 0) goto L1000;
  castpr(); 
  if (i1 != 0) goto L135;
  else goto L850;
 L1000:
  return;
}
	  
void G4BertiniIsobarModel::rou17() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double value3;
  G4double unive;
  G4double a;
  G4double com1;
  G4double ftr;
  if (i3 > 0) goto L2040;
  pt[38]= 0.0; // 1-alpha  part. 6 + 1
  value1 = rlke - 180.0;
  crdet(1, t, value1);
  com2 = crdt[1];
  ftr = massNucleon * rlke * 2.0 * rcpmv + 2.9877156e27; // sqr(e) = sqr(min) + ncnms * rlke * 2 * rcpmv
  univer = sqrt(ftr); // e
 L2005:
  value2 = G4UniformRand();
  G4cout << "    rout17 test2" << G4endl;
  com = value2 * com2; // r-prime
  cagene(b[1]);
  com1 = (sqr(com) + ftr - .501264e26)/(2.0 * univer); // m1r prime)**2 + sqr(e) - 2(massPionCharged)/2e = e alpha
  a = sqr(com1) - sqr(com);
  if (a >= 0) goto L2009;
  pacnt =+ 1.0;
 L2009:
  unive = ((univer - com1) * com1/univer) * sqrt(a);
  // ((e beta * e alpha * p alpha)/e) = f(m, tr)
  crdet(1, rr, value1);
  // (pi, N) fmax(rlke) isobar sampling s.p.
  com1 = G4UniformRand();
  if ((unive/crdt[1]) < com1) goto L2005; // random number <= f(m, tr) / fmax(tr)
  G4cout << "    rout17 test3" << G4endl;
  //  cangid; //:::fix
  massParticle[3] = com;
  massParticle[4] = massPionZero;
  pt[2] = 3.0;
  pt[4] = massPionZero;
  pt[14] = 3.0;
  pt[16]= massPionZero;
  pt[26]= 1.0;
  pt[28]= massNucleon;
  if (isw[9] == 0) {
    if (isw[10] != 0) {
      goto L2030;
    } else {
      goto L2025;
    } 
  }
  if (isw[10] != 0) {
  L135:
    i3 = -1;
    return;
  L2025:
    value1 = 0.4;
    value2 = twoThirds;
    value3 = 0.0;
    goto L2037;
  L2030:
    crdet(2, w, value1); // (pi+-, p) fractional final state with recl. pi1 pi0 l.e.
    value3 = oneThird;
    goto L2036;
  }
  crdet(2, g, value1);
  value3 = strkp; // (pi-, p) fractional final state with recl. pi1 pi0 l.e.
 L2036:
  value1 = crdt[1];
  value2 = crdt[2];
 L2037:
  // g4bertini->alpha(); //:::
 L2040:
  //  g4bertini->ecpl();  //:::
  if (i1 > 0) goto L135; 
  coll(- 1);
  if (col[15] != 0) goto L135; 
  if (pt[38] != 0) { 
    i3 = 0;
    return;
  }
  pt[39] = 0.0;
  pt[3] = ((energy[4] - massParticle[4])/rcpmv) + pt[3];
  i3 = 1;
}

void G4BertiniIsobarModel::rou18() {
  G4cout << "  Entering rou18" << G4endl;
  G4int i;
  G4int j;
  G4int k;
  G4int m;
  switch (i3) {
  case 1:  
    goto L4010;
  case 2:  
    goto L4015;
  case 3:  
    goto L4035;
  case 4:  
    goto L4095;
  case 5:  
    goto L4071;
  case 6: // ::: check if 6 or all other cases  
    goto L4080;
  }
 L4010:
  i = 3;
  col[15] = 1.0;
  k = 27;
  goto L4020;
 L4015:
  i = 3;
  col[15] = 4.0;
  k = 15;
 L4020:
  G4cout << "    rou18 test1" << G4endl;
  pnidk[1] = massParticle[i];
  j = i;
  G4int l;
  for (l = 2; l < 5; l++) {
    pnidk[l] = pxyz[j];
    j += 4;
  }
  pnidk[5] = energy[i];
  pnidk[6] = pt[k - 11];
  //  idk; //::: fix
  if (k != 27) goto L4031;
  pt[15] = pt[15] + ((pnidk[12] - pnidk[6])/rcpmv);
 L4031:
  pt[k] = pt[k] + ((pnidk[13] - massNucleon)/rcpmv);
  i3 = 1;
 L2057:
  iv = k;
  G4cout << "  Leaving rou18" << G4endl;
  return;
 L4035:
  k = 3;
  col[15] = 2.0;
  if (pt[2] < 3.0) goto L4071;
 L4039:
  if (pt[k] > 2500.0) {
    i3 = 5;
    goto L2057;
  }
  if (pt[k] <= 0.0) goto L4050;
  ccofe = eco[1];
  if (pt[k - 1] >= 4.0) {
    ccofe =- ctofe + ctofen;
  }
  if (pt[k] > ccofe) {
  L4050:
    m = static_cast<G4int>(pt[k - 1]); // ::: int added
    pnbc[m] = pnbc[m] + 1.0;
    pt[k] = 0.0;
    i3 = 3;
    goto L2057;
  }
  if (k > 3) goto L4095;
 L4071:
  col[15] = 3.0;
  k = 15;
  if (pt[14] > 2.0) goto L4039;
  i3 = 2;
  goto L2057;
 L4080:
  for (m = 5; m < 8; m++) {
    G4int l = 14;
    pt[m] = pnidk[l];
    pt[m + 12] = pnidk[l + 3];
    ++l;
  }
  pt[11] = pnidk[12];
  pt[12] = pnidk[6];
  i = 4;
  k = 39;
  col[15] = 5.0;
  goto L4020;
 L4095:
  i1 = 3;
 L4100:
  k = 12 * i1 - 33;
  if (i1 == 4) goto L4115;
  if (i1 > 4) goto L4120;
  goto L4120; 
  i2 = - 1;
  goto L4130;
 L4115:
  i2 = 0;
  goto L4130;
 L4120:
  if (i1 < 5) goto L4115; 
  if (i1 != 5) {
    i3 = 4;
    goto L2057;
  }
  i2 = 1;
 L4130:
  if (pt[k] != 0.0) {
    //  g4bertini->pstor();  //:::
  }
  i1++ ;
  goto L4100;
}

void G4BertiniIsobarModel::rou19() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int k;
  G4int m;
  pt[3] = pt[3] + ((pt[11] - pt[12]) / rcpmv);
  // collision allowed
  k = 3;  
 L4135:
  if (pt[k] <= 2500.0) goto L4150;
  i3 = 1;
  return;
 L4150:
  if (pt[k] <= 0.0) goto L4160;
  ccofe = eco[1];
  if (pt[k - 1] < 4.0) goto L4157;
  ccofe =- ctofe + ctofen;
 L4157:
  if (pt[k] > ccofe) goto L4200;
 L4160:
  pt[k] = 0.0;
  if (pt[k - 1] == 3.0) goto L4170;
  if (pt[k - 1] > 3.0) goto L4165; 
  i3 = - 1;
  return;
 L4165:
  if (pt[k - 1] > 5.0) {
    i3 = -1;
    return;}
 L4170:
  m = static_cast<G4int>(pt[k - 1]);
  pnbc[m] = pnbc[m] + 1.0;
  goto L4185;
 L4175:
  i2 = 2;
 L4176:
  i1 = (k/12) + 3;
  // g4bertini->pstor();  //:::
 L4185:
  if (k < 15) goto L4190;
  if (k == 15) goto L4195; 
  else goto L4210;
 L4190:
  k = 15;
  if (pt[15] > 0.0) goto L4175; 
 L4195:
  k = 27;
  pt[27]= pt[27] + ((pnidk[12] - pt[k + 1]) / rcpmv);
  goto L4135;
 L4200:
  if (k < 15) goto L4175; 
  i2 = 0;
  goto L4176;
 L4210:
  if (k < 27) {
    i3 = -1;
    return;
  } 
  if (k <= 27) {
    if (pt[39] > 0) goto L4220;
  }
  i3 = 0;
  return;
 L4220:
  i2 = 1;
  k = 39;
  goto L4176;
}

void G4BertiniIsobarModel::rou21(G4double *v, 
			     G4double *w, 
			     G4double *x, 
			     G4double *y, 
			     G4double *z) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double a;
  G4double com1;
  G4double com4;
  G4double sn;
  G4double ftr;
  G4double ans;
  G4cout << "  Entering rou21" << G4endl;
  G4double value2 = rlke * 4.81633308e24 + 9.0554256e27; // sqr(e(tr)) = rlke * rcpmv * 2 * ncms + 4 * sqr(ncms);
  G4double value3 = sqrt(value2);
  switch (i3) {
  case 1:  
    goto L4330;
  case 2:  
    goto L4360;
  case 3:  
    goto L4365;
  case 4:  
    goto L4410;
  }
 L4330:
  isw[12] = 0;
 L4331:
  pt[38] = 0.0;
  i1 = 0;
  ans = rlke;
 L4333:
  value1 = ans - 300.0;
  crdet(1, v, value1); // (N, N) f(tr) isobar sampling
  ftr = crdt[1];
 L4335:
  sn = G4UniformRand();
  com = sn * ftr; // r prime = f(tr) * random
  cagene(w[1]); // (N, N) mass of isobar s.p. m(r prime)
  // ::: if (i1) 4370, 4336, 4375
  if (i1 < 0) goto L4370;
  if (i1 > 0) goto L4375;
  com1 = (sqr(com) - sqnm + value2)/(2.0 * value3); // e gamma
  a = sqr(com1) - sqr(com);
  if (a >= 0.0) goto L4338;
  pgcnt =+ 1.0;
  goto L4335;
 L4338:
  univer = sqrt(a) * com1 * (value3 - com1)/value3; // f(m, tr) = p gamma * e gamma * e delta/e; // ::: p g + gammaEnergy + ?
  crdet(1, x, value1); // (N, N) fmax(tr) isobar sampling s.p.
  com1 = G4UniformRand();
  if (com1 > (univer / crdt[1])) goto L4335;
  massParticle[4] = massNucleon;
  massParticle[3] = com;
  // g4bertini->angid();  //:::
  pt[4] = massNucleon;
  pt[28]= massNucleon;
  //  g4bertini->alp19(); //:::
 L2040:
  return;
 L4360:
  isw[12] = 2;
  goto L4331;
 L4365:
  isw[13]= 0;
 L4366:
  i1 = -1;
  ans = (sqr(value3 - massPionCharged) - 9.0554256e27) / 4.81633308e24;
  goto L4333;
  // tr prime com1 = rlke prime
 L4370:
  com1 = (sqr(value3 + massNucleon - com) - 9.0554256e27) / 4.81633308e24;
  com2 = com;
  ans = com1;
  com4 = ftr;
  i1 = 1;
  goto L4333;
 L4375:
  com1 = (sqr(com2) - sqr(com) + value2) / (2.0 * value3); // e epsilon
  a = sqr(com1) - sqr(com2);
  if (a < 0.0) {
    pecnt += 1.0;
    goto L4380;
    // f(m1, m2, tr) = p epsilon * e epsilon * e zeta/e
  }
  univer = sqrt(a) * com1 * (value3 - com1)/value3;
  value1 = rlke - 920.0;
  crdet(1, y, value1); // (N, N) fmax(tr) isobar sampling d.p.  fmax(m1, m2, tr)
  value1 = G4UniformRand();
  if (sn > (univer * ftr/(crdt[1] * com4))) {
  L4380:
    ftr = com4;
    i1 = - 1;
    goto L4335;
  }
  value1 = G4UniformRand();
  if (value1 <= 0.5) {
    massParticle[3] = com2;
    massParticle[4] = com;
    goto L4400;
  }
  massParticle[3] = com;
  massParticle[4] = com2;
 L4400:
  //  g4bertini->angid(); //:::
  pt[16] = massNucleon;
  pt[40] = massNucleon;
  if (isw[13] != 0) {
    crdet(1, z, rlke);
    value1 = crdt[1]; // (n, p) fractional final state 3/2 l.e.
  }
  pt[2]  = 3.0;
  pt[4]  = massPionZero;
  pt[14] = 1.0;
  pt[26] = 3.0;
  pt[28] = massPionZero;
  pt[38] = 1.0;
  // g4bertini->alp28(); //:::
  goto L2040;
 L4410:
  isw[13] = 2;
  goto L4366;
}

void G4BertiniIsobarModel::rou22(G4double *v, 
			     G4double *w, 
			     G4double *x, 
			     G4double *y, 
			     G4double *z) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  i5 = static_cast<G4int>(curr[1]);
  for (G4int i = 0; i < 3; i++) {
    coordinate[i] = curr[i + 3];
    dcos[i] = curr[i + 6];
  }
  in = - 1;
  medium = static_cast<G4int>(curr[10]);
  i5 = static_cast<G4int>(curr[1]);
  geo();
  if (i1< 0) goto L5046; 
  if (curr[1] < 2.0) goto L4430; 
  if (curr[1] == 2.0) goto L4690;
  if (curr[1] < 4.0) goto L4695; 
  if (curr[1] == 4.0) goto L4696;
  iv = 1;
 L5046: 
  while(1) { // main loop
    return;
  L4430:
    isw[4] = 1; // p
  L4431:
    abz = 0.0;
    i4 = 1;
    i2 = medium;
    if (i2 < 1) goto L135;
    if (i2 > 1) goto L4470;
    i4 = 6;
    goto L4470;
  L4438:
    isw[6] = 1;
    isw[5] = 1;
    i3 = 0;
    if (com <= 3600.0) goto L4440;
  L135:
    iv = 2;
    continue;
    // com = greatest energy inside nucleus for N only
  L4440:
    if (com <= 560.0) {
      if (com <= 160.0) goto L4460; 
      else goto L4455;
    }
    // g4bertini->dfmax(); // dfmax fills out fmax(1-6) //:::
    i1 = 6;
  L4451:
    if (curr[1] >= 2.0) {
      i3 = 1;
    }
  L4452:
    // g4bertini->store(); //:::
    ex = 0.0;
    signex();
    G4int m;
    G4int lg;
    G4double a;
    G4double com;
    switch (i4) {
    case 1:  
      goto L4479;
    case 2:  
      goto L4535;
    case 3:  
      goto L4610;
    case 4:  
      goto L4630;
    case 5:  
      goto L4544;
    case 6:  
      goto L4665;
    case 7:  
      goto L4810;
    case 8: 
      goto L4825;
    case 9:  
      goto L4855;
    case 10:  
      goto L605;
    case 11:  
      goto L4955;
    }
  L605:
    iv = 3;
    continue;
    //  continue;
  L4455:
    isw[6] = 0;
    //pfmax; // pfmax fills out fmax(1-4) //:::
    i1 = 4;
    goto L4451;
  L4460:
    isw[5] = 0;
    isw[6] = 0;
    //nn; // nn fills out fmax(1-2), fmax(3-6) = 0 //:::
    i1 = 2;
    goto L4451;
  L4470:
    isw[1] = 0;
    isw[2] = 0;
    isw[3] = 0;
    m = medium + static_cast<G4int>(15.0 - 6.0 * curr[1]);
    a = curr[2] - space[m];
  L4472:
    G4int i;
    for (i = 1; i < 4; i++) {
      wkrpn[i    ] = a + space[i + 9];
      wkrpn[i + 3] = a + space[i + 3];
    }
  L4476:
    m = 4 - 3 * isw[4];
    com = wkrpn[m];
    goto L4438;
  L4479:
    switch (i2) {
    case 1:  
      goto L4890;
    case 2:  
      goto L4890;
    case 3:  
      goto L4480;
    }
  L4480:
    if (ex <= d[2]) goto L4495;
    // entry point for cascade pi's after crjab rejection, 
    // prod.reaction <  180 Fermi rejection 
    // for non-absorption reaction also all cascade nuclear rejections 
    // including crjab or rlke < production threshold
    if (d[3] != 0) goto L4660;
  L4490:
    // g4bertini->ccpes(); //:::
    if (i1 > 0 ) goto L135; // goto to punp 
    iv = 4;
    continue;
  L4495:
    if (in > 0) goto L4965;
    bg6ca(3, 0);
    ifcc = 12;
  L4496:
    medium = static_cast<G4int>(clsm);
    iv = 5;
    continue;
  L4500:
    value1 = ex;
    iv = 6;
    continue;
  L4515:
    value1 = ex + d[3];
  L500:
    iv = 7;
    continue;
    // e.p.for rejection following crjab or rlke < 180 in prod. reactions and for Fermi 
    // rejection in scattering and production reactions; for cascade pi's
  L4535:
    if (ex <= d[6]) goto L4495;
    else goto L4490; 
  L4540:
    isw[1] = 1;
  L4543:
    switch (i5) {
    case 1:  
      goto L4544;
    case 2:  
      goto L4544;
    case 3:  
      goto L4544;
    case 4:  
      goto L4825;
    }
    // e.p. for all rejections, crjab, production, and Fermi for cascade N's
  L4544:
    if (ex <= d[3]) goto L4635;
    i2 = 3;
    i4 = 2;
    i3 = 1;
    i1 = 2;
    if (d[4] == 0.0) goto L4476;
    isw[2] = 1;
    isw[3] = 1;
    i4 = 3;
    i2 = 1;
    goto L4476;
  L4610:
    switch (i5) {
    case 1:  
      goto L4611;
    case 2:  
      goto L4611;
    case 3:  
      goto L4611;
    case 4:  
      goto L4855;
    }
  L4611:
    if (ex > d[4]) goto L4625;
    bg6ca(1, 0);
    ifcc = 7;
    goto L4496;
    value1 = ex;
    iv = 8;
    continue;
  L4625:
    i4 = 4;
    i1 = 2;
    i2 = 2;
    i3 = 1;
    goto L4476;
    // e.p. for cascade N and pi0, after crjab and prod. threshold rejections
  L4630:
    if (i5 <= 2) goto L4638;
    if (i5 == 4) goto L4955;
    else goto L135;
    // e.p. for cascade N after all Fermi rejections
  L4638:
    if (ex <= d[5]) {
    L4635:
      bg6ca(2, 0);
      ifcc = 10;
      goto L4496;
    L4640:
      value1 = ex;
      goto L500;
      if (isw[3] == 0 ) goto L4543;
      else goto L4630;
    } 
    i4 = 2;
    i2 = 3;
  L4656:
    i3 = 1;
    i1 = 2;
    goto L4476;
  L4660:
    isw[1] = 1;
    if (curr[1] == 3.0) goto L4815;
    if (curr[1] > 3.0) goto L4930;
    i4 = 5; 
    i2 = 2;
    goto L4656;
    // d.p.
  L4665:
    isw[1] = 1;
    isw[2] = 1;
    isw[3] = 1;
    goto L4479;
    any = fmax[not];
    if (not < 5) goto L4684;
    if (not > 5) goto L4680;
  L1290:
    iv = 9;
    continue;
  L4680:
    if (i5 != 4) goto L4685;
    else goto L1290;
  L4684:
    if (knot >= 15) goto L5000;
    
  L4685:
    switch (i5) {
    case 1:  
      goto L4686;
    case 2:  
      goto L4686;
    case 3:  
      goto L4686;
    case 4:  
      goto L4985;
    }
  L4686:
    if (not == 2) goto L1157;
    else goto L1270;
    iv = 10;
    continue;
  L1157:
    iv = 11;
    continue;
  L1270:
    iv = 12;
    continue;
  L4690:
    isw[4] = 0; // n
    goto L4431;
  L4695:
    isw[11] = 1; // curr[1] = 3 = pi+ (11575)- pi0 = 4(15100)- pi- (14646)
  L4696:
    in = 1;
    isw[1] = 0;
    isw[2] = 0;
    isw[3] = 0;
    isw[5] = 1;
    isw[6]= 0;
    isw[7]= 0;
    isw[8]= 1;
    i6 = i5 - 2;
    i2 = medium;
    com = curr[2] - space[medium + 9];
    switch (medium) { 
    case 1:  
      isw[7] = 1;
      break;
    case 2:  
      isw[8] = 0;
      isw[7] = 1;
    case 3:  
      break;
    }
    // e.p. for pi0 after Fermi rejection
    for (i = 1; i < 4; i++) {
      wkrpn[i]     = com + space[i + 9];
      wkrpn[i + 3] = com + space[i + 3];
    }
    com = com + space[4];
    if (com > 2600.0) goto L135;
    if (com > 100.0) goto L4717;
    lg = 4;
    goto L4910;
  L4717:
    lg = 6;
    //g4bertini->spcn(); //:::
    if (value1 > 0.0) goto L4734;
    com = curr[2] - space[medium + 9];
    if (com >= 360.0) goto L4734;
    switch (i6) {
    case 1:  
      break;
    case 2:  
      goto L4915;
    case 3:  
      break;
    }
    crdet(1, v, com);
    fmax[4] = crdt[1];
  L4734:
    switch (i6) {
    case 1:  
      break;
    case 2:  
      goto L4920;
    case 3:  
      break;
    }
    i1 = 6;
  L4736:
    i4 = 7;
    i2 = 1;
    i3 = 0;
    switch (i6) {
    case 1:  
      break;
    case 2:  
      goto L4760;
    case 3:  
      break;
    }
    if (isw[11] == 0) goto L4790; // pi+
  L4760:
    if (isw[7] == 0) goto L4785; // pi-
    if (isw[8] != 0) goto L4452;
    i2 = 2;
    goto L4452;
  L4785:
    i2 = 3;
    goto L4452;
  L4790:
    i3 = 1;
    goto L4760;
  L4810:
    if (curr[1] > 2.0) goto L4479;
    a = curr[2] - space[medium + 9];
    goto L4472;
  L4815:
    i4 = 8;
  L4816:
    i2 = 2;
  L4818:
    i3 = 0;
    i1 = lg;
    if (curr[1] == 4.0) goto L4827;
    if (isw[11] == 0) i3 = 1;
    if (curr[1] < 3.0) {
    L4822:
      iv = 23;
      continue;
    }
  L4823:
    if (lg < 4) goto L4822; 
    if (lg == 4) goto L4452;
    m = 5 - abs(i5 - 4);
    univer = fmax[m];
    // g4bertini->spcn();//:::
    fmax[m] = univer;
    goto L4452;
  L4827:
    i1++;
    goto L4823;
    // e.p.for pi after crjab rejection also, after rlke <= 180 in prod. reactions 
    // and after Fermi rejection in scattering and prod. reactions
  L4825:
    if (ex > d[3]) goto L4845;
    switch (i6) {
    case 1:  
      break;
    case 2:  
      goto L4940;
    case 3:  
      break;
    }
    iv = 13;
    continue;
    any = fmax[not];
    if (curr[1] >= 3.0) goto L4832;
    ifc = 12;
  L23503:
    iv = 14;
    continue;
  L4832:
    ifcc = static_cast<G4int>((clsm - 2.0) * ((clsm * 5.5) - 8.5) + 12.05);
    goto L23503;
    if (i4 == 10) goto L4865;
    // e.p. for escape prior to choosing reactions - cascade pi+-
    if (clsm == 2.0) goto L4865;
    i4 = 2;
    goto L4816;
  L4845:
    i4 = 2; // e.p. when cascade particle escapes from region 2
    i2 = 3;
    if (d[4] != 0.0) goto L4850;
  L4846:
    switch (i6) {
    case 1:  
      goto L4818;
    case 2:  
      goto L4960;
    case 3:  
      goto L4818;
    }
  L4850:
    isw[2] = 1;
    isw[3] = 1;
    i4 = 9;
    i2 = 1;
    goto L4846;
    // e.p. after all rejections except Fermi in abs. reactions for pi+- escape 
    // from region 1 prior to choosing reaction for pi+-
  L4855:
    if (ex > d[4]) goto L4860;
    switch (i6) {
    case 1:  
      break;
    case 2:  
      goto L4970;
    case 3:  
      break;
    }
    iv = 15;
    continue;
  L4860:
    i4 = 10;
    switch (i6) {
    case 1:  
      goto L4816;
    case 2:  
      goto L4950;
    case 3:  
      goto L4816;
    }
  L4865:
    i4 = 2;
    i2 = 3;
    goto L4818;
    if (in != 0) goto L4875;
    iv = 16;
    continue;
  L4875:
    ifca = 8 * abs(i6 - 2) - 11 * (i6 - 1) * (i6 - 3);
    if (isw[1] == 0) goto L4500;
    if (isw[2] != 0) goto L4515;
    iv = 17;
    continue;
    ifca = 10 * abs(i6 - 2) + 12 * (i6 - 1) * (3 - i6);
    if (isw[3] == 0) goto L4640;
    iv = 18;
    continue;
  L4890:
    if (curr[1] < 3.0) {
      switch (medium)
	{
	case 1:
	  goto L4610;
	case 2:
	  goto L4540;
	}
    } else {
      isw[1] = 1;
    }
    if (medium == 2) goto L4825; 
    isw[2] = 1;
    isw[3] = 1;
    goto L4855;
  L4910:
    isw[5] = 0;
    G4double ans;
    switch (i6) {
    case 1:
    case 2:
      ans = bovera(curr[2], massPionCharged);
      fmax[1] = 0.20e-24 * ans; // (pi0, p) scattering
      fmax[2] = 23.0e-27 * ans; // (pi-, p) scattering
      fmax[3] = 45.1e-27 * ans; // (pi-, p) exchange scattering
      com = curr[2] - space[medium + 9]; // kinetic energy of pi outside nucleus crdet(1, v, com)
      fmax[4] = crdt[1]; // c(pi+, p) abs.
      i1 = 4;
      goto L4736;
    case 3:
      goto L4980;
    }
  L4915:
    crdet(1, w, com);
    fmax[5] = crdt[1];
  L4920:
    i1 = 7;
    goto L4736;
  L4925:
    bg6ca(3, 4);
    ifcc = 24;
  L11410:
    ka = 7;
    medium = static_cast<G4int>(clsm);
    iv = 19;
    continue;
  L4930:
    i4 = 8;
  L4931:
    i2 = 2;
    goto L4818;
  L4940:
    bg6ca(2, 3);
  L4945:
    ifcc = 21;
    goto L11410;
  L4950:
    i4 = 11;
    goto L4931;
  L4955:
    if (ex <= d[5]) goto L4940;
    else goto L4975;
  L4960:
    i1 = 9 - i5;
    goto L4818;
  L4965:
    switch (i6) // ::: not ok fix
      {
      case 1:
      case 2:
	// L231:
	iv = 20;
	continue;
      case 3:
	goto L4925;
      default:
	;
      }
  L4970:
    bg6ca(1, 2);
    goto L4945;
  L4975:
    i1 = 5;
    goto L4865;
  L4980:
    ans = bovera(curr[2], massPionZero);
    fmax[1]   = 89.2e-27 * ans; // (pi0, p) elastic scattering
    fmax[2]   = 45.1e-27 * ans; // (pi0, p) exchange scattering
    fmax[3]   = fmax[1];
    space[48] = fmax[1];
    fmax[4]   = fmax[2];
    space[49] = fmax[2];
    com       = curr[2] - space[medium + 9];
    crdet(1, w, com); // (pi-, p) abs. cross-section
    fmax[5]   = crdt[1]; // (pi0, p) abs.
    space[50] = fmax[5];
    i1 = 5;
    goto L4736;
  L4985:
    if (not < 2) goto L4990;
    if (not == 2) goto L4995;
    iv = 21;
    continue;
  L4990:
    crjab(1, x[1]); // (pi-, p) direct scattering cross-section
  L1170:
    iv = 22;
    continue;
  } // while
 L4995:
  crjab(1, y[1]); // (pi-, p) cross section  scattering cross-section
  goto L1170;
 L5000:
  if (knot >= 16) goto L4995;
  crjab(1, z[1]); // (pi-, n) drct. scattering cross-section
  goto L1170;
}

void G4BertiniIsobarModel::dcintp(G4double *w) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int i;
  G4int m;
  G4int ie;
  G4int lll;
  G4double cst;
  G4double flti;
  z[0] = rlke - 660.0;
  ie = static_cast<G4int>(z[0] / 20.0 + 5.0e-06);
  G4double univ = (z[0] - static_cast<G4double>(ie) * 20.0) / 20.0;
  ++ie;
  G4double unive = (w[ie + 1] - w[ie]) * univ + w[ie];
  // ratio 
  i2 = 0;
  G4int k = 1;
  univ = G4UniformRand();
  if (univ >= unive) goto L190;
  // lt = backward 
  // (n, p) elastic scatering, rlke > 660 
  lll = 179;
  if (w[180] > rlke) {
    i1 = -1;
    return;
  }
  for (i = 0; i < 5; ++i) {
    if (w[k + 180] > rlke) goto L40;
    k += 12;
  }
  i1 = -1;
  // error 
  return;
 L40:
  m = 12;
 L50:
  i1 = m;
  G4int l;
  for (l = 0; l < i1; ++l) {
    z[l + m] = w[lll + k];
    z[l] = w[lll + k - m];
    ++k;
  }
  univ = G4UniformRand();
  unive = (rlke - z[0]) / (z[m] - z[0]);
  i1 = m;
  for (i = 1; i < i1; ++i) {
    G4double p = (z[i + m] - z[i]) * unive + z[i];
    if (p >= univ) goto L80;
  }
  i1 = -1;
  return;
 L80:
  i1 = i;
  univ = G4UniformRand();
  flti = univ + static_cast<G4double>(i1 - 2);
  if (m <= 9) goto L230;
  switch (i1) {
  case 1:  
    i1 = -1;
    return;
  case 2:    
    cst = flti * 0.01 - 1.00;
  case 3:    
    cst = flti * 0.01 - 1.00;
  case 4:    
    cst = univ * 0.02 - 0.98;
  case 5:    
    cst = univ * 0.04 - 0.96;
  case 6:    
    cst = flti * 0.06 - 1.16;
  case 7:    
    cst = flti * 0.06 - 1.16;
  case 8:    
    cst = univ * 0.08 - 0.80;
  case 9:    
    cst = univ * 0.10 - 0.72;
  case 10:    
    cst = univ * 0.12 - 0.62;
  case 11:  
    cst = univ * 0.20 - 0.50;
  case 12:    
    cst = univ * 0.30 - 0.30;
  }
  return;
  // forward 
 L190:
  lll = 143;
  if (w[144] > rlke) {
    i1 = -1;
    return;
  }
  for (i = 0; i < 4; ++i) {
    if (w[k + 144] > rlke) goto L220;
    k += 9;
  }
  i1 = -1;
  return;
 L220:
  m = 9;
  goto L50;
 L230: 
  if (i1 == 1 ) {
    i1 = -1;
    return;
  }
  if (i1 == 2 || i1 == 3 || i1 == 4 || i1 == 5) cst = 1.0 - flti * 0.025;
  if (i1 == 6) cst = univ * 0.05 + 0.85;
  if (i1 == 7) cst = univ * 0.15 + 0.7;
  if (i1 == 8) cst = univ * 0.2 + 0.5;
  if (i1 == 9) cst = univ * 0.5;
  return;
}

void G4BertiniIsobarModel::dcpr(G4int lk) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // pdci has kk = 12, pdch has kk = 11 
  i2 = 0;
  G4int kk = lk;
  G4int k = 1;
  G4double r, cst, flti;
  G4int i;
  for (i = 0; i < 5; ++i) {
    if (w[k] >= rlke) goto L40;
    k += kk;
  }
 L20:
  i2 = 1;
  // error return 
 L30:
  return;  
 L40:
  i1 = kk;
  for (G4int l = 0; l < i1; ++l) {
    z[l + kk] = w[k];
    z[l] = w[k - kk];
    ++k;
  }
  G4double sum = 0.0; 
  G4double fract = (rlke - z[0]) / (z[kk] - z[0]);
  i1 = kk;
  for (i = 1; i < i1; ++i) {
    sum = sum + z[i] + (z[i + kk] - z[i]) * fract;
    if (r < sum) goto L70;
  }
  goto L20;
  // error 
 L70:
  r = G4UniformRand();
  i1 = i;
  if (kk > 11) goto L90;
  if (i1 > 2) goto L80;
  cst = r * 0.4;
  goto L150;
 L80:
  ++i1;
 L90:
  flti = static_cast<G4double>(i1 - 2) + r;
  switch (i1) {
  case 1:  goto L20;
  case 2:  goto L110;
  case 3:  goto L110;
  case 4:  goto L110;
  case 5:  goto L120;
  case 6:  goto L120;
  case 7:  goto L130;
  case 8:  goto L130;
  case 9:  goto L130;
  case 10:  goto L130;
  case 11:  goto L140;
  case 12:  goto L140;
  }
 L110:
  cst = flti * 0.2;
  goto L150;
 L120:
  cst = flti * 0.1 + 0.3;
  goto L150;
 L130:
  cst = flti * 0.04 + 0.6;
  goto L150;
 L140:
  cst = flti * 0.02 + 0.78;
 L150:
  if (G4UniformRand() > 0.5) goto L30;
  cst = -cst;
  goto L30;
}

void G4BertiniIsobarModel::pinst() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int k;
  i1 = 0;

  medium = static_cast<G4int>(clsm);
  if (inc != 0) goto L10;
  else goto L140;
 L10:
  inc = 0;
  i1 = medium - 1;  //:: verify
  if (i1 < 0) {
    goto L20;
  } else if (i1 == 0) {
    goto L50;
  } else {
    goto L40;
  }
 L20:
  i1 = 1;
 L30:
  return;
 L40:
  i1 = medium - 3;
  if (i1 < 0) {
    goto L60;
  } else if (i1 == 0) {
    goto L90;
  } else {
    goto L20;
  }
 L50:
  ++ipec[11];
  // incident particle on region 1 collision in region 1 
  goto L30;
 L60:
  if (d[3] != 0.0) {
    goto L80;
  }
  ++ipec[7];
  goto L30;
  // inc. particle on region 2 collision in region 2 
 L80:
  ++ipec[8];
  // inc. particle on region 1 collision in region 2 
  goto L30;
 L90:
  if (d[2] != 0.0) goto L100;
  else goto L120;
 L100:
  if (d[3] != 0.0) {
    goto L110;
  } else {
    goto L130;
  }
 L110:
  ++ipec[4];
  // inc. particle on region 1 collision in region 3 
  goto L30;
 L120:
  ++ipec[2];
  //     inc. particle on region 3 collision in region 3 
  goto L30;
 L130:
  ++ipec[3];
  //     inc. particle on region 2 collision in region 3 
  goto L30;
 L140:
  k = static_cast<G4int>(curr[10]);
  k = (medium - 1) * 3 + k;
  cc[k - 1] += 1.0;
  goto L30;
  // collision region 1 particle origin k 
}

void G4BertiniIsobarModel::rou11() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  switch (i3) {
  case 1:  goto L10;
  case 2:  goto L30;
  case 3:  goto L40;
  case 4:  goto L20;
  }
 L10:
  pt[1] = 3.0;
 L20:
  ik = it;
  pt[13] = 1.0;
 L30:
  massParticle[2] = massPionCharged;
 L40:
  if (340.0 >= rlke) goto L70;
  i3 = 1;
 L60:
  return;
 L70:
  snt = mud();
  if ((i1 = ik - 3) < 0) { //::: fix
    goto L100;
  } else if (i1 == 0) {
    goto L80;
  } else {
    goto L90;
  }
 L80:
  i3 = 2;
  goto L60;
 L90:
  i3 = 3;
  goto L60;
 L100:
  crdet(51, t, rlke);
  i3 = 4;
  goto L60;
}

void G4BertiniIsobarModel::rou15() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  switch (i3) {
  case 1:  goto L10;
  case 2:  goto L80;
  case 3:  goto L20;
  case 4:  goto L100;
  }
 L10:
  if (rlke <= 2500.0) goto L40;
 L20:
  i3 = 4;
 L30:
  goto L160;
 L40:
  if (ik < 3) {
    goto L70;
  } else if (ik == 3) {
    goto L50;
  } else {
    goto L60;
  }
 L50:
  i3 = 1;
  goto L30;
 L60:
  i3 = 2;
  goto L30;
 L70:
  snt = mud();
  crdet(51, t, rlke);
  i2 = 0;
  cst = crdt[1] - fabs(snt * (crdt[1] - crdt[0])); 
 L80:
  if (i2 <= 0) {
    goto L90;
  } else {
    goto L20;
  }
 L90:
  i3 = 3;
  goto L30;
 L100:
  if (G4UniformRand() < crdt[0]) goto L150;
  value2 = 1.0;
  // scatt. forward 
  value1 = G4UniformRand();
  if (value1 >= crdt[3]) goto L140;
 L120:
  value1 = G4UniformRand();
  // to sample from uniform dist. 
 L130:
  cst = value1 * value2;
  goto L90;
 L140:
  com = 0.0;
  value1 = big7(com, i1);
  goto L130;
 L150:
  value2 = -1.0;
  value1 = G4UniformRand();
  i1 += i2;
  if (value1 >= crdt[1]) {
    goto L140;
  } else {
    goto L120;
  }
 L160:
  return ;
}


void G4BertiniIsobarModel::rou20(G4double *t, 
			     G4double *b,  
			     G4double *rr, 
			     G4double *w,  
			     G4double *g) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  switch (i3) {
  case 1:  goto L10;
  case 2:  goto L50;
  case 3:  goto L70;
  case 4:  goto L80;
  case 5:  goto L190;
  case 6:  goto L200;
  case 7:  goto L290;
  case 8:  goto L60;
  }
 L10:
  massParticle[3] = massNucleon;
  pinst();
  if (i1 != 0) goto L20;
  else goto L40;
 L20:
  i3 = 1;
 L30:
  return;
 L40:
  i3 = 2;
  goto L30;
 L50:
  isw[8] = 0;
 L60:
  isw[9] = 2;
  i3 = 3;
  goto L30;
  // pi - single production 
 L70:
  pt[1] = 2.0;
  i3 = 4;
  goto L30;
  // collision particle pi- 
 L80:
  pt[1] = 2.0;
  pt[13] = 1.0;
 L90:
  massParticle[2] = massNucleon;
  if (740.0 >= rlke) goto L100;
  else goto L140;
 L100:
  if (300.0 >= rlke) goto L130;
  value1 = rlke - 300.0;
  crdet(5, t, value1);
  // (n-p) diff.crs. intermediate energy 
  i1 = 3;
  i2 = 3;
 L120:
  i3 = 5;
  goto L30;
 L130:
  crdet(5, b, rlke);
  // (n-p) diff.crs. low energy 
  i1 = 3;
  i2 = 1;
  goto L120;
 L140:
  if (3500.0 >= rlke) goto L150;
  else goto L20;
 L150:
  if (it >= 17) goto L20;
  dcintp(rr); // (n - p) differential cross-section high energy 
 L170:
  if (i1 >= 0) {
    goto L180;
  } else {
    goto L20;
  }
 L180:
  i3 = 6;
  goto L30;
 L190:
  pt[1] = 1.0;
  pt[13] = 2.0;
  goto L90;
 L200:
  pt[1] = 2.0;
  pt[13] = 2.0;
 L210:
  massParticle[2] = massNucleon;
  if (500.0 >= rlke) goto L220;
  else goto L230;
 L220:
  i3 = 7;
  goto L30;
 L230:
  if (1000.0 >= rlke) goto L240;
  else goto L270;
 L240:
  dcpr(12);
  if (i2 == 1) goto L20;
  // sample + mu in cst 
  snt = sqrt(1.0 - cst * cst); // (p - p) scattering 
  i3 = 8;
  goto L30; //  -, scattering backward, mu < 0 
 L270:
  if (3500.0 >= rlke) {
    goto L280;
  } else {
    goto L20;
  }
 L280:
  dcpr(11);
  if (i2 == 1) {
    goto L20;
  }
  // (p-p) differential cross-section high energy 
  goto L170;
 L290:
  pt[1] = 1.0;
  pt[13] = 1.0;
  goto L210; // no pi production possible 
}

void G4BertiniIsobarModel::isob(){
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int kase;
  G4int loop;
  G4int nwds;
  G4int j; 
  G4int k;
  G4int m;
  G4int nwrit;
  G4int ka;
  G4int ii;
  G4int ne;
  G4int lk;
  G4int jp;
  G4int  mm;
  G4int nof;
  if (nter != 0) goto L20;
  nter = 1;
  jp = 0;
  ik = 0;
  casesn = 1.0;
  ln = 2;
  sqnm = 2.2638564e27;
  rcpmv = 5.0613e10;
  G4int i;
  for (i = 0; i < 161; ++i) {
    cdd = frinn[i] * 100.0;
    ii = static_cast<G4int>(cdd + 1.0);
  }
 L20:
  loop = 0;
  inpt = 0;
  for (i = 0; i < 48; ++i) pt[i] = 0.0;
  nof = 0;
  begru = 0.0;
  nor = 0;
  //     nor=number of record 
  ++nof;
  //     einc=inc.part.en. 
  //     casesn=no.of inc.part. 
  //     andit=ang.dist. 
  //     prtin=inc.part. 
  //     trsym=struck part. 
  //     randi=input rand no. 
  //     ip=type of reaction(-1=s.p.,0=d.p.,1=choice) 
  //     ik=0,struck part.at rest--otherwise assume energy 
  //     jp=0, no print out 
  kase = 0;
  // sentinel for getting cross sections when struck particle has energy 
  for (i = 0; i < 3; ++i) calcin[i] = 0.0;
  strkp = -1.0;
  switch (no) {
  case 1:  goto L60;
  case 2:  goto L60;
  case 3:  goto L70;
  case 4:  goto L80;
  case 5:  goto L70;
  }
 L60:
  massParticle[0] = massNucleon;
  goto L120;
 L70:
  massParticle[0] = massPionCharged;
  goto L90;
 L80:
  massParticle[0] = massPionZero;
 L90:
  if (no >= 4) goto L110;
  isw[10] = 1;
  goto L120;
 L110:
  isw[10] = 0;
 L120:
  massParticle[1] = massNucleon;
  ka = 1;
  energy[0] = einc * rcpmv + massParticle[0];
  pxyz[0] = 0.0;
  pxyz[4] = 0.0;
  pxyz[8] = sqrt(sqr(energy[0]) - sqr(massParticle[0]));
  rlke = einc;
  if (ik != 0) goto L150;
  p2 = 0.0;
  pxyz[1] = 0.0;
  pxyz[5] = 0.0;
  pxyz[9] = 0.0;
  energy[1] = massNucleon;
  goto L150;
 L140:
 L150:
  if (no >= 3) goto L170;
  if (rlke <= 3500.0) goto L240;
  else goto L200;
 L170:
  if (rlke <= 180.0) goto L200;
  if (rlke <= 2500.0) goto L190;
  else goto L200;
 L190:
  if (kase != 0) goto L980;
  else goto L240;
 L200:
  if (ik != 0) goto L210;
  else goto L220;
 L210:
  if (kase != 0) goto L230;
 L220:
  ne = 3;
  // rlke in error 
  goto L1060;
 L230:
  ++loop;
  if (loop <= 1000) goto L140;
  else goto L220;
 L240:
  if (ip >= 0) goto L470;
  if (no >= 3) goto L380;
  else goto L470;
 L260:
  lk = 2;
  i3 = 1;
 L270:
  qrdet(1, pspcl, calcin[0]);
  // argument = pspcl, (p, p), (n, n), s.p. 
  calcin[0] = crdt[0];
  if (lk >= 6) goto L330;
 L280:
  qrdet(1, pec, rlke); // elastic cross section 
 L290:
  calcin[2] = crdt[0] + calcin[0] + calcin[1];
  // total 
  // (p, p), (n, n), elastic 
  goto L640;
 L300:
  lk = 3;
  i3 = 2;
 L310:
  qrdet(1, spcln, calcin[0]);
  // argument = nspcl, (n, p), (p, n), s.p. 
  calcin[0] = crdt[0];
  if (lk >= 6) goto L330;
 L320:
  qrdet(1, ecn, rlke);
  goto L290;
  // (p, n), (n, p), elastic 
 L330:
  cratio = calcin[0] / (calcin[0] + calcin[1]);
  value1 = G4UniformRand();
  if (value1 <= cratio) {
    goto L340;
  } else {
    goto L360;
  }
 L340:
  if (ka != no) goto L370;
  i3 = 1;
  goto L280;
 L360:
  if (ka != no) {
    goto L320;
  } else {
    goto L280;
  }
 L370:
  i3 = 2;
  goto L320;
 L380:
  isw[9] = 0;
  isw[8] = 0;
  // 9 set for pi0, 10 set for (G4int pi+, n), (pi-, p), (pi0, p), and (pi0, n) 
  calcin[0] = rlke - 180.0;
  lk = 1;
  // pi 
  i3 = 0;
  if (no < 4) {
    goto L390;
  } else if (no == 4) {
    goto L430;
  } else {
    goto L460;
  }
 L390:
  switch (ka) {
  case 1:  goto L400;
  case 2:  goto L410;
  }
 L400:
  qrdet(1, ppscl, calcin[0]);
  calcin[0] = crdt[0];
  qrdet(1, ppec, rlke);
  goto L290;
  // (pi+, p), s.p. or (pi-, n) 
 L410:
  qrdet(1, pmscl, calcin[0]);
  calcin[0] = crdt[0];
  qrdet(1, pmec, rlke);
 L420:
  calcin[2] = crdt[0];
  isw[9] = 2;
  qrdet(1, pmxc, rlke);
  calcin[2] += crdt[0] + calcin[0] + calcin[1];
  // (pi+, n), s.p or (pi-, p) 
  goto L640;
 L430:
  isw[8] = 2;
  switch (ka) {
  case 1:  goto L440;
  case 2:  goto L450;
  }
 L440:
  qrdet(1, pnscl, calcin[0]); //::: check paremeters (all)
  calcin[0] = crdt[0];
  qrdet(1, pnec, rlke);
  goto L420;
  // (pi0, p), s.p. 
 L450:
  qrdet(1, pnnsl, calcin[0]);
  calcin[0] = crdt[0];
  qrdet(1, pnnec, rlke);
  goto L420;
  // (pi0, n), s.p. 
 L460:
  switch (ka) {
  case 1:  goto L410;
  case 2:  goto L400;
  }
  // (pi-, p), (pi-, n), s.p. 
 L470:
  if (rlke <= 920.0) {
    goto L480;
  } else {
    goto L540;
  }
 L480:
  if (ip != 0) {
    goto L490;
  } else {
    goto L200;
  }
 L490:
  if (rlke <= 360.0) goto L200;
  if (kase != 0) goto L980;
  calcin[0] = rlke - 360.0;
  switch (no) {
  case 1:  goto L520;
  case 2:  goto L530;
  }
 L520:
  isw[3] = 1;
  switch (ka) {
  case 1:  goto L580;
  case 2:  goto L610;
  }
 L530:
  isw[3] = 0;
  switch (ka) {
  case 1:  goto L610;
  case 2:  goto L580;
  }
 L540:
  if (kase != 0) goto L980;
  calcin[1] = rlke - 920.0;
  calcin[0] = rlke - 360.0;
  switch (no) {
  case 1:  goto L560;
  case 2:  goto L630;
  }
 L560:
  isw[3] = 1;
  switch (ka) {
  case 1:  goto L570;
  case 2:  goto L600;
  }
 L570:
  qrdet(1, pdpcl, calcin[1]);
  // argument = pdpcl, (p, p), (n, n), d.p. 
  calcin[1] = crdt[0];
  // (n, n) or (p, p) d.p. 
 L580:
  if (ip >= 0) {
    goto L590;
  } else {
    goto L260;
  }
 L590:
  lk = 2 * static_cast<G4int>(ip) + 4;
  i3 = 3;
  goto L270;
 L600:
  qrdet(1, dpcln, calcin[1]);
  // argument = ndpcl, (n, p), (p, n), d.p. 
  calcin[1] = crdt[0];
 L610:
  if (ip >= 0) {
    goto L620;
  } else {
    goto L300;
  }
 L620:
  lk = static_cast<G4int>(ip) + 5;
  i3 = 4;
  goto L310;
  // (n, p) or (p, n) d.p. 
 L630:
  isw[3] = 0;
  switch (ka) {
  case 1:  goto L600;
  case 2:  goto L570;
  }
 L640:
  if (jp < 0) {
    goto L1060;
  } else if (jp == 0) {
    goto L660;
  }
  // jp = 0, no print out 
  // write(io,425)einc,no,ka,casesn,(rands(i),i=2,4),(calcin(i),i=1,3) 
  // 1, ip,nof,ik 
  // 425 format('   einc      no       ka     casesn     randi      cal 
  // 1cin ip nof ik'/1h ,d11.3,2i8,d11.3,1x,3z4, 
  // 2 3d11.3,3i8) 
 L660:
  kase = 1;
  if (ik != 0) goto L140;
  switch (no) {
  case 1:  goto L1000;
  case 2:  goto L1000;
  case 3:  goto L690;
  case 4:  goto L690;
  case 5:  goto L690;
  }
 L680:
  i3 = 0;
 L690:
  qou17(fripn, pnmi, fmxsp, pcfsl, pnfsl);
  // 
  // ::: pattern if(i3)3050,2084,2081 huomaa viimeinen linkki seuraavalle riville
  // 2081 if(col(15)-1.0)2084,4035,2082
  if (i3 < 0) {
    goto L1100;
  } else if (i3 == 0) {
    goto L740;
  }
 L700:
  if (col[14] < 1.0) { 
    goto L740;
  } else  if (col[14] == 1.0) {
    goto L800;
  }
  if (col[14] < 3.0) {
    goto L780;
  } else if (col[14] == 3.0) {
    goto L770;
  }
  if (col[14] < 5.0) {
    goto L790;
  } else if (col[14] == 5.0) {
    goto L1030;
  }
  ne = 11;
  goto L1060;
 L740:
  qollm();
  if (pt[37] == 0.0) {
    i3 = 1;
    goto L810;
  }
  i3 = 2;
  goto L810;
 L770:
  i3 = 4;
  goto L810;
 L780:
  i3 = 5;
  goto L810;
 L790:
  i3 = 6;
  goto L810;
 L800:
  i3 = 3;
 L810:
  qou18();
  i3 = i3;
  switch (i3) {
  case 1:  goto L700;
  case 2:  goto L820;
  case 3:  goto L700;
  case 4:  goto L830;
  }
 L820:
  ne = 12;
  goto L1060;
 L830:
  massParticle[3] = massNucleon;
  k = static_cast<G4int>(ip);
  j = 2;
  if (ip < 0) {
    goto L840;
  } else if (ip == 0) {
    goto L870;
  } else {
    goto L850;
  }
 L840:
  nwrit = 22;
  mm = 26;
  goto L880;
 L850:
  if (pt[37] != 0.0) {
    goto L860;
  } else {
    goto L840;
  }
 L860:
  ++k;
 L870:
  nwrit = 29;
  mm = 38;
 L880:
  out[0] = k;
  i1 = mm;
  for (i = 1; i < i1; i += 12) {
    m = i;
    k = m + 3;
    for (G4int n = 0; n < 2; ++n) {
      i2 = k;
      for (G4int l = m; l < i2; ++l) {
	out[j] = pt[l];
	++j;
      }
      m = i + 6;
      k = m + 2;
    }
  }
  ++nor;
  if (jp < 0) {
    goto L1060;
  } else if (jp == 0) {
    goto L940;
  }
  //     jp=0, no print out 
  //:::  dolio(c3, c1, (char *)nor, (ftnlen)sizeof(G4int));
  i1 = nwrit;
  for (i = 0; i < i1; ++i) {
    //:::	dolio(c5, c1, (char *)out[i - 1], (ftnlen)sizeof());
  }
  //:::  dolio(c5, c1, (char *)e[1], (ftnlen)sizeof());
  for (i = 1; i < 10; i += 4) {
    //:::	dolio(c5, c1, (char *)pxyz[i - 1], (ftnlen)sizeof());
  }
  // 5027 format(i7/(7d14.5)) 
 L940:
  nwds = nwrit + 4;
  // number of words in record 
  for (i = 0; i < 48; ++i) pt[i] = 0.0;
  pgcnt = 0.0;
  pacnt = 0.0;
  loop = 0;
  begru += 1.0;
  if (casesn < begru) {
    goto L1060;
  } else if (casesn == begru) {
    goto L1090;
  }
  if (ik != 0) goto L140;
 L980:
  switch (lk) {
  case 1:  goto L680;
  case 2:  goto L990;
  case 3:  goto L1010;
  case 4:  goto L1020;
  case 5:  goto L1050;
  case 6:  goto L1110;
  }
 L990:
  i3 = 1;
 L1000:
  qou21(frinn, dmin, fmxsn, fmxdn, fsln);
  if (icoun1 > icoumx || icoun2 > icoumx) {
    nopart = -1;
    return;
  }
  goto L690;
 L1010:
  i3 = 2;
  goto L1000;
 L1020:
  i3 = 3;
  goto L1000;
 L1030:
  qou19();
  if (i3 >= 0) goto L830;
  ne = 13;
  goto L1060;
 L1050:
  i3 = 4;
  goto L1000;
 L1060:
  //     dolio(c3, c1, (char *)ne, (ftnlen)sizeof(G4int));
  //      dolio(c3, c1, (char *)nor, (ftnlen)sizeof(G4int));
  //      dolio(c3, c1, (char *)no, (ftnlen)sizeof(G4int));
  //      dolio(c3, c1, (char *)nof, (ftnlen)sizeof(G4int));
  //      dolio(c5, c1, (char *)begru, (ftnlen)sizeof());
  //      dolio(c5, c1, (char *)strkp, (ftnlen)sizeof());
  //      dolio(c5, c1, (char *)rlke, (ftnlen)sizeof());
  for (i = 0; i < 4; ++i) {
    // dolio(c3, c1, (char *)rands[i - 1], (ftnlen)sizeof(G4int))
    ;
  }
  // 136 format(' error at 135'/1h ,'  ne     nor    no     nof 
  // 1 begru         strkp          rlke          rands(1-4)'/4i7,3d13.42,1x,4z4) 
  if (ne  >= 3) {
    goto L1070;
  } else {
    goto L1090;
  }
 L1070:
  if (calcin[0] != 0.0) {
    goto L1080;
  } else {
    goto L1090;
  }
 L1080:
  ++nor;
 L1090:
  return;
 L1100:
  ne = 10;
  goto L1060;
 L1110:
  value1 = G4UniformRand();
  if (value1 <= cratio) {
    goto L1120;
  } else {
    goto L1130;
  }
 L1120:
  if (ka != no) {
    goto L1010;
  } else {
    goto L990;
  }
 L1130:
  if (ka != no) {
    goto L1140;
  } else {
    goto L1020;
  }
 L1140:
  goto L1050;
}

void G4BertiniIsobarModel::qlp19() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double univ = G4UniformRand();
  pt[1] = 1.0;
  pt[25] = 1.0;
  pt[13] = 3.0;
  pt[15] = massPionZero;
  if (isw[11] <= 0) {
    goto L10;
  } else {
    goto L140;
  }
 L10:
  if (univ <= 0.25) {
    goto L20;
  } else {
    goto L110;
  }
 L20:
  if (isw[3] != 0) {
    goto L40;
  }
  pt[1] = 2.0;
 L40:
  univ = G4UniformRand();
  if (univ <= twoThirds) {
    goto L50;
  } else {
    goto L90;
  }
 L50:
  pt[13] = 4.0;
 L60:
  if (isw[3] != 0) {
    goto L80;
  }
 L70:
  pt[25] = 2.0;
 L80:
  return ;
 L90:
  pt[15] = massPionCharged;
  if (isw[3] != 0) goto L70;
 L100:
  pt[13] = 5.0;
  goto L80;
 L110:
  if (isw[3] != 0) goto L130;
  pt[25] = 2.0;
  goto L90;
 L130:
  pt[1] = 2.0;
  pt[15] = massPionCharged;
  goto L80;
 L140:
  if (univ <= 0.5) {
    goto L150;
  } else {
    goto L190;
  }
 L150:
  if (isw[3] != 0) {
    goto L160;
  } else {
    goto L170;
  }
 L160:
  pt[1] = 2.0;
 L170:
  univ = G4UniformRand();
  if (univ <= oneThird) goto L90;
  pt[13] = 4.0;
  goto L60;
 L190:
  if (isw[3] != 0) goto L210;
  pt[1] = 2.0;
 L210:
  univ = G4UniformRand();
  if (univ <= twoThirds) {
    goto L220;
  } else {
    goto L230;
  }
 L220:
  pt[13] = 4.0;
  if (isw[3] != 0) {
    goto L70;
  } else {
    goto L80;
  }
 L230:
  pt[15] = massPionCharged;
  if (isw[3] != 0) {
    goto L240;
  } else {
    goto L70;
  }
 L240:
  goto L100;
}

void G4BertiniIsobarModel::qlp28() { 
  // decription: 
  // parameters:
  // uses:
  // changes:
  if (isw[12] != 0) goto L230;
  r = G4UniformRand();
  if (r <= 0.6) {
    goto L20;
  } else {
    goto L120;
  }
 L20:
  pt[3] = massPionCharged;
  r = G4UniformRand();
  if (isw[3] != 0) {
    goto L30;
  } else {
    goto L90;
  }
 L30:
  if (r <= oneThird) {
    goto L40;
  } else {
    goto L70;
  }
 L40:
  pt[25] = 5.0;
 L50:
  pt[27] = massPionCharged;
 L60:
  return ;
 L70:
  pt[25] = 4.0;
 L80:
  pt[37] = 2.0;
  goto L60;
 L90:
  pt[1] = 5.0;
  pt[13] = 2.0;
  if (r <= oneThird) {
    goto L100;
  } else {
    goto L110;
  }
 L100:
  pt[27] = massPionCharged;
  goto L80;
 L110:
  pt[25] = 4.0;
  goto L60;
 L120:
  r = G4UniformRand();
  if (isw[3] != 0) {
    goto L130;
  } else {
    goto L180;
  }
 L130:
  if (r <= twoThirds) {
    goto L140;
  } else {
    goto L160;
  }
 L140:
  pt[1] = 4.0;
 L150:
  r = G4UniformRand();
  if (r <= twoThirds) {
    goto L110;
  } else {
    goto L100;
  }
 L160:
  pt[13] = 2.0;
 L170:
  pt[3] = massPionCharged;
  goto L150;
 L180:
  if (r <= twoThirds) {
    goto L190;
  } else {
    goto L220;
  }
 L190:
  pt[1] = 4.0;
 L200:
  pt[13] = 2.0;
 L210:
  r = G4UniformRand();
  if (r <= twoThirds) {
    goto L70;
  } else {
    goto L40;
  }
 L220:
  pt[1] = 5.0;
  pt[3] = massPionCharged;
  goto L210;
 L230:
  if (r <= value1) {
    goto L240;
  } else {
    goto L270;
  }
 L240:
  pt[3] = massPionCharged;
  if (isw[3] != 0) {
    goto L260;
  }
  pt[1] = 5.0;
  pt[13] = 2.0;
  goto L50;
 L260:
  pt[37] = 2.0;
  goto L40;
 L270:
  r = G4UniformRand();
  if (isw[3] != 0) {
    goto L280;
  } else {
    goto L310;
  }
 L280:
  if (r <= oneThird) { //::: tarkista
    goto L290;
  } else {
    goto L300;
  }
 L290:
  pt[3] = massPionCharged;
  goto L200;
 L300:
  pt[1] = 4.0;
  goto L210;
 L310:
  if (r <= oneThird) {
    goto L320;
  } else {
    goto L330;
  }
 L320:
  pt[1] = 5.0;
  goto L170;
 L330:
  pt[13] = 2.0;
  goto L140;
}

void G4BertiniIsobarModel::qlpha() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double univ = G4UniformRand();
  if (value3 < 0.0) {
    goto L300;
  } else if (value3 == 0) {
    goto L10;
  } else {
    goto L140;
  }
 L10:
  if (univ <= value1) {
    goto L20;
  } else {
    goto L120;
  }
 L20:
  if (isw[10] != 0) goto L40;
  pt[1] = 5.0;
  pt[25] = 2.0;
 L40:
  pt[3] = massPionCharged;
  massParticle[3] = massPionCharged;
  univ = G4UniformRand();
  if (univ <= value2) {
    goto L50;
  } else {
    goto L70;
  }
 L50:
  pt[13] = 4.0;
 L60:
  return ;
 L70:
  if (isw[10] != 0) goto L110;
  pt[25] = 1.0;
 L90:
  pt[13] = 5.0;
 L100:
  pt[15] = massPionCharged;
  goto L60;
 L110:
  pt[25] = 2.0;
  goto L100;
 L120:
  pt[1] = 4.0;
  if (isw[10] != 0) goto L100;
  pt[13] = 5.0;
  goto L110;
 L140:
  if (univ <= value1) {
    goto L150;
  } else {
    goto L200;
  }
 L150:
  massParticle[3] = massPionCharged;
  if (isw[10] != 0) {
    goto L160;
  } else {
    goto L190;
  }
 L160:
  pt[1] = 5.0;
 L170:
  pt[15] = massPionCharged;
  pt[3] = massPionCharged;
  goto L60;
 L190:
  pt[13] = 5.0;
  pt[25] = 2.0;
  goto L170;
 L200:
  if (univ <= value2) {
    goto L210;
  } else {
    goto L250;
  }
 L210:
  pt[1] = 4.0;
  univ = G4UniformRand();
  if (univ <= value3) goto L240;
  if (isw[10] != 0) goto L50;
 L230:
  pt[25] = 2.0;
  goto L50;
 L240:
  if (isw[10] != 0) {
    goto L110;
  } else {
    goto L90;
  }
 L250:
  massParticle[3] = massPionCharged;
  pt[3] = massPionCharged;
  univ = G4UniformRand();
  if (univ <= twoThirds) {
    goto L260;
  } else {
    goto L280;
  }
 L260:
  if (isw[10] != 0) goto L230;
  pt[1] = 5.0;
  goto L50;
 L280:
  if (isw[10] != 0) goto L90;
  pt[25] = 2.0;
  goto L160;
 L300:
  if (univ <= value1) {
    goto L310;
  } else {
    goto L340;
  }
 L310:
  massParticle[3] = massPionCharged;
  pt[3] = massPionCharged;
  univ = G4UniformRand();
  if (value3 != -1.0) goto L330;
  if (univ <= oneThird) {
    goto L90;
  } else {
    goto L230;
  }
 L330:
  pt[1] = 5.0;
  if (univ <= oneThird) {
    goto L110;
  } else {
    goto L50;
  }
 L340:
  if (univ <= value2) {
    goto L350;
  } else {
    goto L380;
  }
 L350:
  pt[1] = 4.0;
  univ = G4UniformRand();
  if (value3 != -1.0) goto L370;
  if (univ <= twoThirds) {
    goto L50;
  } else {
    goto L110;
  }
 L370:
  if (univ <= twoThirds) {
    goto L230;
  } else {
    goto L90;
  }
 L380:
  massParticle[3] = massPionCharged;
  pt[3] = massPionCharged;
  if (value3 + 1.0 != 0.0) {
    goto L390;
  } else {
    goto L160;
  }
 L390:
  pt[13] = 5.0;
  goto L110;
} 

void G4BertiniIsobarModel::qngid() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // calculates cos and sin theta, sin and cos phi
  G4int icurr;
  icurr = static_cast<G4int>(curr[0]);
  // curr(1) has not been set;therefore,icurr=0 and one goes only to 1 
  // this is equivalend to permitting only single production. 
  // if(icurr.gt.0)write(io,1)icurr,curr(1),it 
  // 1 format(1h ,2x,'qngid:icurr=',i5,2x,'curr(1)=',d12.6,2x,'it=',i5) 

  switch (icurr) {
  case 1:  goto L10;
  case 2:  goto L10;
  case 3:  goto L30;
  case 4:  goto L30;
  case 5:  goto L30;
  }
  // incident particle - N 
 L10:
  if (it == 21 || it == 22) {
    goto L20;
  }
  // single production 
  if (rlke > 3500.0) G4cerr << "qngid1" << G4endl;
  if (rlke < 500.0) goto L70;
  tesiso = 0.75;
  if (rlke < 1000.0) goto L50;
  tesiso = 0.5;
  if (rlke < 1300.0) goto L50;
  tesiso = 0.25;
  if (rlke < 2500.0) goto L50;
  goto L60;
  // double production 
 L20:
  if (rlke > 3500.0) G4cerr << "qngid2" << G4endl;
  goto L60;
  // incident particle - pi 
 L30:
  r = G4UniformRand();
  if (rlke > 2500.0) G4cerr << "qngid3" << G4endl;
  cst = -0.9999995;
  snt = 0.003162;
  if (it != 11) goto L40;
  if (r <= 0.75) goto L70;
  goto L80;
  // (pi+, p), (pi-, n) 
  // (pi0, n), (pi0, p) 
 L40:
  if (it != 12 && it != 28) G4cerr << "qngid4" << G4endl;
  if (rlke < 500.0) cst = -cst;
  if (r <= 0.8) {
    goto L70;
  }
  goto L80;
 L50:
  r = G4UniformRand();
  if (r <= tesiso) goto L70;
  // backward/forward 
 L60:
  r = G4UniformRand();
  // test for direction 
  cst = 0.9999995;
  snt = 0.003162;
  if (r <= 0.5) goto L80;
  cst = -0.9999995;
  goto L80;
  // isotropic 
 L70:
  pol1(cst, snt);
  // calculates cos, sin phi 
 L80:
  azio(sopc, sops);
  return;
}

void G4BertiniIsobarModel::qoll() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double a = sqr(massParticle[3]);
  col[14] = 0.0;
  col[0] = energy[0] + energy[1];
  // total energy particles 1 and 2 
  for (G4int i = 0; i < 3; ++i) col[i] = sqr(massParticle[i]);  
  col[4] = col[2] + col[1] + 
    (energy[0] * energy[1] - 
     (pxyz[0] * pxyz[1] + pxyz[4] * pxyz[5] + pxyz[8] * pxyz[9])
     ) * 2.0;
  col[5] = sqrt(col[4]);
  col[6] = col[5] / col[0];
  // gam 
  col[7] = col[5] * 2.0;
  col[8] = (col[3] + col[4] - a) / col[7];
  com2 = col[8] * col[8];

  if (col[3] <= 2.9882156e27f) {
    goto L30;
  } else {
    goto L50;
  }
  // >, pm(3) = isobar--lte,test for roundoff range,(min)sqd+or-5e23 
 L30:
  if (col[3] >= 2.9872156e27f) {
    goto L40;
  } else {
    goto L60;
  }
  // < , pi or N mass = pm(3) 
 L40:
  col[3] = 2.9877156e27;
  massParticle[2] = 5.466005e13;
 L50:
  if (com2 >= col[3]) {
    goto L90;
  } else {
    goto L70;
  }
 L60:
  if (col[3] <= sqnm) {
    goto L50;
  } else {
    goto L120;
  }
  // <=, have N or p--gt, go to error 
 L70:
  if (com2 >= col[3] * 0.99) {
    goto L80;
  } else {
    goto L90;
  }
 L80:
  com2 = col[3];
  col[8] = massParticle[2];
 L90:
  col[9] = sqrt(com2 - col[3]);
  // p3 prime 
  if (ik != 0) goto L100;
  else goto L170;
 L100:
  col[10] = (col[4] + col[1] - col[2]) / col[7]; // e1 prime 
  col[11] = sqrt(sqr(col[10]) - col[1]); // p1 prime 
  col[12] = (col[6] * energy[0] - col[10]) / col[11]; // beta 
  com = 1.0 - (sqr(col[12]) + sqr(col[6]));
  if (com >= 5e-6f) goto L150;
  if (com >= 5e-6f) goto L140;
 L120:
  col[14] = 1.0; // error 
 L130:
  return;
 L140:
  col[13] = 0.002236067977;
  goto L160;
 L150:
  col[13] = sqrt(com); // alpha 
 L160:
  energy[2] = (col[8] + col[9] * (col[12] * cst + col[13] * sopc * snt)) / col[6];
  energy[3] = col[0] - energy[2];
  goto L130;
 L170:
  for (i = 10; i < 14; ++i) col[i] = 0.0;
  energy[2] = 0.0;
  energy[3] = 0.0;
  goto L130;
  // for part at rest e1 and p1 prime,alpha and beta not 
  // calculated. e(3) ande(4) done later 
}

void G4BertiniIsobarModel::qene(G4double *z) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int i;
  G4int i101;
  G4double bz;
  G4double cz;
  G4double xz;
  G4double disc;
  G4double ra;
  G4double rb;
  G4double rc;
  G4double sqba;
  G4double sqac;
  G4double sqa;
  G4double sba;
  G4double sca;
  G4double plus;
  G4double aminus;
  G4double cd = com * 100.0;
  i = static_cast<G4int>(cd + 1.0);
  G4double az = z[i];
  if (i == 1) goto L150;
  bz = z[i + 1];
  i101 = 101 - (i + 1);
  // if(i101.lt.0)write(io,3)i101,bz 
  // 3 format(1h ,2x,'qene:i101=',i5,2x,'bz=',d15.9/) 
  // if ((i1 = 101 - (i + 1)) < 0) {  ::: ks. siisititty versio
  if (i > 100) {
    goto L70;
  } else if (i == 0) {
    goto L20;
  } else {
    goto L30;
  }
 L20:
  cz = bz + (bz - az) * 0.5;
  goto L40;
 L30:
  cz = z[i + 2];
 L40:
  xz = cd - static_cast<G4double>(i - 1);
  sca = cz - az; // f(2) - f(0) 
 L50:
  sba = bz - az; // f(1) - f(0) 
  sqa = sqr(az); // sqr(f(0)) 
  sqac = sqa - cz * cz; // sqr(f(0))-sqr(f(2)) 
  sqba = bz * bz - sqa; // sqr(f(1))-sqr(f(0)) 
  rb = sqac + 2 * sqba; // (asq-csq) + 2 * (bsq-asq) 
  rc = az * cz * sca - sba * (az * 2.0 * bz + xz * (bz - cz) * sca);
  ra = sca - sba - sba; // (c - a) - 2 * (b - a) 
  if (ra != 0.0) goto L60;
  com = az + xz * sba;
  goto L80;
 L60:
  disc = sqr(rb) - 4.0 * ra * rc;
  if (disc >= 0.0) {
    goto L90;
  } else {
    goto L70;
  }
  // sqr(b) - 4 * a * c 
 L70:
  G4cerr << "qene1" << G4endl;
 L80:
  return;
 L90:
  disc = sqrt(disc);
  plus = (disc - rb) / (ra + ra);
  aminus = (- rb - disc) / (ra + ra);
  if (i == 1) goto L160;
 L100:
  if (plus > bz) goto L120;
  if (plus < az) goto L120;
  if (aminus > bz) goto L110;
  if (aminus >= az) goto L140;
 L110:
  com = plus;
  goto L80;
 L120:
  if (aminus > bz) goto L70;
  if (aminus < az) goto L70;
 L130:
  com = aminus;
  goto L80;
 L140:
  ra = xz * sba + az;
  rb = fabs(ra - aminus);
  rc = fabs(ra - plus);
  if (rb > rc) goto L110;
  goto L130;
 L150:
  cz = z[i + 1];
  sca = cz - az;
  bz = az + sca * 0.7071067812;
  xz = cd + cd;
  goto L50;
  // (cz - az) (cz-az) = c,cz = mass for r = 1, az = mass for r = 0, 
  // c = const.0 or parabola 
  // (m - az) (m - az) = 0.5 * c, determines mass, bz, for r = 1/2 
 L160:
  bz = cz;
  xz -= cd;
  sba = cz - az;
  goto L100;
}

void G4BertiniIsobarModel::qdk() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double univ = sqr(pnidk[5]); // m(p1) squared decay pion mass squared 
  pnidk[6] = (pnidk[0] * pnidk[0] + univ - sqnm) / (pnidk[0] * 2.0); // e(pi)prime  decay pion energy prime 
  pnidk[7] = sqrt(pnidk[6] * pnidk[6] - univ); // decay pi momentum prime  p(d) 
  pol1(pnidk[19], pnidk[20]); // cos theta, sin theta 
  azio(pnidk[21], pnidk[22]); // cos phi, sin phi 
  pnidk[8] = pnidk[21] * pnidk[20] * pnidk[7]; // decay pi x momentum component prime 
  pnidk[9] = pnidk[20] * pnidk[22] * pnidk[7]; // p(p1) prime y 
  pnidk[10] = pnidk[19] * pnidk[7];  // p(p1) prime z 
  univ = pnidk[8] * pnidk[1] + pnidk[9] * pnidk[2] + pnidk[10] * pnidk[3];  // p p1 prime dot p 
  pnidk[11] = (pnidk[6] * pnidk[4] + univ) / pnidk[0]; // decay pi energy  e(pi) 
  pnidk[12] = pnidk[4] - pnidk[11];
  univ = (pnidk[4] / pnidk[0] - 1.0) * univ / 
    (pnidk[1] * pnidk[1] + pnidk[2] * pnidk[2] + pnidk[3] * pnidk[3]); // (e/m - 1.0) * p(p1)prime dot p/p squared 
  G4double unive = pnidk[6] / pnidk[0]; // e pi prime over mass 
  for (G4int i = 1; i < 4; ++i) {
    pnidk[i + 12] = pnidk[i] * (univ + unive) + pnidk[i + 7];
    pnidk[i + 14] = pnidk[i] - pnidk[i + 12];
  }
  return;
  // pion momentum components and nucleon momentum 
  // components 
} 

void G4BertiniIsobarModel::qstor() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int  j;
  G4int k;
  G4int l;
  G4int m;
  G4int jj;
  G4double univ;

  l = i1 * 12 - 28;
  if (i2 < 0) {
    goto L10;
  } else if (i2 == 0) {
    goto L60;
  } else {
    goto L70;
  }
 L10:
  jj = 0;
  if (massParticle[2] <= massNucleon) goto L30;
  ++i1;
  jj = 1;

 L30:
  // x-y-z-coordinates of collision point 
  univ = sqrt(pxyz[i1 - 1] * pxyz[i1 - 1] + pxyz[i1 + 3] * pxyz[i1 + 3] + pxyz[i1 + 7] * pxyz[i1 + 7]);
  k = i1 + 7;
  i1 = k;
  G4int i;
  for (i = i1; i < i1; i += 4) {
    pt[l] = pxyz[i] / univ;
    ++l;
  }
  i1 -= jj;
 L50:
  return;
 L60:
  k = 14;
  goto L90;
 L70:
  if (i2 >= 2) goto L110;
  k = 17;
 L90:
  univ = sqrt(sqr(pnidk[k - 1]) + sqr(pnidk[k]) + sqr(pnidk[k + 1]));
  pt[l - 4] = 1.0;
  j = k + 2;
  i1 = j;
  for (i = k; i < i1; ++i) {
    pt[l] = pnidk[i] / univ;
    ++l;
  }
  goto L50;
 L110:
  univ = sqrt(sqr(pt[l - 4]) + sqr(pt[l - 3]) + sqr(pt[l - 2]));
  k = l - 1;
  m = l - 2;
  i1 = k;
  for (i = m; i < i1; ++i) {
    pt[l] = pt[i] / univ;
    ++l;
  }
  pt[m - 1] = 1.0;
  goto L50;
} 

void G4BertiniIsobarModel::qollm() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int k;
  G4double univ;
  G4double unive;
  G4double a;
  G4double b;
  G4double z;
  // description:
  if (ik != 0) goto L10;
  else goto L60;
 L10:
  univ = energy[1] + col[5] - col[10];
  unive = energy[0] + col[10];
  univer = col[0] + col[5];
  k = 16;
  G4int i;
  for (i = 0; i < 9; i += 4) {
    col[k] = (pxyz[i] * univ - pxyz[i] * unive) / univer;
    col[k + 3] = (pxyz[i] + pxyz[i]) / col[0];
    // vx 
    ++k;
  }
  col[21] = -pxyz[8] * pxyz[5] / col[0]; // qx, abbreviated form since pxyz(1) = pxyz(5) = 0.0 
  col[22] = pxyz[1] * pxyz[8] / col[0]; // qy, abbreviated form since pxyz(1) = pxyz(5) = 0.0 
  col[23] = 0.0; // abbreviated form since pxyz(1) = pxyz(5) = 0.0 
  a = snt / col[13];
  b = a * col[9];
  // (-beta * cos phi * sin theta / alpha + cos theta) / p1p * p3p 
  univ = col[9] * (cst - a * sopc * col[12]) / col[11];
  unive = b * sops / col[11];
  // p3p * sin phi * sin theta / p1p * alpha 
  univer = sopc * b + (energy[2] + col[8]) / (col[6] + 1.0);
  // cos phi * sin theta * p3p / alpha  +  (e3 + e3p) / (1.0 + gamma) 
  k = 19;
  for (i = 2; i < 11; i += 4) {
    pxyz[i] = col[k] * univer + col[k + 3] * unive + col[k - 3] * univ;
    ++k;
  }
  for (i = 0; i < 9; i += 4) pxyz[i + 3] = pxyz[i] + pxyz[i + 1 ] - pxyz[i + 2];
 L50:
  return;
 L60:
  for (i = 15; i < 17; ++i) {
    col[i] = 0.0;
    col[i + 3] = 0.0;
  }
  col[17] = pxyz[8] * massParticle[1] / col[5];  // p(1), z * m2/m = p(bar prime) 1, z 
  col[20] = pxyz[8] / col[0];  // vz velocity 
  for (i = 21; i < 24; ++i) col[i] = 0.0;
  // cross product p1 prime x v 
  pxyz[2] = col[9] * snt * sopc;
  // x component p3 bar = p3 prime x sin theta x cos phi 
  pxyz[6] = col[9] * snt * sops;
  // y comp. p3 bar = p3 prime x sin theta x sin phi 
  pxyz[10] = col[9] * cst;
  z = pxyz[8] / col[5];
  pxyz[10] = pxyz[10] + z * pxyz[8] * pxyz[10] / (col[0] + col[5]) + z * col[8];
  // z comp. p3 bar = p3 prime cos theta + (p1z sq * p3z (prime)cos 
  // theta / (e prime * (e + e prime)) + p1z * e3 prime / e prime 
  energy[2] = sqrt(sqr(pxyz[2]) + sqr(pxyz[6]) + sqr(pxyz[10]) + sqr(massParticle[2]));
  for (i = 0; i < 9; i += 4) pxyz[i + 3] = pxyz[i] - pxyz[i + 2];
  energy[3] = sqrt(sqr(pxyz[3]) + sqr(pxyz[7]) + sqr(pxyz[11]) + sqr(massParticle[3]));
  if (pt[37] != 0.0) goto L50;
  pt[2] = (energy[3] - massParticle[3]) / rcpmv + pt[2];
  goto L50;
}

void G4BertiniIsobarModel::qou17(G4double *t, 
			     G4double *b, 
			     G4double *r, 
			     G4double *w, 
			     G4double *g) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double unive;
  G4double a;
  G4double com1;
  G4double ftr;
  if (i3 <= 0) goto L10;
  else goto L160;
 L10:
  pt[37] = 0.0;
  // 1-alpha  particle 6+1 
  value1 = rlke - 180.0;
  qrdet(1, t, value1); //::: was t[1] check 
  com2 = crdt[0];
  ftr = 2 * massNucleon * rlke * rcpmv + 2.9877156e27;
  // sqr(e) = sqr(min) + ncnms * rlke * 2 * rcpmv 
  univer = sqrt(ftr); // e 
 L30:
  value2 = G4UniformRand();
  com = value2 * com2;
  // r - prime 
  qene(b); //::: 
  com1 = (sqr(com) + ftr - 5.01264e25) / (univer * 2.0);
  // m1r prime)^2 + sqr(e) - 2(massPionCharged)/2e = e alpha 
  a = sqr(com1) - sqr(com);
  if (a >= 0.0) goto L50;
  pacnt += 1.0;
  goto L30;
 L50:
  unive = (univer - com1) * com1 * (sqrt(a) / univer); // ((e beta * e alpha * p alpha) / e) = f(m, tr) 
  qrdet(1, rr, value1); // (pi, N) fmax(rlke)isobar sampling s.p. 
  com1 = G4UniformRand(); 
  if (unive / crdt[0] >= com1) goto L60; // random number <= f(m, tr) / fmax(tr) 
  else goto L30;
  
 L60:
  // if(curr(1).gt.0.0)write(io,2011)curr(1) 
  // 2011 format(1h ,2x,'qou17:calling qngid at 2010:curr(1)=',f10.2/) 
  qngid();
  massParticle[3] = massPionZero;
  massParticle[2] = com;
  pt[1 ] = 3.0;
  pt[3 ] = massPionZero;
  pt[13] = 3.0;
  pt[15] = massPionZero;
  pt[25] = 1.0;
  pt[27] = massNucleon;
  if (isw[8] == 0) {
    if (isw[9] != 0) goto L120;
    else goto L110;
  }
  if (isw[9] != 0) goto L130;
 L90:
  i3 = -1;
  return;
 L110:
  value1 = 0.4;
  value2 = twoThirds;
  value3 = 0.0;
  goto L150;
 L120:
  qrdet(2, w, value1); // (pi+-, p) fractional final state with recl. pi1 pi0 l.e. 
  value3 = oneThird;
  goto L140;
 L130:
  qrdet(2, g, value1); // (pi-, p) fract. 0 in.sta. with recl. pi1 pio l.e. 
  value3 = strkp; 
 L140:
  value1 = crdt[0];
  value2 = crdt[1];
 L150:
  qlpha();
 L160:
  pt[2 ] = 0.0;
  pt[14] = 0.0;
  pt[26] = 0.0;
  pt[38] = 0.0;
  qoll();
  if (col[14] != 0.0) goto L90;
  if (pt[37] != 0.0) {
    i3 = 0;
    return;
  }
  pt[38] = 0.0;
  if (ik != 0) {
    pt[2] = (energy[3] - massParticle[3]) / rcpmv + pt[2];
  }
  i3 = 1;
  return;
}

void G4BertiniIsobarModel::qou18() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int i; 
  G4int j;
  G4int k; 
  G4int l;
  switch (i3) {
  case 1:  goto L10;
  case 2:  goto L20;
  case 3:  goto L80;
  case 4:  goto L120;
  case 5:  goto L80;
  case 6:  goto L100;
  }
 L10:
  i = 3;
  col[14] = 1.0;
  k = 27;
  goto L30;
 L20:
  i = 3;
  col[14] = 4.0;
  k = 15;
 L30:
  pnidk[0] = massParticle[i - 1];
  j = i;
  for (l = 1; l < 4; ++l) {
    pnidk[l] = pxyz[j ];
    j += 4;
  }
  pnidk[4] = energy[i - 1];
  pnidk[5] = pt[k - 12];
  qdk();
  if (k != 27) goto L60;
  pt[14] += (pnidk[11] - pnidk[5]) / rcpmv;
 L60:
  pt[k - 1] += (pnidk[12] - massNucleon) / rcpmv;
  i3 = 1;
 L70:
  iv = k;
  return;
 L80:
  col[14] = 3.0;
  k = 15;
  if (pt[13] <= 2.0) goto L90;
  else goto L120;
 L90:
  i3 = 2;
  goto L70;
 L100:
  l = 14;
  G4int m;
  for (m = 4; m < 7; ++m) {
    pt[m] = pnidk[l];
    pt[m + 12] = pnidk[l + 3];
    ++l;
  }
  pt[10] = pnidk[11];
  pt[11] = pnidk[5];
  i = 4;
  k = 39;
  col[14] = 5.0;
  goto L30;
 L120:
  i1 = 3;
 L130:
  k = i1 * 12 - 33;
  if (i1 > 4) goto L170;
  if (i1 == 4) goto L160;
  i2 = -1;
  goto L200;
 L160:
  i2 = 0;
  goto L200;
 L170:
  if (i1 < 5) {
    goto L160;
  } else if (i1 == 5) {
    goto L190;
  }
  i3 = 4;
  goto L70;
 L190:
  i2 = 1;
 L200:
  if (pt[k - 1] == 0.0) goto L220;
  qstor();
 L220:
  ++i1;
  goto L130;
}

void G4BertiniIsobarModel:: qou19() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int k;
  pt[2] += (pt[10] - pt[11]) / rcpmv;
  // collision allowed 
  k = 3;
  goto L100;
 L40:
  i2 = 2;
 L50:
  i1 = k / 12 + 3;
  qstor();
  if (k > 15) goto L120;
  if (k != 15) {
    k = 15;
    goto L40;
  }
  k = 27;
  pt[26] += (pnidk[11] - pt[k]) / rcpmv;
 L100:
  if (k < 15) goto L40;
  i2 = 0;
  goto L50;
 L120:
  if (k < 27) {
    i3 = -1;
    return;
  }
  if (k > 27) goto L140;
  if (pt[38] <= 0.0) {
  L140:
    i3 = 0;
    return;
  }
  i2 = 1;
  k = 39;
  goto L50;
}

void G4BertiniIsobarModel::qou21(G4double *v,
			     G4double *w,
			     G4double *x,
			     G4double *y,
			     G4double *z) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double ans;
  G4double ftr;
  G4double sn;
  G4double com1;
  G4double com4;
  G4double a;

  value2 = rlke * 4.81633308e24 + 9.0554256e27;
  G4double value3 = sqrt(value2);
  switch (i3) {
  case 1:  goto L10;
  case 2:  goto L130;
  case 3:  goto L140;
  case 4:  goto L270;
  }
 L10:
  isw[11] = 0;
 L20:
  pt[37] = 0.0;
  i1 = 0;
  ans = rlke;
 L40:
  value1 = ans - 300.0;
  qrdet(1, v, value1); // (N, N) f(tr) isobar sampling 
  ftr = crdt[0];
  ++icoun1;
 L50:
  sn = G4UniformRand();
  ++icoun2;
  if (icoun2 > 50) {
    nopart = -1;
    return;
  }
  com = sn * ftr;
  qene(w); // (N, N) mass of isobar s.p. m(r prime) 
  if (i1 < 0) goto L160;
  if (i1 > 0) goto L170;
  com1 = (sqr(com) - sqnm + value2) / (value3 * 2.0); // e gamma 
  a = sqr(com1) - sqr(com);
  if (a >= 0.0) goto L80;
  pgcnt += 1.0;
  goto L50;
 L80:
  univer = sqrt(a) * com1 * (value3 - com1) / value3;
  qrdet(1, x, value1);
  com1 = G4UniformRand();
  if (com1 - univer / crdt[0] > 0.0) goto L50;
  massParticle[3] = massNucleon;
  massParticle[2] = com;
  // if(curr(1) > 0)write(io,4344)curr(1) 
  // 4344  format(1h ,2x,'qou21:calling qngid at 4345:curr(1)=',f10.2/) 
  qngid();
  pt[3]  = massNucleon;
  pt[27] = massNucleon;
  qlp19();
  return;
 L130:
  isw[11] = 2;
  goto L20;
 L140:
  isw[12] = 0;
 L150:
  i1 = -1;
  ans = (sqr(static_cast<G4double>(value3 - 7.08e12)) - 9.0554256e27) / 4.81633308e24;
  goto L40;
  // tr prime com1 = rlke prime 
 L160:
  com1 = (sqr(value3 + massNucleon - com) - 9.0554256e27) / 4.81633308e24;
  com2 = com;
  ans = com1;
  com4 = ftr;
  i1 = 1;
  goto L40;
 L170:
  com1 = (sqr(com2) - sqr(com) + value2) / (value3 * 2.0); // e epsilon 
  a = sqr(com1) - sqr(com2);
  if (a < 0.0) {
    pecnt += 1.0;
    goto L200;
    // f(m1, m2, tr) = p epsilon * e epsilon * e zeta/e 
  }
  univer = sqrt(a) * com1 * (value3 - com1) / value3;
  value1 = rlke - 920.0;
  qrdet(1, y, value1);
  sn = G4UniformRand();
  if (sn - univer * ftr / (crdt[0] * com4) <= 0.0) {
    goto L210;
  }
 L200:
  ftr = com4;
  i1 = -1;
  goto L50;
 L210:
  value1 = G4UniformRand();
  if (value1 <= 0.5) {
    massParticle[2] = com2;
    massParticle[3] = com;
    goto L240;
  }
  massParticle[2] = com;
  massParticle[3] = com2;
 L240:
  // if(curr(1).gt.0)write(io,4399)curr(1) 
  // 4399  format(1h ,2x,'qou21:calling qngid at 4399:curr(1)=',f10.2/) 
  qngid();
  // 4400 call qngid 
  pt[15] = massNucleon;
  pt[39] = massNucleon;
  if (isw[12] != 0) {
    qrdet(1, z, rlke);
    value1 = crdt[0]; // (n, p) fract. 0 in. sta. 3/2 l.e. 
  }
  pt[1] = 3.0;
  pt[3] = massPionZero;
  pt[13] = 1.0;
  pt[25] = 3.0;
  pt[27] = massPionZero;
  pt[37] = 1.0;
  qlp28();
  return;
 L270:
  isw[12] = 2; 
  goto L150;
}

//---------------------------------------------------------------------------------------
// methods ready for testing and further development

void G4BertiniIsobarModel::qrdet(G4int     nodata,  
			     G4double *data,  
			     G4double  energy) {
  // decription : 
  // parameters : nodata = data per energy interval, *data , energy
  // uses       :
  // changes    : crdt
  G4int ie = static_cast<G4int>(fabs(energy / 20.0)); // energy interval
  G4double univ = (energy - static_cast<G4double>(ie) * 20.0) / 20.0;
  for (G4int i = 0; i < 25; ++i) crdt[i] = 0.0;
  G4int k = nodata * ie + 1;
  G4int n;
  if (inpt != 0) {
    k = inpt - 1 + k;
    n = 2;
  } else  n = nodata;  // whole interval considered
  G4int l = k + nodata;
  for (i = 0; i < n; ++i) {
    crdt[i] = (data[l] - data[k]) * univ + data[k];
    ++k;
    ++l;
  }
  inpt = 0;
  return;
}

void G4BertiniIsobarModel::signex() {
  // decription: random exponential divided by sigma ci region i, ex = distance in sampling routine  
  // parameters:
  // uses:
  // changes:
  ex += exprn() / sign;
}

void G4BertiniIsobarModel::pol1(G4double &cs, 
			    G4double &si) {
  // decription : calculate uniformly random cos(phi), sin(phi) 
  // parameters : cs, si
  // uses       : -
  // changes    : cs, ci
  cs = G4UniformRand();
  if (G4UniformRand() < 0.5) cs = -cs;
  si = sqrt(1.0 - (cs * cs));
  return;
}

void G4BertiniIsobarModel::rout3() { 
  // decription: set isw[10] and call undis()
  // parameters:
  // uses:
  // changes:
  if (no < 4) isw[10] = 1;
  else isw[10] = 0;
  // undis(); // :::
  // g4bertini->undis();//:::: 
  inc = 1; 
}

void G4BertiniIsobarModel::rout4() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  geo();
  if (i1 >= 0) {
    curr[2] = massPionCharged;
    curr[0] = no;
    partin();
    //  g4bertini->spac32(32);//:::
  }
}

void G4BertiniIsobarModel::rout5(G4double *t, 
			     G4double *b,
			     G4double *rr) {
  // decription: calculate cross sections
  // parameters:
  // uses:
  // changes:
  if (not <  2) crjab(1,  t[1]);
  if (not == 2) crjab(1,  b[1]); // (pi+, p) elastic scattering cross-section
  if (not >  2) crjab(1, rr[1]); // (pi-, p)  direct scattering cross-section
  return;
}
	  
void G4BertiniIsobarModel::rout6a() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  for (;;){
    i1    = 0;
    spisom();
    strkp = - 2.0;
    i1    = 1;
    spisom();
    strkp = - 1.0;
    com = (awd[medium] - 7.0) * 2.0 * rcpmv;
    if (com > energy[2]) break;
  }
  massParticle[2] = 2.0 * massNucleon;
  massParticle[3] = massNucleon;
  energy[2]       = massParticle[2] + energy[2];
  return;
}

void G4BertiniIsobarModel::rou10() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  i3 = 0;
  if (ex > d[4]) {
    //    g4bertini->spac32(31); //:::
    i3 = - 1;
  } else {
    curr[2] = out[15];
    wkrpn[1] = out[15]; // kinetic energy for p and n in region 1
    wkrpn[4] = out[18]; 
  }
  return;
}

void G4BertiniIsobarModel::rou12() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  i3 = 0;
  azio(sopc, sops);
  cacoll(0);
  if (col[15] != 0) {
    i3 = - 1;
    return;
  } else if (ke > 0) {
    com = (energy[4] - massNucleon) / rcpmv;
  } else if (ke == 0) {
    cole4(abz);
  } else {
    i3=-1;
    return;
  }
  i1 =  - 1;
  value1 = com;
  if (pt[14] < 2.0) i3 = 1;
}

G4double G4BertiniIsobarModel::mud() {
  // decription: generate random variable sine 
  // parameters: -
  // uses: -
  // changes: inpt, sine
  // sine = (0.02 n - r) / .02 = n - r / .02   n = inpt   r/.02 = (n - 1) + x 
  G4double sine = 50.0 * G4UniformRand();
  inpt = static_cast<G4int>(sine + 1.0); 
  sine = static_cast<G4double>(inpt) - sine;
  return sine;
}

void G4BertiniIsobarModel::cole4(G4double abz) {  //::: 
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int k;
  k = static_cast<G4int>(clsm);
  com = (energy[3] - massNucleon) / rcpmv; // kinetic energy of particle 4 in MeV 
  // type of particle 1 - 5  current nucleon 
  if (curr[0] < 3.0) return;
  if (abz != 0) {
    com += absec;
    return;
  }  
  if (it != 5) { 
    univ = 1.0;
  } else if (it == 24) {
    univ = 0.0;
  } else  if (it == 6) {
    univ = 0.0;
  } else if (it == 26) {
    univ = 1.0;
  } else {
    return;
  }
  // neutron with wrong energy 
  unive = space[k + 2] - space[k + 8]; // region i (n, p) well depth difference 
  if (univ != 0.0) com -= unive;
  else com += unive;
}

G4double G4BertiniIsobarModel::big7(G4double c, 
				G4int ix) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // returns the biggest of the ix+1 random numbers
  // c is a initial value for the biggest number
  // note: we replace original dflran function with G4UniformRand 
  G4double great = c; 
  for (G4int i = 1; i < (ix + 1); ++i) {
    c = G4UniformRand();
    if (c >= great) great = c;
  }
  return great;
} 
G4double G4BertiniIsobarModel::angle[15] = 
{
  0.0           , 0.0           , 1.0  * degree,  
  2.0  * degree , 4.0  * degree , 8.0  * degree, 
  12.0 * degree , 20.0 * degree , 45.0 * degree,
  0.0           , 0.0           , 0.0          ,
  0.0           , 0.0           , 0.0
};
											 
G4double G4BertiniIsobarModel::pmxxx[12] = 
{
  4.5, 		7.5, 	10.5, 
  14.0, 	22.0, 	34.0, 
  50.0, 	90.0, 	170.0,
  290.0, 	490.0, 1000.0
};
								
G4double G4BertiniIsobarModel::dndpip[12][3] = 
{
  {-5.3242445e-03,   2.2346497e-02,    4.5267868e-01},
  {1.1044145e-03,   -3.2094955e-02,    5.6748962e-01},
  {8.9269876e-04,   -2.8545380e-02,    5.5277443e-01},
  {5.0222874e-04,   -2.0401001e-02,    5.1032925e-01},
  {2.0570680e-04,   -1.1921406e-02,    4.4970989e-01},
  {6.6768611e-05,   -5.8776736e-03,    3.8399315e-01},
  {2.2491440e-05,   -2.9177666e-03,    3.3454037e-01},
  {5.6584831e-06,   -1.2033470e-03,    2.9090124e-01},
  {9.4826100e-07,   -3.7099421e-04,    2.5414348e-01},
  {1.6872946e-07,   -1.1652336e-04,    2.3341179e-01},
  {3.4642653e-08,   -4.0553510e-05,    2.2265625e-01},
  {5.3050826e-09,   -1.1835014e-05,    2.1562809e-01}
};

G4double G4BertiniIsobarModel::dndpim[12][3] = 
{
  {-3.7355423e-03,   2.9491067e-02,    2.0408440e-01},
  {-1.3875589e-04,  -1.4313459e-03,    2.7039909e-01},
  {1.4499575e-04,   -5.3753257e-03,    2.8402138e-01},
  {1.0721385e-04,   -4.5639873e-03,    2.7966785e-01},
  {4.8756599e-05,   -2.8809309e-03,    2.6756477e-01},
  {1.6727019e-05,   -1.4842749e-03,    2.5234318e-01},
  {5.7825819e-06,   -7.5179338e-04,    2.4009037e-01},
  {1.4668331e-06,   -3.1205267e-04,    2.2888833e-01},
  {2.4605833e-07,   -9.6268952e-05,    2.1935964e-01},
  {4.3786713e-08,   -3.0238181e-05,    2.1398067e-01},
  {8.9794412e-09,   -1.0516495e-05,    2.1118736e-01},
  {1.3733370e-09,   -3.0666124e-06,    2.0936286e-01}
};

G4double G4BertiniIsobarModel::pz[12] =  
{
  4.5,		7.5,	10.5, 
  14.0,		22.0,	34.0, 
  50.0,		90.0,	170.0, 
  290.0,	490.0,	1000.0
};

G4double G4BertiniIsobarModel::cpnu[12][3] = 
{
  {-1.3039768e-02,    1.3335037e-01,    6.7083645e-01},
  {-2.2712350e-03,    4.0141106e-02,    8.7221718e-01}, 
  {-6.1076880e-04,    1.6270638e-02,    9.5784855e-01}, 
  {-2.3752451e-04,    8.5802078e-03,    9.9742126e-01}, 
  {-7.3343515e-05,    3.9319992e-03,    1.0303373    }, 
  {-1.8943101e-05,    1.5830994e-03,    1.0556831    }, 
  {-5.4761767e-06,    6.8706274e-04,    1.0705557    }, 
  {-1.2256205e-06,    2.5695562e-04,    1.0814734    }, 
  {-1.9255094e-07,    7.5042248e-05,    1.0894747    }, 
  {-3.3469405e-08,    2.3089349e-05,    1.0937023    }, 
  {-6.9267116e-09,    8.0280006e-06,    1.0958271    }, 
  {-1.0377335e-09,    2.3245811e-06,    1.0972347    }
};

G4double G4BertiniIsobarModel::cpk[12][3] =
{
  {-1.4113426e-02,    5.3243160e-01,    3.4706497e-01},
  {-3.8183928e-03,    4.4261360e-01,    5.4272461e-01}, 
  {-1.6031265e-03,    4.1049194e-01,    6.5904236e-01}, 
  {-8.7887049e-04,    3.9549255e-01,    7.3689270e-01}, 
  {-4.1687489e-04,    3.8226700e-01,    8.3128357e-01}, 
  {-1.7768145e-04,    3.7179470e-01,    9.4566345e-01}, 
  {-8.0645084e-05,    3.6527634e-01,    1.0555573    }, 
  {-3.0577183e-05,    3.6011600e-01,    1.1883698    }, 
  {-9.1493130e-06,    3.5625267e-01,    1.3625183    }, 
  {-2.9318035e-06,    3.5418224e-01,    1.5346680    }, 
  {-1.0281801e-06,    3.5308743e-01,    1.6921387    }, 
  {-2.9010698e-07,    3.5234547e-01,    1.8784180    }
};

G4double G4BertiniIsobarModel::cpipnu[12][3] = 
{
  {-1.3925135e-02,    1.9439220e-01,    2.6330948e-02},
  {-4.9353242e-03,    1.1407471e-01,    2.0571613e-01}, 
  {-1.9524097e-03,    7.0723534e-02,    3.6305332e-01}, 
  {-9.7447634e-04,    5.0484657e-02,    4.6774292e-01}, 
  {-4.1649491e-04,    3.4583092e-02,    5.8101559e-01}, 
  {-1.6780943e-04,    2.3748994e-02,    6.9900227e-01}, 
  {-7.8253448e-05,    1.7732978e-02,    7.9998779e-01}, 
  {-3.3177203e-05,    1.3085365e-02,    9.1971684e-01}, 
  {-1.2346776e-05,    9.3094110e-03,    1.0908327    }, 
  {-5.0680246e-06,    6.8699718e-03,    1.2951612    }, 
  {-2.2703025e-06,    5.2463412e-03,    1.5307465    }, 
  {-8.7339140e-07,    3.8301349e-03,    1.8892488    }
};

G4double G4BertiniIsobarModel::cpipk[12][3] = 
{
  {-5.6201220e-04,    1.8275261e-01,   -1.3130379e-01},
  {-9.0676546e-04,    1.8500519e-01,   -1.3445663e-01}, 
  {-4.1097403e-04,    1.7770004e-01,   -1.0755920e-01}, 
  {-1.9192696e-04,    1.7316246e-01,   -8.4091187e-02}, 
  {-6.6518784e-05,    1.6959763e-01,   -5.8746338e-02}, 
  {-1.8179417e-05,    1.6750813e-01,   -3.6132813e-02}, 
  {-5.3644180e-06,    1.6665268e-01,   -2.1896362e-02}, 
  {-1.1995435e-06,    1.6622925e-01,   -1.1184692e-02}, 
  {-1.8626451e-07,    1.6604996e-01,   -3.3721924e-03}, 
  {-2.9802322e-08,    1.6599941e-01,    4.8828125e-04}, 
  {-3.7252903e-09,    1.6598701e-01,    2.6855469e-03}, 
  {-1.1641532e-09,    1.6598129e-01,    3.6621094e-03}
};

G4double G4BertiniIsobarModel::cpimnu[12][3] = 
{
  {-6.8021417e-03,    1.0507488e-01,    2.5606155e-02},
  {-2.4644732e-03,    6.6451073e-02,    1.1157703e-01}, 
  {-1.0727644e-03,    4.6208382e-02,    1.8512058e-01}, 
  {-5.9020519e-04,    3.6210060e-02,    2.3689651e-01}, 
  {-2.8741732e-04,    2.7553260e-02,    2.9872513e-01}, 
  {-1.3258681e-04,    2.0779014e-02,    3.7283039e-01}, 
  {-6.7789108e-05,    1.6413033e-02,    4.4637489e-01}, 
  {-3.0794181e-05,    1.2584507e-02,    5.4531479e-01}, 
  {-1.1968659e-05,    9.1618299e-03,    7.0087528e-01}, 
  {-5.0016679e-06,    6.8241954e-03,    8.9691067e-01}, 
  {-2.2565946e-06,    5.2303076e-03,    1.1282806    }, 
  {-8.7132503e-07,    3.8254932e-03,    1.4840240    }
};

G4double G4BertiniIsobarModel::cpimk[12][3] = 
{
  {-9.9664927e-04,    1.0531235e-01,   -6.1017036e-02},
  {-4.5859814e-04,    1.0030270e-01,   -4.9365997e-02}, 
  {-1.6081333e-04,    9.5970154e-02,   -3.3605576e-02}, 
  {-6.6816807e-05,    9.4027519e-02,   -2.3590088e-02}, 
  {-2.1457672e-05,    9.2741013e-02,   -1.4461517e-02}, 
  {-5.5730343e-06,    9.2055321e-02,   -7.0800781e-03}, 
  {-1.6018748e-06,    9.1790199e-02,   -2.6702881e-03}, 
  {-3.5762787e-07,    9.1663361e-02,    5.1879883e-04}, 
  {-5.2154064e-08,    9.1610909e-02,    2.8533936e-03}, 
  {-5.8207661e-09,    9.1595650e-02,    4.2114258e-03}, 
  {-1.8626451e-09,    9.1592789e-02,    4.7149658e-03}, 
  {-4.6566129e-10,    9.1590405e-02,    4.7454834e-03}
};






