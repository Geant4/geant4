#include "G4InuclSpecialFunctions.hh"
#include "Randomize.hh"

G4double G4InuclSpecialFunctions::getAL(G4double A) {
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::getAL" << G4endl;
  }

  return 0.76 + 2.2 / pow(A, 0.333333);
}

G4double G4InuclSpecialFunctions::csNN(G4double e) {
  G4int verboseLevel = 2;
 
  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::csNN" << G4endl;
  }

  G4double snn;

  if (e < 40.0) {
    snn = -1174.8 / (e * e) + 3088.5 / e + 5.3107;

  } else {
    snn = 93074.0 / (e * e) - 11.148 / e + 22.429;
  };

  return snn; 
}

G4double G4InuclSpecialFunctions::csPN(G4double e) {
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::csPN" << G4endl;
  }

  G4double spn;

  if (e < 40.0) {
    spn = -5057.4 / (e * e) + 9069.2 / e + 6.9466;

  } else {
    spn = 239380.0 / (e * e) + 1802.0 / e + 27.147;
  };

  return spn; 
}

G4double G4InuclSpecialFunctions::FermiEnergy(G4double A, G4double Z, G4int ntype) {
  // calculates the nuclei Fermi energy for 0 - neutron and 1 - proton
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::FermiEnergy" << G4endl;
  }

  const G4double C = 55.4;
  G4double Ef;

  if (ntype == 0) {
    Ef = C * pow((A - Z) / A, 0.666667);

  } else {
    Ef = C * pow(Z / A, 0.666667);
  };

  return Ef; 
}

G4double G4InuclSpecialFunctions::inuclRndm() { 
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::inuclRndm" << G4endl;
  }

  G4double rnd = G4UniformRand(); 

  return rnd;
} 

G4double G4InuclSpecialFunctions::randomGauss(G4double sigma) {
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::randomGauss" << G4endl;
  }

  const G4double eps = 1.0e-6;
  const G4double twopi = 6.2831854;
  G4double r1 = inuclRndm();
  r1 = r1 > eps ? r1 : eps;
  G4double r2 = inuclRndm();
  r2 = r2 > eps ? r2 : eps;
  r2 = r2 < 1.0 - eps ? r2 : 1.0 - eps; 

  return sigma * sin(twopi * r1) * sqrt(-2.0 * log(r2)); 
} 

G4double G4InuclSpecialFunctions::randomPHI() { 
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::randomPHI" << G4endl;
  }

  const G4double twopi = 6.2831853;

  return twopi * inuclRndm();
} 

G4std::pair<G4double, G4double> G4InuclSpecialFunctions::randomCOS_SIN() {
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::randomCOS_SIN" << G4endl;
  }

  G4double CT = 1.0 - 2.0 * inuclRndm();

  return G4std::pair<G4double, G4double>(CT, sqrt(1.0 - CT * CT));
}

G4std::vector<G4double> G4InuclSpecialFunctions::generateWithFixedTheta(G4double ct,
									G4double p) {
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::generateWithFixedTheta" << G4endl;
  }

  G4std::vector<G4double> momr(4);
  G4double phi = randomPHI();
  G4double pt = p * sqrt(fabs(1.0 - ct * ct));
  G4std::vector<G4double> mom1(4);
  momr[1] = pt * cos(phi);
  momr[2] = pt * sin(phi);
  momr[3] = p * ct;

  return momr;
}
