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
#include "G4BertiniElasticCollision.hh"

G4BertiniElasticCollision::G4BertiniElasticCollision() {
  ;
}

G4BertiniElasticCollision::~G4BertiniElasticCollision(){
  ;
}

void G4BertiniElasticCollision::interpolateElasticNeutronData(G4int medium, 
							      G4int kdd, 
							      G4double theEnergy) {
  // decription: interpolates on ec to set fone, the coefficient of the 
  //             first Legendre polynomial in the scattering angle distribution
  //             for the elastic neutron scattering 
  // parameters: medium, kdd, theEnergy
  // uses:
  // changes: totels, fone den[10][17] not given in parameter list
  G4int loc;
  G4int npts;
  G4int i;
  G4double eff;
  G4double sg;
  G4double ess;
  if (kdd > 0) {
    idd = kdd;
    loc = locf1[idd];
    npts = static_cast<G4int>(locf1[idd + 1] - locf1[idd]);
    geti(ef[loc], npts, theEnergy, i); 
    loc += i - 1;
    fone = f1[loc];
    eff  = ef[loc];
    fone += (f1[loc + 1] - fone) * (theEnergy - eff) / (ef[loc + 1] - eff);
  } else {
    nn = noel[medium];
    for (G4int l = 0; l < nn; ++l) {
      idd  = id[l][medium];
      loc  = locsig[idd];
      npts = locsig[idd + 1] - locsig[idd];
      geti(es1[loc], npts, theEnergy, i);
      loc += i - 1;
      sg   = sige[loc];
      ess  = es1[loc];
      sg  += (sige[loc + 1] - sg) * (theEnergy - ess) / (es1[loc + 1] - ess);
      sgels[l] = sg * den[l][medium];
      totels += sgels[l];
    }
  }
}




void G4BertiniElasticCollision::geti(G4double es, 
				     G4int npts, 
				     G4double theEnergy,
				     G4int i) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // status: verify algorithm
  if (es1[0] < theEnergy) {
    if (verboseLevel > 1) {
      G4cout << "error in geti"      << G4endl;
      G4cout << "es1[0] :" << es1[0] << G4endl;
      G4cout << "e      :" << theEnergy      << G4endl;
    }
    G4Exception("error in geti");
    return;
  } else if (es1[0] == theEnergy) {
    i = 1; 
    return;
  } else if (es1[0] > theEnergy) {
    G4int itop = 1;
    G4int ibot = npts;
    for(;;){
      if ((ibot - itop) == 1) {
	i = itop;
	return;
      } else {
	G4int i = (ibot + itop) / 2;
	if (es1[i] < theEnergy) {
	  ibot = i; 
	} else if (es1[i] == theEnergy) {
          --i;
	  return;
	} else if (es1[i] > theEnergy) {
	  itop = i;
	}
      }
    }
  }
}

void G4BertiniElasticCollision::scatteringWithHydrogen() {
  // decription : elastic scattering with hydrogen, struck particle at rest
  // parameters :
  // uses       :
  // changes    :
  G4double sopc;
  G4double sops;
  azio(sopc, sops);
  massParticle[3]   = massNucleon;
  pt[0]   = 0.0;
  for (G4int i =  2; i < 13; ++i)  pt[i] = 0.0;
  for (i = 14; i < 48; ++i)  pt[i] = 0.0;
  for (i =  0; i < 24; ++i) col[i] = 0.0;
  G4double a = sqr(massParticle[3]);
  col[0]  = energy[0] + energy[1]; // total energy particles 1 and 2
  for (i = 0; i < 3; ++i) col[i + 1] = massParticle[i] * massParticle[i]; // mass particle i sqd.
  col[5]  = col[3] + col[2] + 2.0 * (energy[1] * energy[2] - (pxyz[1] * pxyz[2] + pxyz[5] * pxyz[6] + pxyz[9] * pxyz[10]));
  col[6]  = sqrt(col[5]);
  col[7]  = col[6]   / col[1];         // gamma
  col[8]  = 2.0 * col[6];
  col[9]  = (col[4] + col[5] - a) / col[8];
  G4double com2 = sqr(col[9]);
  col[10] = sqrt(com2 - col[4]);       // p3 prime
  col[18] = pxyz[9]  * massParticle[2] / col[6]; // p[1], z * m2/m = p[bar prime] 1, z
  col[21] = pxyz[9]  / col[1];         // vz velocity
  pxyz[3] = col[10]  * snt * sopc;     // x component p3 bar =p3 prime x sin theta x cos phi
  pxyz[7] = col[10]  * snt * sops;     // y comp. p3 bar =p3 prime x sin theta x sin phi
  pxyz[11] = col[10] * cst;
  G4double z = pxyz[9] / col[6];
  pxyz[11] = pxyz[11] + (z * pxyz[9] * pxyz[11]) / (col[1] + col[6]) + z * col[9];
  // z comp. p3 bar = p3 prime cos theta + (p1z sq * p3z [prime] cos
  // theta / (e prime * (e + e prime]] + p1z * e3 prime / e prime
  energy[3] = sqrt(sqr(pxyz[3]) + sqr(pxyz[7]) + sqr(pxyz[11]) + sqr(massParticle[3]));
  for (i = 0; i < 9; i += 4) pxyz[i + 3] = pxyz[i] - pxyz[i + 2];
  energy[4] = sqrt(sqr(pxyz[4]) + sqr(pxyz[8]) + sqr(pxyz[12]) + sqr(massParticle[4])); // kinetic energy of 3 and 4
  pt[3]   = (energy[3] - massParticle[3]) / rcpmv;   // momentum of 3 and 4, in laboratory
  pt[15]  = (energy[4] - massParticle[4]) / rcpmv;     
  G4double p3 = sqrt(sqr(pxyz[3]) + sqr(pxyz[7]) + sqr(pxyz[11])); // direct cosines of 3 and 4
  G4double p4 = sqrt(sqr(pxyz[4]) + sqr(pxyz[8]) + sqr(pxyz[12])); 
  pt[8]   = pxyz[3]  / p3;
  pt[9]   = pxyz[7]  / p3;
  pt[10]  = pxyz[11] / p3;
  pt[20]  = pxyz[4]  / p4;
  pt[21]  = pxyz[8]  / p4;
  pt[22]  = pxyz[12] / p4;
}








