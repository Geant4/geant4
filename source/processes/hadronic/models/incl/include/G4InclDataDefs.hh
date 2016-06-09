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
// $Id: G4InclDataDefs.hh,v 1.5 2008/06/25 17:20:04 kaitanie Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// All data structures needed by INCL4 are defined here.

#ifndef InclDataDefs_hh
#define InclDataDefs_hh 1

#define FSIZE 15
/**
 * Initial values of a hadronic cascade problem.
 */
class G4Calincl {
public:
  G4Calincl() {};
  ~G4Calincl() {};
  
  /**
   * Here f is an array containing the following initial values:
   * - f[0] : target mass number
   * - f[1] : target charge number
   * - f[2] : bullet energy
   * - f[3] : minimum proton energy to leave the target (default: 0.0)
   * - f[4] : nuclear potential (default: 45.0 MeV)
   * - f[5] : time scale (default: 1.0)
   * - f[6] : bullet type (1: proton, 2: neutron, 3: pi+, 4: pi0 5: pi-, 6:H2, 7: H3, 8: He3, 9: He4
   * - f[7] : minimum neutron energy to leave the target (default: 0.0)
   * - f[8] : target material identifier (G4Mat)
   * - f[9] : not used
   * - f[10] : not used
   * - f[11] : not used
   * - f[12] : not used
   * - f[13] : not used
   * - f[14] : not used
   */
  G4double f[FSIZE];

  /**
   * Number of events to be processed.
   */
  G4int icoup;
};

#define IGRAINESIZE 19
/**
 * Random seeds used by internal random number generators.
 * @see G4Incl::standardRandom
 * @see G4Incl::gaussianRandom
 */
class G4Hazard{
public:
  G4Hazard() {};
  ~G4Hazard() {};

  /**
   * Random seed
   */ 
  G4long ial;

  /**
   * An array of random seeds.
   */
  G4long igraine[IGRAINESIZE];
};

#define MATSIZE 500
#define MATGEOSIZE 6
/**
 * Target nuclei to be taken into account in the cascade problem.
 */
class G4Mat {
public:
  G4Mat() { };
  ~G4Mat() { };

  /**
   * Charge numbers.
   */
  G4int zmat[MATSIZE];

  /**
   * Mass number
   */ 
  G4int amat[MATSIZE];

  /**
   *
   */
  G4double bmax_geo[MATGEOSIZE][MATSIZE];

  /**
   * Number of materials.
   */
  G4int nbmat;
};

#define LGNSIZE 9
/**
 * Properties of light nucleus used as a bullet.
 */
class G4LightGausNuc {
public:
  G4LightGausNuc() {};
  ~G4LightGausNuc() {};
  
  G4double rms1t[LGNSIZE];
  G4double pf1t[LGNSIZE];
  G4double pfln[LGNSIZE];
  G4double tfln[LGNSIZE];
  G4double vnuc[LGNSIZE];
};

#define LNSIZE 30
/**
 * Data of light nuclei.
 */
class G4LightNuc {
public:
  G4LightNuc() {};
  ~G4LightNuc() {};

  /**
   * r
   */
  G4double r[LNSIZE];

  /**
   * a
   */
  G4double a[LNSIZE];
};

#define SAXWROWS 30 
#define SAXWCOLS 500
/**
 * Woods-Saxon density and its first derivative.
 */
class G4Saxw {
public:
  G4Saxw() {};
  ~G4Saxw() {};
  
  /**
   * x
   */
  G4double x[SAXWROWS][SAXWCOLS]; 

  /**
   * y
   */
  G4double y[SAXWROWS][SAXWCOLS]; 

  /**
   * s
   */
  G4double s[SAXWROWS][SAXWCOLS]; 

  /**
   * imat
   */
  G4int imat;

  /**
   * n
   */
  G4int n;

  /**
   * k
   */
  G4int k;
};

/**
 * Parameters for INCL4 model.
 */
class G4Ws {
public:
  G4Ws() {};
  ~G4Ws() {};
  
  /**
   * r0
   */
  G4double r0;

  /**
   * adif
   */
  G4double adif;

  /**
   * Maximum radius of the nucleus
   */
  G4double rmaxws;

  /**
   * drws
   */
  G4double drws;

  /**
   * Shape of the surface of the nucleus:
   * - -1: Woods-Saxon density with impact parameter dependence
   * -  0: Woods-Saxon density without impact parameter dependence
   * -  1: Sharp surface (hard sphere)
   */
  G4double nosurf;

  /**
   * Parameter related to the maximum radius of the nucleus.
   * 
   * rmaxws = r0 + xfoisa*A
   */
  G4double xfoisa;

  /**
   * Pauli blocking used in the simulation:
   * - 0: statistic Pauli blocking
   * - 1: strict Pauli blocking
   * - 2: no Pauli blocking
   */
  G4double npaulstr;

  /**
   * Maximum impact parameter
   */
  G4double bmax;
};

#define DTONSIZE 13
/**
 * Random seeds used by internal random number generators.
 * @see G4Incl::standardRandom
 * @see G4Incl::gaussianRandom
 */
class G4Dton {
public:
  G4Dton() {};
  ~G4Dton() {};
  
  G4double c[DTONSIZE];
  G4double d[DTONSIZE];
  G4double fn;
};

#define SPL2SIZE 100
/**
 * Random seeds used by internal random number generators.
 * @see G4Incl::standardRandom
 * @see G4Incl::gaussianRandom
 */
class G4Spl2 {
public:
  G4Spl2() {};
  ~G4Spl2() {};
  
  G4double x[SPL2SIZE];
  G4double y[SPL2SIZE];
  G4double a[SPL2SIZE];
  G4double b[SPL2SIZE];
  G4double c[SPL2SIZE];
  G4int n;
};

// incl4.2.cc:

//#define BL1SIZE 300
#define BL1SIZE 3000
/**
 * Random seeds used by internal random number generators.
 * @see G4Incl::standardRandom
 * @see G4Incl::gaussianRandom
 */
class G4Bl1 {
public:
  G4Bl1() {};
  ~G4Bl1() {};
  
  G4double p1[BL1SIZE],p2[BL1SIZE],p3[BL1SIZE];
  G4double eps[BL1SIZE];
  G4int ind1[BL1SIZE],ind2[BL1SIZE];
  G4double ta;
};

#define BL2CROISSIZE 19900
#define BL2INDSIZE 19900
/**
 * 
 */
class G4Bl2 {
public:
  G4Bl2() {};
  ~G4Bl2() {};
  
  /**
   * 
   */
  G4double crois[BL2CROISSIZE];

  /**
   *
   */
  G4int k;

  /**
   *
   */
  G4int ind[BL2INDSIZE];

  /**
   *
   */
  G4int jnd[BL2INDSIZE];
};

//#define BL3SIZE 300
#define BL3SIZE 3000
/**
 *
 */
class G4Bl3 {
public:
  G4Bl3() {};
  ~G4Bl3() {};
  
  /**
   * r1 and r2
   */
  G4double r1,r2;

  /**
   * Nucleon positions
   */
  G4double x1[BL3SIZE], x2[BL3SIZE],x3[BL3SIZE];

  /**
   * Mass numbers
   */
  G4int ia1,ia2;

  /**
   * rab2
   */
  G4double rab2;
};

/**
 * G4Bl4
 */
class G4Bl4 {
public:
  G4Bl4() {};
  ~G4Bl4() {};

  /**
   * tmax5
   */
  G4double tmax5;
};

//#define BL5SIZE 300
#define BL5SIZE 3000
/**
 * G4Bl5
 */
class G4Bl5 {
public:
  G4Bl5() {};
  ~G4Bl5() {};
  
  /**
   * tlg
   */
  G4double tlg[BL5SIZE];

  /**
   * nesc
   */
  G4int nesc[BL5SIZE];
};

/**
 * G4Bl6
 */
class G4Bl6 {
public:
  G4Bl6() {};
  ~G4Bl6() {};
  
  /**
   * xx10
   */
  G4double xx10;

  /**
   * isa
   */
  G4double isa;
};

/**
 * G4Bl8
 */
class G4Bl8 {
public:
  G4Bl8() {};
  ~G4Bl8() {};

  /**
   * rathr
   */
  G4double rathr;

  /**
   * ramass
   */
  G4double ramass;
};

//#define BL9SIZE 300
#define BL9SIZE 3000
/**
 * G4Bl9
 */
class G4Bl9 {
public:
  G4Bl9() {
    l1 = 0;
    l2 = 0;
  };
  ~G4Bl9() {};

  /**
   * hel
   */
  G4double hel[BL9SIZE];

  /**
   * l1 and l2
   */
  G4int l1,l2;
};

/**
 * G4Bl10
 */
class G4Bl10 {
public:
  G4Bl10() {};
  ~G4Bl10() {};

  /**
   * ri4, rs4, r2i, r2s, pdummy, pf
   */
  G4double ri4,rs4,r2i,r2s,pdummy,pf;
};

/**
 * G4Kind
 */
class G4Kind {
public:
  G4Kind() {};
  ~G4Kind() {};

  /**
   * kindf7
   */
  G4int kindf7;
};

#define VARSIZE 3
#define VAEPSSIZE 250
#define VAAVM 1000
/**
 * Extra information on collisions between nucleons.
 */
class G4VarAvat {
public:
  G4VarAvat() {};
  ~G4VarAvat() {};

  /**
   *
   */
  G4int kveux;

  /**
   *
   */
  G4double bavat;

  /**
   *
   */
  G4int nopartavat,ncolavat;

  /**
   *
   */
  G4double r1_in[VARSIZE],r1_first_avat[VARSIZE];

  /**
   *
   */
  G4double epsd[VAEPSSIZE],eps2[VAEPSSIZE],eps4[VAEPSSIZE],eps6[VAEPSSIZE],epsf[VAEPSSIZE];

  /**
   *
   */
  G4int nb_avat;

  /**
   *
   */
  G4double timeavat[VAAVM],l1avat[VAAVM],l2avat[VAAVM],jpartl1[VAAVM],jpartl2[VAAVM];

  /**
   *
   */
  G4double del1avat[VAAVM],del2avat[VAAVM],energyavat[VAAVM];

  /**
   *
   */
  G4double bloc_paul[VAAVM],bloc_cdpp[VAAVM],go_out[VAAVM];
};

#define VARNTPSIZE 255
class G4VarNtp {
public:
  G4VarNtp() {};
  ~G4VarNtp() {};

  /**
   * Clear and initialize all variables and arrays.
   */
  void clear() {
    particleIndex = 0;
    projType = 0;
    projEnergy = 0.0;
    targetA = 0;
    targetZ = 0;
    massini = 0;
    mzini = 0;
    exini = 0;
    pcorem = 0;
    mcorem = 0;
    pxrem = 0;
    pyrem = 0;
    pzrem = 0;
    mulncasc = 0;
    mulnevap = 0;
    mulntot = 0;
    bimpact = 0.0;
    jremn = 0;
    kfis = 0;
    estfis = 0;
    izfis = 0;
    iafis = 0;
    ntrack = 0;
    for(G4int i = 0; i < VARNTPSIZE; i++) {
      itypcasc[i] = 0;
      avv[i] = 0;
      zvv[i] = 0;
      enerj[i] = 0.0;
      plab[i] = 0.0;
      tetlab[i] = 0.0;
      philab[i] = 0.0;
      full[i] = false;
    }
  }

  void addParticle(G4double A, G4double Z, G4double E, G4double P, G4double theta, G4double phi) {
    if(full[particleIndex]) {
      G4cout <<"G4VarNtp: Error. Index i = " << particleIndex << " is already occupied by particle:" << G4endl;
      G4cout <<"A = " << avv[particleIndex] << " Z = " << zvv[particleIndex] << G4endl;
      G4cout <<"Tried to replace it with:" << G4endl;
      G4cout <<"A = " << Z << " Z = " << Z << G4endl;
    } else {
      avv[particleIndex] = (int) A;
      zvv[particleIndex] = (int) Z;
      enerj[particleIndex] = E;
      plab[particleIndex] = P;
      tetlab[particleIndex] = theta;
      philab[particleIndex] = phi;
      full[particleIndex] = true;
      ntrack = particleIndex + 1;
      particleIndex++;
    }
  }

  /**
   * Baryon number conservation check.
   */
  G4int getTotalBaryonNumber() {
    G4int baryonNumber = 0;
    for(G4int i = 0; i < ntrack; i++) {
      if(avv[i] > 0) {
	baryonNumber += avv[i];
      }
    }
    return baryonNumber;
  }

  /**
   * Return total energy.
   */
  G4double getTotalEnergy() {
    G4double energy = 0.0;
    for(G4int i = 0; i < ntrack; i++) {
      energy += std::sqrt(std::pow(plab[i], 2) + std::pow(getMass(i), 2)); // E^2 = p^2 + m^2
    }

    return energy;
  }

  /**
   * Return total three momentum.
   */
  G4double getTotalThreeMomentum() {
    G4double momentum = 0;
    for(G4int i = 0; i < ntrack; i++) {
      momentum += plab[i];
    }
    return momentum;
  }

  G4double getMomentumSum() {
    G4double momentum = 0;
    for(G4int i = 0; i < ntrack; i++) {
      momentum += plab[i];
    }
    return momentum;
  }

  G4double getMass(G4int particle) {
    const G4double protonMass = 938.272;
    const G4double neutronMass = 939.565;
    const G4double pionMass = 139.57;

    G4double mass = 0.0;
    if(avv[particle] ==  1 && zvv[particle] ==  1) mass = protonMass;
    if(avv[particle] ==  1 && zvv[particle] ==  0) mass = neutronMass;
    if(avv[particle] == -1)                        mass = pionMass;
    if(avv[particle] > 1)
      mass = avv[particle] * protonMass + zvv[particle] * neutronMass;
    return mass;
  }

  /**
   * Dump debugging output.
   */
  void dump() {
    G4int nProton = 0, nNeutron = 0;
    G4int nPiPlus = 0, nPiZero = 0, nPiMinus = 0;
    G4int nH2 = 0, nHe3 = 0, nAlpha = 0;
    G4int nFragments = 0;
    G4int nParticles = 0;
    G4cout <<"Particles produced in the event (" << ntrack << "):" << G4endl;
    G4cout <<"A \t Z \t Ekin \t Ptot \t Theta \t Phi" << G4endl;
    for(G4int i = 0; i < ntrack; i++) {
      nParticles++;
      if(avv[i] ==  1 && zvv[i] ==  1) nProton++;  // Count multiplicities
      if(avv[i] ==  1 && zvv[i] ==  0) nNeutron++;
      if(avv[i] == -1 && zvv[i] ==  1) nPiPlus++;
      if(avv[i] == -1 && zvv[i] ==  0) nPiZero++;
      if(avv[i] == -1 && zvv[i] == -1) nPiMinus++;
      if(avv[i] ==  2 && zvv[i] ==  1) nH2++;
      if(avv[i] ==  3 && zvv[i] ==  2) nHe3++;
      if(avv[i] ==  4 && zvv[i] ==  2) nAlpha++;
      if(                zvv[i] >   2) nFragments++;

      G4cout << i << " \t " << avv[i] << " \t " << zvv[i] << " \t " << enerj[i] << " \t " 
             << plab[i] << " \t " << tetlab[i] << " \t " << philab[i] << G4endl;
    }

    G4cout <<"Summary of event: " << G4endl;
    G4cout <<"Projectile type: " << projType <<" Energy: " << projEnergy << G4endl;
    G4cout <<"Target A = " << targetA << " Z = " << targetZ << G4endl;
    G4cout <<"Remnant from cascade: " << G4endl;
    G4cout <<"A = " << massini << " Z = " << mzini << " excitation E = " << exini << G4endl;
    G4cout <<"Particle multiplicities:" << G4endl;
    G4cout <<"Protons: " << nProton << " Neutrons:  " << nNeutron << G4endl;
    G4cout <<"pi+: " << nPiPlus << " pi0: " << nPiZero << " pi-: " << nPiMinus << G4endl;
    G4cout <<"H2: " << nH2 << " He3: " << nHe3 << " Alpha: " << nAlpha << G4endl;
    G4cout <<"Nucleus fragments = " << nFragments << G4endl;
    G4cout <<"Conservation laws:" << G4endl;
    G4cout <<"Baryon number = " <<  getTotalBaryonNumber() << G4endl;
    G4cout <<"Number of particles = " << nParticles << G4endl;
  }

  /**
   * Projectile type.
   */
  G4int projType;

  /**
   * Projectile energy.
   */
  G4double projEnergy;

  /**
   * Target mass number.
   */
  G4int targetA;

  /**
   * Target charge number.
   */
  G4int targetZ;

  /**
   * A of the remnant.
   */
  G4double massini;

  /**
   * Z of the remnant.
   */
  G4double mzini;

  /**
   * Excitation energy.
   */
  G4double exini;

  G4double pcorem, mcorem, pxrem, pyrem, pzrem;

  /**
   * Cascade n multip.
   */
  G4int mulncasc;

  /**
   * Evaporation n multip.
   */
  G4int mulnevap;

  /**
   * Total n multip.
   */
  G4int mulntot;

  /**
   * Impact parameter.
   */
  G4double bimpact;

  /**
   * Remnant Intrinsic Spin.
   */
  G4int jremn;

  /**
   * Fission 1/0=Y/N.
   */
  G4int kfis;

  /**
   * Excit energy at fis.
   */
  G4double estfis;

  /**
   * Z of fiss nucleus.
   */
  G4int izfis;

  /**
   * A of fiss nucleus.
   */
  G4int iafis;

  /**
   * Number of particles.
   */
  G4int ntrack;

  /**
   * The state of the index:
   * true = reserved
   * false = free
   */
  G4bool full[VARNTPSIZE];

  /**
   * emitted in cascade (0) or evaporation (1).
   */
  G4int itypcasc[VARNTPSIZE];

  
  /**
   * A (-1 for pions).
   */
  G4int avv[VARNTPSIZE];

  /**
   * Z
   */
  G4int zvv[VARNTPSIZE];

  /**
   * Kinetic energy.
   */
  G4double enerj[VARNTPSIZE];

  /**
   * Momentum.
   */
  G4double plab[VARNTPSIZE];

  /**
   * Theta angle.
   */
  G4double tetlab[VARNTPSIZE];

  /**
   * Phi angle.
   */
  G4double philab[VARNTPSIZE];

private:
  G4int particleIndex;
};

/**
 * Pauli blocking.
 */
class G4Paul {
public:
  G4Paul() {};
  ~G4Paul() {};
  
  /**
   *
   */
  G4double ct0,ct1,ct2,ct3,ct4,ct5,ct6,pr,pr2,xrr,xrr2;

  /**
   *
   */
  G4double cp0,cp1,cp2,cp3,cp4,cp5,cp6;
};


#endif
