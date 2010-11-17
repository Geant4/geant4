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
// $Id: G4InclDataDefs.hh,v 1.12 2010-11-17 20:19:09 kaitanie Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// All data structures needed by INCL4 are defined here.

#ifndef InclDataDefs_hh
#define InclDataDefs_hh 1

#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"

class G4InclFermi {
public:
  G4InclFermi() {
    G4double hc = 197.328;
    G4double fmp = 938.2796;
    pf=1.37*hc;
    pf2=pf*pf;
    tf=std::sqrt(pf*pf+fmp*fmp)-fmp;
  };
  ~G4InclFermi() {};

  G4double tf,pf,pf2;
};

#define max_a_proj 61

/**
 * (eps_c,p1_s,p2_s,p3_s,eps_c used to store the kinematics of
 * nucleons for composit projectiles before entering the potential)
 */
class G4QuadvectProjo {
public:
  G4QuadvectProjo() {
    for(G4int i = 0; i < max_a_proj; ++i) {
      eps_c[i] = 0.0;
      t_c[i] = 0.0;
      p3_c[i] = 0.0;
      p1_s[i] = 0.0;
      p2_s[i] = 0.0;
      p3_s[i] = 0.0;
    }
  };

  ~G4QuadvectProjo() {};

  G4double eps_c[max_a_proj],p3_c[max_a_proj],
    p1_s[max_a_proj],p2_s[max_a_proj],p3_s[max_a_proj],
    t_c[max_a_proj];
};

class G4VBe {
public:
  G4VBe()
    :ia_be(0), iz_be(0),
     rms_be(0.0), pms_be(0.0), bind_be(0.0)
  { };

  ~G4VBe() {};

  G4int ia_be, iz_be;
  G4double rms_be, pms_be, bind_be;
};

/**
 * Projectile spectator
 */
class G4InclProjSpect {
public:
  G4InclProjSpect() {
    //    G4cout <<"Projectile spectator data structure created!" << G4endl;
    clear();
  };
  ~G4InclProjSpect() {};

  void clear() {
    for(G4int i = 0; i < 21; i++) tab[i] = 0.0;
    for(G4int i = 0; i < 61; i++) n_projspec[i] = 0;
    a_projspec = 0;
    z_projspec = 0;
    t_projspec = 0.0;
    ex_projspec = 0.0;
    p1_projspec = 0.0;
    p2_projspec = 0.0;
    p3_projspec = 0.0;
    m_projspec = 0.0;
  };

  G4double tab[21];
  G4int n_projspec[61];
  G4int a_projspec,z_projspec;
  G4double ex_projspec,t_projspec, p1_projspec, p2_projspec, p3_projspec, m_projspec;
};

#define IGRAINESIZE 19
/**
 * Random seeds used by internal random number generators.
 * @see G4Incl::standardRandom
 * @see G4Incl::gaussianRandom
 */
class G4Hazard{
public:
  G4Hazard() {
    ial = 0;
    for(G4int i = 0; i < IGRAINESIZE; ++i) {
      igraine[i] = 0;
    }
  };

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
  G4Mat() {
    nbmat = 0;

    for(G4int i = 0; i < MATSIZE; ++i) {
      zmat[i] = 0;
      amat[i] = 0;
    }

    for(G4int i = 0; i < MATGEOSIZE; ++i) {
      for(G4int j = 0; j < MATSIZE; ++j) {
	bmax_geo[i][j] = 0;
      }
    }
};

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
  G4LightGausNuc() {
    for(G4int i = 0; i < LGNSIZE; ++i) {
      rms1t[i] = 0.0;
      pf1t[i] = 0.0;
      pfln[i] = 0.0;
      tfln[i] = 0.0;
      vnuc[i] = 0.0;
    }
  };

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
  G4LightNuc() {
    for(G4int i = 0; i < LNSIZE; ++i) {
      r[i] = 0.0;
      a[i] = 0.0;
    }
  };

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
  G4Saxw() {
    for(G4int i = 0; i < SAXWROWS; ++i) {
      for(G4int j = 0; j < SAXWCOLS; ++j) {
	x[i][j] = 0.0;
	y[i][j] = 0.0;
	s[i][j] = 0.0;
      }
    }
    imat = 0; n = 0; k = 0;
  };

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
  G4Ws() {
    fneck = 0.0;
    r0 = 0.0;
    adif = 0.0;
    rmaxws = 0.0;
    drws = 0.0;
    nosurf = 0.0;
    xfoisa = 0.0;
    bmax = 0.0;
    npaulstr = 0.0;
  };

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

  G4double fneck;
};

#define DTONSIZE 13
/**
 * Random seeds used by internal random number generators.
 * @see G4Incl::standardRandom
 * @see G4Incl::gaussianRandom
 */
class G4Dton {
public:
  G4Dton() {
    fn = 0.0;
    for(G4int i = 0; i < DTONSIZE; ++i) {
      c[i] = 0.0;
      d[i] = 0.0;
    }
  };

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
  G4Spl2() {
    for(G4int i = 0; i < SPL2SIZE; ++i) {
      x[i] = 0.0; y[i] = 0.0;
      a[i] = 0.0; b[i] = 0.0; c[i] = 0.0;
    }
    n = 0;
  };

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
  G4Bl1() {
    ta = 0.0;
    for(G4int i = 0; i < BL1SIZE; ++i) {
      p1[i] = 0.0; p2[i] = 0.0; p3[i] = 0.0; eps[i] = 0.0;
      ind1[i] = 0; ind2[i] = 0;
    }
  };

  ~G4Bl1() {};
  
  G4double p1[BL1SIZE],p2[BL1SIZE],p3[BL1SIZE];
  G4double eps[BL1SIZE];
  G4int ind1[BL1SIZE],ind2[BL1SIZE];
  G4double ta;

  void dump(G4int numberOfParticles) {
    static G4int dumpNumber = 0;
    G4cout <<"Dump number" << dumpNumber << " of particle 4-momenta (G4Bl1):" << G4endl;
    G4cout <<"ta = " << ta << G4endl;
    for(G4int i = 0; i < numberOfParticles; i++) {
      G4cout <<"i = " << i << "   p1 = " << p1[i] << "   p2 = " << p2[i] << "   p3 = " << p3[i] << "   eps = " << eps[i] << G4endl;
    }
    dumpNumber++;
  }
};

#define BL2SIZE 19900
/**
 * 
 */
class G4Bl2 {
public:
  G4Bl2() {
    k = 0;
    for(G4int i = 0; i < BL2SIZE; ++i) {
      crois[i] = 0.0;
      ind[i] = 0;
      jnd[i] = 0;
    }    
  };

  ~G4Bl2() {};
  
  void dump() {
    G4cout <<"Avatars: (number of avatars = " << k << ")" << G4endl;
    for(G4int i = 0; i <= k; i++) {
      G4cout <<"i = " << i << G4endl;
      G4cout <<"crois[" << i << "] = " << crois[i] << G4endl;
      G4cout <<"ind[" << i << "] = " << ind[i] << G4endl;
      G4cout <<"jnd[" << i << "] = " << jnd[i] << G4endl;
    }
  }

  /**
   * 
   */
  G4double crois[BL2SIZE];

  /**
   *
   */
  G4int k;

  /**
   *
   */
  G4int ind[BL2SIZE];

  /**
   *
   */
  G4int jnd[BL2SIZE];
};

//#define BL3SIZE 300
#define BL3SIZE 3000
/**
 *
 */
class G4Bl3 {
public:
  G4Bl3() {
    r1 = 0.0; r2 = 0.0;
    ia1 = 0; ia2 = 0;
    rab2 = 0.0;

    for(G4int i = 0; i < BL3SIZE; ++i) {
      x1[i] = 0.0; x2[i] = 0.0; x3[i] = 0.0;
    }
  };

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

  void dump() {
    static G4int dumpNumber = 0;
    G4cout <<"Dump number" << dumpNumber << " of particle positions (G4Bl3):" << G4endl;
    G4cout <<" ia1 = " << ia1 << G4endl;
    G4cout <<" ia2 = " << ia2 << G4endl;
    G4cout <<" rab2 = " << rab2 << G4endl;
    for(G4int i = 0; i <= (ia1 + ia2); i++) {
      G4cout <<"i = " << i << "   x1 = " << x1[i] << "   x2 = " << x2[i] << "   x3 = " << x3[i] << G4endl;
    }
    dumpNumber++;
  }
};

/**
 * G4Bl4
 */
class G4Bl4 {
public:
  G4Bl4()
    :tmax5(0.0)
  {};

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
  G4Bl5() {
    for(G4int i = 0; i < BL5SIZE; ++i) {
      tlg[i] = 0.0;
      nesc[i] = 0;
    }
  };

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
  G4Bl6()
    :xx10(0.0), isa(0.0)
  {};

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
  G4Bl8()
    :rathr(0.0), ramass(0.0)
  {};

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
    for(G4int i = 0; i < BL9SIZE; ++i) {
      hel[i] = 0.0;
    }
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
  G4Bl10()
    :ri4(0.0), rs4(0.0), r2i(0.0), r2s(0.0), pdummy(0.0), pf(0.0)
  {};
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
  G4Kind()
    :kindf7(0)
  {};
  ~G4Kind() {};

  /**
   * kindf7
   */
  G4int kindf7;
};

/**
 * Projectile parameters.
 */
class G4Bev {
public:
  /**
   * Initialize all variables to zero.
   */
  G4Bev() {
    ia_be = 0;
    iz_be = 0;
    rms_be = 0.0;
    pms_be = 0.0;
    bind_be = 0.0;
  };
  ~G4Bev() {};

  /**
   * Mass number.
   */
  G4int ia_be;

  /**
   * Charge number.
   */
  G4int iz_be;

  /**
   * rms
   */
  G4double rms_be;

  /**
   * pms
   */
  G4double pms_be;

  /**
   * bind
   */
  G4double bind_be;
};

#define VARSIZE 3
#define VAEPSSIZE 250
#define VAAVM 1000
/**
 * Extra information on collisions between nucleons.
 */
class G4VarAvat {
public:
  G4VarAvat() {
    kveux = 0;
    bavat = 0.0;
    nopartavat = 0; ncolavat = 0;
    nb_avat = 0;

    for(G4int i = 0; i < VARSIZE; ++i) {
      r1_in[i] = 0.0;
      r1_first_avat[i] = 0.0;
    }

    for(G4int i = 0; i < VAEPSSIZE; ++i) {
      epsd[i] = 0.0;
      eps2[i] = 0.0;
      eps4[i] = 0.0;
      eps6[i] = 0.0;
      epsf[i] = 0.0;
    }

    for(G4int i = 0; i < VAAVM; ++i) {
      timeavat[i] = 0.0;
      l1avat[i] = 0.0;
      l2avat[i] = 0.0;
      jpartl1[i] = 0.0;
      jpartl2[i] = 0.0;
      del1avat[i] = 0.0;
      del2avat[i] = 0.0;
      energyavat[i] = 0.0;
      bloc_paul[i] = 0.0;
      bloc_cdpp[i] = 0.0;
      go_out[i] = 0.0;
    }
  };

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

#define VARNTPSIZE 301
class G4VarNtp {
public:
  G4VarNtp() {
    clear();
  };

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
    masp = 0.0; mzsp = 0.0; exsp = 0.0; mrem = 0.0;
    // To be deleted?
    spectatorA = 0;
    spectatorZ = 0;
    spectatorEx = 0.0;
    spectatorM = 0.0;
    spectatorT = 0.0;
    spectatorP1 = 0.0;
    spectatorP2 = 0.0;
    spectatorP3 = 0.0;
    massini = 0;
    mzini = 0;
    exini = 0;
    pcorem = 0;
    mcorem = 0;
    pxrem = 0;
    pyrem = 0;
    pzrem = 0;
    erecrem = 0;
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
    needsFermiBreakup = false;
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

  /**
   * Add a particle to the INCL/ABLA final output.
   */
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
   * Projectile spectator A, Z, Eex;
   */
  G4double masp, mzsp, exsp, mrem;

  /**
   * Spectator nucleus mass number for light ion projectile support.
   */
  G4int spectatorA;

  /**
   * Spectator nucleus charge number for light ion projectile support.
   */
  G4int spectatorZ;

  /**
   * Spectator nucleus excitation energy for light ion projectile support.
   */
  G4double spectatorEx;

  /**
   * Spectator nucleus mass.
   */
  G4double spectatorM;

  /**
   * Spectator nucleus kinetic energy.
   */
  G4double spectatorT;

  /**
   * Spectator nucleus momentum x-component.
   */
  G4double spectatorP1;

  /**
   * Spectator nucleus momentum y-component.
   */
  G4double spectatorP2;

  /**
   * Spectator nucleus momentum z-component.
   */
  G4double spectatorP3;

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

  G4double pcorem, mcorem, pxrem, pyrem, pzrem, erecrem;

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
   * Does this nucleus require Fermi break-up treatment? Only
   * applicable when used together with Geant4.
   * true = do fermi break-up (and skip ABLA part)
   * false = use ABLA
   */
  G4bool needsFermiBreakup;

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
  G4Paul()
    :ct0(0.0), ct1(0.0), ct2(0.0), ct3(0.0), ct4(0.0), ct5(0.0), ct6(0.0),
     pr(0.0), pr2(0.0), xrr(0.0), xrr2(0.0),
     cp0(0.0), cp1(0.0), cp2(0.0), cp3(0.0), cp4(0.0), cp5(0.0), cp6(0.0)
  {};

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
