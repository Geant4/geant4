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
// ABLAXX statistical de-excitation model
// Jose Luis Rodriguez, GSI (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (initial translation of ablav3p)
// Aleksandra Kelic, GSI (ABLA07 code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//
#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"

// Data structures needed by ABLA evaporation code.

#ifndef G4AblaDataDefs_hh
#define G4AblaDataDefs_hh 1

#ifdef ABLAXX_IN_GEANT4_MODE
#include "globals.hh"
#else
#include "G4INCLGeant4Compat.hh"
#endif

#include <cmath>

// ABLA

class G4Nevent {
public:
  G4Nevent() {};
  ~G4Nevent() {};
  
  G4int ii;
};

// ABLA
#define PACESIZEROWS 500
#define PACESIZECOLS 500
/**
 * Masses.
 */

class G4Pace {

public:
  G4Pace() {};

  ~G4Pace() {};
  
  G4double dm[PACESIZEROWS][PACESIZECOLS];
};

#define MASSIZEROWS 154
#define MASSIZECOLS 13

class G4Mexp {

public:
  G4Mexp() {};

  ~G4Mexp() {};
  
  G4double massexp[MASSIZEROWS][MASSIZECOLS];
  G4double bind[MASSIZEROWS][MASSIZECOLS];
  G4int mexpiop[MASSIZEROWS][MASSIZECOLS];
};

#define EC2SUBROWS 154
#define EC2SUBCOLS 99
  /**
   *
   */

class G4Ec2sub {
public:
  G4Ec2sub() {};

  ~G4Ec2sub() {};

  G4double ecnz[EC2SUBROWS][EC2SUBCOLS]; 

  /**
   * Dump the contents of the ecnz data table.
   */
  void dump() {
    for(G4int i = 0; i < EC2SUBROWS; i++) {
      for(G4int j = 0; j < EC2SUBCOLS; j++) {
	//G4cout << ecnz[i][j] << " ";
      }
      //      G4cout << G4endl;
    }
  }
};

class G4Ald {
public:
  /**
   * 
   */
  G4Ald()
    :av(0.0), as(0.0), ak(0.0), optafan(0.0)
  {};
  ~G4Ald() {};
  
  G4double av,as,ak,optafan;
};

#define ECLDROWS 154
#define ECLDCOLS 99

#define ECLDROWSbeta 251
#define ECLDCOLSbeta 137
/**
 * Shell corrections and deformations.
 */

class G4Ecld {

public:
  G4Ecld() {};
  ~G4Ecld() {};

  /**
   * Ground state shell correction frldm for a spherical ground state.
   */
  G4double ecgnz[ECLDROWS][ECLDCOLS];

  /**
   * Shell correction for the saddle point (now: == 0).
   */
  G4double ecfnz[ECLDROWS][ECLDCOLS];

  /**
   * Difference between deformed ground state and ldm value.
   */
  G4double vgsld[ECLDROWS][ECLDCOLS]; 

  /**
   * Alpha ground state deformation (this is not beta2!)       
   * beta2 = std::sqrt(5/(4pi)) * alpha 
   */
  G4double alpha[ECLDROWS][ECLDCOLS];

  /**
   * RMS function for lcp emission barriers
   */
  G4double rms[ECLDROWS][ECLDCOLS];

  /**
   * Beta2 deformations
   */
  G4double beta2[ECLDROWSbeta][ECLDCOLSbeta];

  /**
   * Beta4 deformations
   */
  G4double beta4[ECLDROWSbeta][ECLDCOLSbeta];
};

class G4Fiss {
  /**
   * Options and parameters for fission channel.
   */

public:
  G4Fiss()
    :bet(0.0), ifis(0.0), ucr(0.0), dcr(0.0), optshp(0), optxfis(0), optct(0), optcol(0), 
     at(0), zt(0)
  {};
  ~G4Fiss() {};
  
  G4double bet,ifis,ucr,dcr;
  G4int optshp, optxfis,optct,optcol,at,zt;
};

#define FBROWS 101
#define FBCOLS 161
/**
 * Fission barriers.
 */
 
class G4Fb {
  
public:
  G4Fb() {};
  ~G4Fb() {;}
  
  //  G4double efa[FBROWS][FBCOLS];
  G4double efa[FBCOLS][FBROWS];
};

/**
 * Options
 */

class G4Opt {

public:
  G4Opt()
    :optemd(0), optcha(0), optshpimf(0), optimfallowed(0), nblan0(0)
  {};
  ~G4Opt() {};
  
  G4int optemd,optcha,optshpimf,optimfallowed,nblan0;                                 
};

#define EENUCSIZE 2002
#define XHESIZE 50
class G4Eenuc {
public:
  G4Eenuc() {
    for(G4int i = 0; i < EENUCSIZE; ++i) {
      she[i] = 0.0;
    }
    for(G4int i = 0; i < XHESIZE; ++i) {
      for(G4int j = 0; j < EENUCSIZE; ++j) {
	xhe[i][j] = 0.0;
      }
    }
  };
  ~G4Eenuc() {};
  
  G4double she[EENUCSIZE],xhe[XHESIZE][EENUCSIZE];                                            
};

//#define VOLANTSIZE 200
#define VOLANTSIZE 301
/**
 * Evaporation and fission output data.
 */

class G4Volant {
  
public:
  G4Volant()
  {
    clear();
  }

  ~G4Volant() {};

  void clear()
  {
    for(G4int i = 0; i < VOLANTSIZE; i++) {
      copied[i] = false;
      acv[i] = 0;
      zpcv[i] = 0;
      pcv[i] = 0;
      xcv[i] = 0;
      ycv[i] = 0;
      zcv[i] = 0;
      iv = 0;
    }
  }

  G4double getTotalMass()
  {
    G4double total = 0.0;
    for(G4int i = 0; i <= iv; i++) {
      total += acv[i];
    }
    return total;
  }

  void dump()
  {
/*
    G4double totA = 0.0, totZ = 0.0, totP = 0.0;
    //    G4cout <<"i \t ACV \t ZPCV \t PCV" << G4endl; 
    for(G4int i = 0; i <= iv; i++) {
      if(i == 0 && acv[i] != 0) {
	//	G4cout <<"G4Volant: Particle stored at index " << i << G4endl;
      }
      totA += acv[i];
      totZ += zpcv[i];
      totP += pcv[i];
      //      G4cout << "volant" << i << "\t" << acv[i] << " \t " << zpcv[i] << " \t " << pcv[i] << G4endl;
    }
    //    G4cout <<"Particle count index (iv) = " << iv << G4endl;
    //    G4cout <<"ABLA Total: A = " << totA << " Z = " << totZ <<  " momentum = " << totP << G4endl;
*/
  }

  G4double acv[VOLANTSIZE],zpcv[VOLANTSIZE],pcv[VOLANTSIZE],xcv[VOLANTSIZE];
  G4double ycv[VOLANTSIZE],zcv[VOLANTSIZE];
  G4bool copied[VOLANTSIZE];
  G4int iv; 
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
      svv[i] = 0;
      enerj[i] = 0.0;
      pxlab[i] = 0.0;
      pylab[i] = 0.0;
      pzlab[i] = 0.0;
      full[i] = false;
    }
  }

  /**
   * Add a particle to the INCL/ABLA final output.
   */
  void addParticle(G4double A, G4double Z, G4double E, G4double P, G4double theta, G4double phi) {
    if(full[particleIndex]) {
      //      G4cout <<"A = " << Z << " Z = " << Z << G4endl;
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
  void dump()
  {
/*
    G4int nProton = 0, nNeutron = 0;
    G4int nPiPlus = 0, nPiZero = 0, nPiMinus = 0;
    G4int nH2 = 0, nHe3 = 0, nAlpha = 0;
    G4int nGamma=0;
    G4int nFragments = 0;
    G4int nParticles = 0;
    for(G4int i = 0; i < ntrack; i++) {
      nParticles++;
      if(avv[i] ==  1 && zvv[i] ==  1) nProton++;  // Count multiplicities
      if(avv[i] ==  1 && zvv[i] ==  0) nNeutron++;
      if(avv[i] ==  0 && zvv[i] ==  0) nGamma++;
      if(avv[i] == -1 && zvv[i] ==  1) nPiPlus++;
      if(avv[i] == -1 && zvv[i] ==  0) nPiZero++;
      if(avv[i] == -1 && zvv[i] == -1) nPiMinus++;
      if(avv[i] ==  2 && zvv[i] ==  1) nH2++;
      if(avv[i] ==  3 && zvv[i] ==  2) nHe3++;
      if(avv[i] ==  4 && zvv[i] ==  2) nAlpha++;
      if(                zvv[i] >   2) nFragments++;
    }
*/
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
   * S (-1 for lambda_0).
   */
  G4int svv[VARNTPSIZE];

  /**
   * Kinetic energy.
   */
  G4double enerj[VARNTPSIZE];

  /**
   * Momentum.
   */
  G4double plab[VARNTPSIZE];
  G4double pxlab[VARNTPSIZE];
  G4double pylab[VARNTPSIZE];
  G4double pzlab[VARNTPSIZE];

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

#endif
