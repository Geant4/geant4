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
// $Id: G4AblaDataDefs.hh,v 1.13 2010-11-13 00:08:36 kaitanie Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Data structures needed by ABLA evaporation code.

#ifndef G4AblaDataDefs_hh
#define G4AblaDataDefs_hh 1

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
	G4cout << ecnz[i][j] << " ";
      }
      G4cout << G4endl;
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
};

class G4Fiss {
  /**
   * Options and parameters for fission channel.
   */

public:
  G4Fiss()
    :akap(0.0), bet(0.0), homega(0.0), koeff(0.0), ifis(0.0),
     optshp(0), optxfis(0), optles(0), optcol(0)
  {};
  ~G4Fiss() {};
  
  G4double akap,bet,homega,koeff,ifis;
  G4int optshp, optxfis,optles,optcol;
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
    :optemd(0), optcha(0), eefac(0.0)
  {};
  ~G4Opt() {};
  
  G4int optemd,optcha;
  G4double eefac;                                  
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
    G4double totA = 0.0, totZ = 0.0, totP = 0.0;
    G4cout <<"i \t ACV \t ZPCV \t PCV" << G4endl; 
    for(G4int i = 0; i <= iv; i++) {
      if(i == 0 && acv[i] != 0) {
	G4cout <<"G4Volant: Particle stored at index " << i << G4endl;
      }
      totA += acv[i];
      totZ += zpcv[i];
      totP += pcv[i];
      G4cout << "volant" << i << "\t" << acv[i] << " \t " << zpcv[i] << " \t " << pcv[i] << G4endl;
    }
    G4cout <<"Particle count index (iv) = " << iv << G4endl;
    G4cout <<"ABLA Total: A = " << totA << " Z = " << totZ <<  " momentum = " << totP << G4endl;
  }

  G4double acv[VOLANTSIZE],zpcv[VOLANTSIZE],pcv[VOLANTSIZE],xcv[VOLANTSIZE];
  G4double ycv[VOLANTSIZE],zcv[VOLANTSIZE];
  G4bool copied[VOLANTSIZE];
  G4int iv; 
};

#endif
