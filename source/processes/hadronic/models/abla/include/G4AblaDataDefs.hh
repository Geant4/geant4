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
// $Id: G4AblaDataDefs.hh,v 1.1 2008-02-27 18:31:11 miheikki Exp $ 
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
};

class G4Ald {
public:
  /**
   * 
   */
  G4Ald() {};
  ~G4Ald() {};
  
  G4double av,as,ak,optafan;
};

class G4Ablamain {
public:
  G4Ablamain() {};
  ~G4Ablamain() {};
  
  G4double ap,zp,at,zt,eap,beta,bmaxnuc,crtot,crnuc,r_0, r_p,r_t,pi,bfpro,snpro,sppro,shell;
  G4int imax, inum;
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
  G4Fiss() {};
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
  G4Opt() {};
  ~G4Opt() {};
  
  G4int optemd,optcha;
  G4double eefac;                                  
};

#define EENUCSIZE 2002
#define XHESIZE 50
class G4Eenuc {
public:
  G4Eenuc() {};
  ~G4Eenuc() {};
  
  G4double she[EENUCSIZE],xhe[XHESIZE][EENUCSIZE];                                            
};

#define EMDPARSIZE 1000
/**
 * Energies widths and cross sections for em excitation.
 */

class G4Emdpar {

public:
  G4Emdpar() {};
  ~G4Emdpar() {};
  
  G4double egdr,egqr,fwhmgdr,fwhmgqr,cremde1,cremde2;                  
  G4double ae1[EMDPARSIZE],be1[EMDPARSIZE],ce1[EMDPARSIZE],ae2[EMDPARSIZE];      
  G4double be2[EMDPARSIZE],ce2[EMDPARSIZE],sre1[EMDPARSIZE],sre2[EMDPARSIZE];
  G4double xre1[EMDPARSIZE],xre2[EMDPARSIZE],ds1,ds2;                             
};

//#define VOLANTSIZE 200
#define VOLANTSIZE 2000
/**
 * Evaporation and fission output data.
 */

class G4Volant {
  
public:
  G4Volant() {};
  ~G4Volant() {};

  void dump()
  {
    G4cout <<"i \t ACV \t ZPCV \t PCV" << G4endl; 
    for(G4int i = 0; i <= iv; i++) {
      G4cout << "volant" << i << "\t" << acv[i] << " \t " << zpcv[i] << " \t " << pcv[i] << G4endl;
    }
  }

  G4double acv[VOLANTSIZE],zpcv[VOLANTSIZE],pcv[VOLANTSIZE],xcv[VOLANTSIZE];
  G4double ycv[VOLANTSIZE],zcv[VOLANTSIZE];
  G4int iv; 
};

#endif
