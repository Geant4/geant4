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
// $Id: G4AblaDataDefs.hh,v 1.1 2007-05-25 05:32:03 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Data structures needed by ABLA evaporation code.

#ifndef G4AblaDataDefs_hh
#define G4AblaDataDefs_hh 1

// ABLA

typedef struct {
	G4int ii;
} G4Nevent;

// ABLA
#define PACESIZEROWS 500
#define PACESIZECOLS 500
typedef struct {
	G4double dm[PACESIZEROWS][PACESIZECOLS];
} G4Pace;

#define EC2SUBROWS 154
#define EC2SUBCOLS 99
typedef struct {
	G4double ecnz[EC2SUBROWS][EC2SUBCOLS]; 
} G4Ec2sub;

typedef struct {
	G4double av,as,ak,optafan;
} G4Ald;

typedef struct {
	G4double ap,zp,at,zt,eap,beta,bmaxnuc,crtot,crnuc,r_0, r_p,r_t,pi,bfpro,snpro,sppro,shell;
	G4int imax, inum;
} G4Ablamain;

#define ECLDROWS 154
#define ECLDCOLS 99
typedef struct {
	G4double ecgnz[ECLDROWS][ECLDCOLS];
	G4double ecfnz[ECLDROWS][ECLDCOLS];
	G4double vgsld[ECLDROWS][ECLDCOLS]; 
	G4double alpha[ECLDROWS][ECLDCOLS];
} G4Ecld;

typedef struct {
	G4double akap,bet,homega,koeff,ifis;
	G4int optshp, optxfis,optles,optcol;
} G4Fiss;

#define FBROWS 101
#define FBCOLS 161
typedef struct {
	G4double efa[FBROWS][FBCOLS];
} G4Fb;

typedef struct {
  G4int optemd,optcha;
  G4double eefac;                                  
} G4Opt;

#define EENUCSIZE 2002
#define XHESIZE 50
typedef struct {
  G4double she[EENUCSIZE],xhe[XHESIZE][EENUCSIZE];                                            
} G4Eenuc;

#define EMDPARSIZE 1000
typedef struct {
  G4double egdr,egqr,fwhmgdr,fwhmgqr,cremde1,cremde2;                  
  G4double ae1[EMDPARSIZE],be1[EMDPARSIZE],ce1[EMDPARSIZE],ae2[EMDPARSIZE];      
  G4double be2[EMDPARSIZE],ce2[EMDPARSIZE],sre1[EMDPARSIZE],sre2[EMDPARSIZE];
  G4double xre1[EMDPARSIZE],xre2[EMDPARSIZE],ds1,ds2;                             
} G4Emdpar;

#define VOLANTSIZE 200
typedef struct {
  G4double acv[VOLANTSIZE],zpcv[VOLANTSIZE],pcv[VOLANTSIZE],xcv[VOLANTSIZE];
  G4double ycv[VOLANTSIZE],zcv[VOLANTSIZE];
  G4int iv; 
} G4Volant;

#endif
