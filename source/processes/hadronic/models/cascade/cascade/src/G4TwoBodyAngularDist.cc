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
// $Id$
// Author:  Michael Kelsey (SLAC)
// Date:    21 February 2013
//
// Description: Singleton class to evaluate two-body angular distribution
//		functions based on intial/final state codes.
//
// 20130307  M. Kelsey -- Add verbosity interface

#include "G4TwoBodyAngularDist.hh"
#include "G4GamP2NPipAngDst.hh"
#include "G4GamP2PPi0AngDst.hh"
#include "G4GammaNuclAngDst.hh"
#include "G4HadNElastic1AngDst.hh"
#include "G4HadNElastic2AngDst.hh"
#include "G4InuclParticleNames.hh"
#include "G4NuclNuclAngDst.hh"
#include "G4PiNInelasticAngDst.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

// Constructor and destructor

const G4TwoBodyAngularDist G4TwoBodyAngularDist::theInstance;

G4TwoBodyAngularDist::G4TwoBodyAngularDist()
  : gp_npip(new G4GamP2NPipAngDst), gp_ppi0(new G4GamP2PPi0AngDst),
    nnAngDst(new G4NuclNuclAngDst), qxAngDst(new G4PiNInelasticAngDst),
    hn1AngDst(new G4HadNElastic1AngDst), hn2AngDst(new G4HadNElastic2AngDst),
    gnAngDst(new G4GammaNuclAngDst) {;}

G4TwoBodyAngularDist::~G4TwoBodyAngularDist() {
  delete gp_npip;
  delete gp_ppi0;
  delete nnAngDst;
  delete qxAngDst;
  delete hn1AngDst;
  delete hn2AngDst;
  delete gnAngDst;
}


// Set verbosity for all generators (const-cast required)

void G4TwoBodyAngularDist::setVerboseLevel(G4int verbose) {
  const_cast<G4TwoBodyAngularDist&>(theInstance).passVerbose(verbose);
}

void G4TwoBodyAngularDist::passVerbose(G4int verbose) {
  if (gp_npip)   gp_npip->setVerboseLevel(verbose);
  if (gp_ppi0)   gp_ppi0->setVerboseLevel(verbose);
  if (nnAngDst)  nnAngDst->setVerboseLevel(verbose);
  if (qxAngDst)  qxAngDst->setVerboseLevel(verbose);
  if (hn1AngDst) hn1AngDst->setVerboseLevel(verbose);
  if (hn2AngDst) hn2AngDst->setVerboseLevel(verbose);
  if (gnAngDst)  gnAngDst->setVerboseLevel(verbose);
}


// Return appropriate distribution generator for specified interaction

const G4VTwoBodyAngDst* 
G4TwoBodyAngularDist::ChooseDist(G4int is, G4int fs, G4int kw) const {
  // gamma-nucleon -> nucleon pi0
  if ((is == gam*pro && fs == pro*pi0) ||
      (is == gam*neu && fs == neu*pi0)) {
    return gp_ppi0;
  } 

  // gamma-nucleon charge exchange
  if ((is == gam*pro && fs == neu*pip) ||
      (is == gam*neu && fs == pro*pim)) {
    return gp_npip;
  } 

    // nucleon-nucleon or hyperon-nucleon
  if (is == pro*pro || is == pro*neu || is == neu*neu ||
      is == pro*lam || is == pro*sp  || is == pro*s0  ||
      is == pro*sm  || is == pro*xi0 || is == pro*xim ||
      is == pro*om  ||
      is == neu*lam || is == neu*sp  || is == neu*s0  ||
      is == neu*sm  || is == neu*xi0 || is == neu*xim ||
      is == neu*om) {
    return nnAngDst;
  }

  // gamma p -> pi+ n, gamma p -> pi0 p, gamma p -> K Y (and isospin variants)
  if (kw == 2 && (is == pro*gam || is == neu*gam)) {
    return gnAngDst;
  }

  // pion-nucleon charge/strangeness exchange
  if (kw == 2) {
    return qxAngDst;
  }

  // pi+p, pi0p, gamma p, k+p, k0bp, pi-n, pi0n, gamma n, k-n, or k0n
  if (is == pro*pip || is == pro*pi0 || is == pro*gam ||
      is == pro*kpl || is == pro*k0b ||
      is == neu*pim || is == neu*pi0 || is == neu*gam ||
      is == neu*kmi || is == neu*k0) {
    return hn1AngDst;
  }

  // pi-p, pi+n, k-p, k0bn, k+n, or k0p
  if (is == pro*pim || is == pro*kmi || is == pro*k0 ||
      is == neu*pip || is == neu*kpl || is == neu*k0b) {
    return hn2AngDst;
  }

  // Invalid interaction
  return 0;
}