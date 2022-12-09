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
// Author:  Michael Kelsey (SLAC)
// Date:    21 February 2013
//
// Description: Singleton class to evaluate two-body angular distribution
//		functions based on initial/final state codes.
//
// 20130307  M. Kelsey -- Add verbosity interface
// 20130422  M. Kelsey -- Add three-body distributions, for temporary use
// 20130619  Change singleton instance to be thread-local, to avoid collisions.
// 20141121  Use G4AutoDelete to avoid end-of-thread memory leaks

#include "G4TwoBodyAngularDist.hh"
#include "G4AutoDelete.hh"
#include "G4GamP2NPipAngDst.hh"
#include "G4GamP2PPi0AngDst.hh"
#include "G4GammaNuclAngDst.hh"
#include "G4PP2PPAngDst.hh"
#include "G4NP2NPAngDst.hh"
#include "G4Pi0P2Pi0PAngDst.hh"
#include "G4PimP2Pi0NAngDst.hh"
#include "G4PimP2PimPAngDst.hh"
#include "G4PipP2PipPAngDst.hh"
#include "G4HadNElastic1AngDst.hh"
#include "G4HadNElastic2AngDst.hh"
#include "G4InuclParticleNames.hh"
#include "G4NuclNuclAngDst.hh"
#include "G4HadNucl3BodyAngDst.hh"
#include "G4NuclNucl3BodyAngDst.hh"
#include "G4PiNInelasticAngDst.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

// Singleton is created at first invocation

G4ThreadLocal G4TwoBodyAngularDist* G4TwoBodyAngularDist::theInstance = 0;

const G4TwoBodyAngularDist* G4TwoBodyAngularDist::GetInstance() {
  if (!theInstance) {
    theInstance = new G4TwoBodyAngularDist;
    G4AutoDelete::Register(theInstance);
  }

  return theInstance;
}

// Constructor and destructor

G4TwoBodyAngularDist::G4TwoBodyAngularDist()
  : gp_npip(new G4GamP2NPipAngDst), gp_ppi0(new G4GamP2PPi0AngDst),
    ppAngDst(new G4PP2PPAngDst), npAngDst(new G4NP2NPAngDst),
    nnAngDst(new G4NuclNuclAngDst), pi0pAngDst(new G4Pi0P2Pi0PAngDst),
    pipCXAngDst(new G4PimP2Pi0NAngDst), pimpAngDst(new G4PimP2PimPAngDst),
    pippAngDst(new G4PipP2PipPAngDst), qxAngDst(new G4PiNInelasticAngDst),
    hn1AngDst(new G4HadNElastic1AngDst), hn2AngDst(new G4HadNElastic2AngDst),
    gnAngDst(new G4GammaNuclAngDst), hn3BodyDst(new G4HadNucl3BodyAngDst),
    nn3BodyDst(new G4NuclNucl3BodyAngDst)
{;}

G4TwoBodyAngularDist::~G4TwoBodyAngularDist() {
  delete gp_npip;
  delete gp_ppi0;
  delete ppAngDst;
  delete nnAngDst;
  delete pi0pAngDst;
  delete pipCXAngDst;
  delete pimpAngDst;
  delete pippAngDst;
  delete qxAngDst;
  delete hn1AngDst;
  delete hn2AngDst;
  delete gnAngDst;
  delete npAngDst;
  delete hn3BodyDst;
  delete nn3BodyDst;
}


// Set verbosity for all generators (const-cast required)

void G4TwoBodyAngularDist::setVerboseLevel(G4int verbose) {
  const_cast<G4TwoBodyAngularDist*>(GetInstance())->passVerbose(verbose);
}

void G4TwoBodyAngularDist::passVerbose(G4int verbose) {
  if (gp_npip)   gp_npip->setVerboseLevel(verbose);
  if (gp_ppi0)   gp_ppi0->setVerboseLevel(verbose);
  if (ppAngDst)  ppAngDst->setVerboseLevel(verbose);
  if (nnAngDst)  nnAngDst->setVerboseLevel(verbose);
  if (pi0pAngDst) pi0pAngDst->setVerboseLevel(verbose);
  if (pipCXAngDst) pipCXAngDst->setVerboseLevel(verbose);
  if (pimpAngDst) pimpAngDst->setVerboseLevel(verbose);
  if (pippAngDst) pippAngDst->setVerboseLevel(verbose);
  if (qxAngDst)  qxAngDst->setVerboseLevel(verbose);
  if (hn1AngDst) hn1AngDst->setVerboseLevel(verbose);
  if (hn2AngDst) hn2AngDst->setVerboseLevel(verbose);
  if (gnAngDst)  gnAngDst->setVerboseLevel(verbose);
  if (npAngDst)  npAngDst->setVerboseLevel(verbose);
  if (hn3BodyDst) hn3BodyDst->setVerboseLevel(verbose);
  if (nn3BodyDst) nn3BodyDst->setVerboseLevel(verbose);
}


// Return appropriate distribution generator for specified interaction

const G4VTwoBodyAngDst* 
G4TwoBodyAngularDist::ChooseDist(G4int is, G4int fs, G4int kw) const {
  // TEMPORARY:  Three-body distributions for hN/NN
  if (fs==0 && kw==0) {
    if (is == pro*pro || is == pro*neu || is == neu*neu) return nn3BodyDst;
    else return hn3BodyDst;
  }

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

  // pp and nn elastic
  if (is == pro*pro || is == neu*neu) return ppAngDst;

  // np and pn elastic
  if (is == pro*neu) return npAngDst;

  // pi+ p and pi- n elastic
  if ((fs == is) && (is == pip*pro || is == pim*neu) ) return pippAngDst;

  // pi- p and pi+ n elastic
  if ((fs == is) && (is == pim*pro || is == pip*neu) ) return pimpAngDst;

  // pi0 p and pi0 n elastic
  if ((fs == is) && (is == pi0*pro || is == pi0*neu) ) return pi0pAngDst;

  // pi- p -> pi0 n, pi+ n -> pi0 p, pi0 p -> pi+ n, pi0 n -> pi- p
  if ((is == pim*pro && fs == pi0*neu) || (is == pip*neu && fs == pi0*pro) ||
      (is == pi0*pro && fs == pip*neu) || (is == pi0*neu && fs == pim*pro) )
    return pipCXAngDst;

  // hyperon-nucleon
  if (is == pro*lam || is == pro*sp  || is == pro*s0  ||
      is == pro*sm  || is == pro*xi0 || is == pro*xim ||
      is == pro*om  ||
      is == neu*lam || is == neu*sp  || is == neu*s0  ||
      is == neu*sm  || is == neu*xi0 || is == neu*xim ||
      is == neu*om) {
    return nnAngDst;
  }

  // gamma p -> K Y (and isospin variants)
  if (kw == 2 && (is == pro*gam || is == neu*gam)) {
    return gnAngDst;
  }

  // pion-nucleon strangeness production
  if (kw == 2) {
    return qxAngDst;
  }

  // gamma p, k+p, k0bp, gamma n, k-n, or k0n
  if (is == pro*gam ||
      is == pro*kpl || is == pro*k0b ||
      is == neu*gam ||
      is == neu*kmi || is == neu*k0) {
    return hn1AngDst;
  }

  // k-p, k0bn, k+n, or k0p
  if (is == pro*kmi || is == pro*k0 ||
      is == neu*kpl || is == neu*k0b) {
    return hn2AngDst;
  }

  // Invalid interaction
  return 0;
}
