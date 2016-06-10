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
// $Id: G4MultiBodyMomentumDist.cc 71652 2013-06-19 17:20:45Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: Singleton class to evaluate multi-body momentum distribution
//		functions based on intial state codes and multiplicity.
//
// 20130308  Use envvar to enable/disable use of 3-body generators.
// 20130619  Change singleton instance to be thread-local, to avoid collisions.
// 20141121  Use G4AutoDelete to avoid end-of-thread memory leaks

#include "G4MultiBodyMomentumDist.hh"
#include "G4AutoDelete.hh"
#include "G4CascadeParameters.hh"
#include "G4NuclNucl3BodyMomDst.hh"
#include "G4NuclNucl4BodyMomDst.hh"
#include "G4HadNucl3BodyMomDst.hh"
#include "G4HadNucl4BodyMomDst.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;


// Singleton is created at first invocation

G4ThreadLocal G4MultiBodyMomentumDist* G4MultiBodyMomentumDist::theInstance = 0;

const G4MultiBodyMomentumDist* G4MultiBodyMomentumDist::GetInstance() {
  if (!theInstance) {
    theInstance = new G4MultiBodyMomentumDist;
    G4AutoDelete::Register(theInstance);
  }

  return theInstance;
}

// Constructor and destructor

G4MultiBodyMomentumDist::G4MultiBodyMomentumDist()
  : nn3BodyDst(new G4NuclNucl3BodyMomDst),
    nn4BodyDst(new G4NuclNucl4BodyMomDst),
    hn3BodyDst(new G4HadNucl3BodyMomDst),
    hn4BodyDst(new G4HadNucl4BodyMomDst) {;}

G4MultiBodyMomentumDist::~G4MultiBodyMomentumDist() {
  delete nn3BodyDst;
  delete nn4BodyDst;
  delete hn3BodyDst;
  delete hn4BodyDst;
}


// Set verbosity for all generators (const-cast required)

void G4MultiBodyMomentumDist::setVerboseLevel(G4int verbose) {
  const_cast<G4MultiBodyMomentumDist*>(GetInstance())->passVerbose(verbose);
}

void G4MultiBodyMomentumDist::passVerbose(G4int verbose) {
  if (nn3BodyDst) nn3BodyDst->setVerboseLevel(verbose);
  if (nn4BodyDst) nn4BodyDst->setVerboseLevel(verbose);
  if (hn3BodyDst) hn3BodyDst->setVerboseLevel(verbose);
  if (hn4BodyDst) hn4BodyDst->setVerboseLevel(verbose);
}


// Return appropriate distribution generator for specified interaction

const G4VMultiBodyMomDst* 
G4MultiBodyMomentumDist::ChooseDist(G4int is, G4int mult) const {
  if (is == pro*pro || is == pro*neu || is == neu*neu) {
    //***** REMOVED BY VLADIMIR UZHINSKY 18 JULY 2011
    if (G4CascadeParameters::use3BodyMom() && mult==3) return nn3BodyDst;
    return nn4BodyDst;
  }

  else {	// FIXME:  All other initial states use pi-N scattering
    //***** REMOVED BY VLADIMIR UZHINSKY 18 JULY 2011
    if (G4CascadeParameters::use3BodyMom() && mult==3) return hn3BodyDst;
    return hn4BodyDst;
  }

  // Invalid interaction
  return 0;
}
