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
// G4CascadeParameters.cc
// Encapsulate all user-configurable parameters with associated envvars
//
// 20120912  M. Kelsey -- Add interface to support UI commands

#include "G4CascadeParameters.hh"
#include "G4CascadeParamMessenger.hh"
#include <stdlib.h>
#include <iostream>
using std::endl;

G4CascadeParameters* G4CascadeParameters::fpInstance = 0;

// Singleton accessor

G4CascadeParameters* G4CascadeParameters::Instance() {
  if (!fpInstance)
  {
    fpInstance = new G4CascadeParameters;
  }
  return fpInstance;
}


// Constructor initializes everything once

#define OLD_RADIUS_UNITS (3.3836/1.2)		// Used with NucModel params

G4CascadeParameters::G4CascadeParameters()
  : G4CASCADE_VERBOSE(getenv("G4CASCADE_VERBOSE")),
    G4CASCADE_USE_PRECOMPOUND(getenv("G4CASCADE_USE_PRECOMPOUND")),
    G4CASCADE_DO_COALESCENCE(getenv("G4CASCADE_DO_COALESCENCE")),
    G4CASCADE_RANDOM_FILE(getenv("G4CASCADE_RANDOM_FILE")),
    G4NUCMODEL_USE_BEST(getenv("G4NUCMODEL_USE_BEST")),
    G4NUCMODEL_RAD_2PAR(getenv("G4NUCMODEL_RAD_2PAR")),
    G4NUCMODEL_RAD_SCALE(getenv("G4NUCMODEL_RAD_SCALE")),
    G4NUCMODEL_RAD_SMALL(getenv("G4NUCMODEL_RAD_SMALL")),
    G4NUCMODEL_RAD_ALPHA(getenv("G4NUCMODEL_RAD_ALPHA")),
    G4NUCMODEL_RAD_TRAILING(getenv("G4NUCMODEL_RAD_TRAILING")),
    G4NUCMODEL_FERMI_SCALE(getenv("G4NUCMODEL_FERMI_SCALE")),
    G4NUCMODEL_XSEC_SCALE(getenv("G4NUCMODEL_XSEC_SCALE")),
    G4NUCMODEL_GAMMAQD(getenv("G4NUCMODEL_GAMMAQD")),
    DPMAX_2CLUSTER(getenv("DPMAX_2CLUSTER")),
    DPMAX_3CLUSTER(getenv("DPMAX_3CLUSTER")),
    DPMAX_4CLUSTER(getenv("DPMAX_4CLUSTER")),
    messenger(0) {
  messenger = new G4CascadeParamMessenger(this);
  Initialize();
}

void G4CascadeParameters::Initialize() {
  VERBOSE_LEVEL = (G4CASCADE_VERBOSE ? atoi(G4CASCADE_VERBOSE) : 0);
  USE_PRECOMPOUND = (0!=G4CASCADE_USE_PRECOMPOUND);
  DO_COALESCENCE = (0!=G4CASCADE_DO_COALESCENCE);
  RANDOM_FILE = (G4CASCADE_RANDOM_FILE ? G4CASCADE_RANDOM_FILE : "");
  BEST_PAR = (0!=G4NUCMODEL_USE_BEST);
  TWOPARAM_RADIUS = (0!=G4NUCMODEL_RAD_2PAR);
  RADIUS_SCALE = (G4NUCMODEL_RAD_SCALE ? strtod(G4NUCMODEL_RAD_SCALE,0)
		  : (BEST_PAR?1.0:OLD_RADIUS_UNITS));
  RADIUS_SMALL = ((G4NUCMODEL_RAD_SMALL ? strtod(G4NUCMODEL_RAD_SMALL,0)
		   : (BEST_PAR?1.992:(8.0/OLD_RADIUS_UNITS))) * RADIUS_SCALE);
  RADIUS_ALPHA = (G4NUCMODEL_RAD_ALPHA ? strtod(G4NUCMODEL_RAD_ALPHA,0)
		  : (BEST_PAR?0.84:0.70));
  RADIUS_TRAILING = ((G4NUCMODEL_RAD_TRAILING ? strtod(G4NUCMODEL_RAD_TRAILING,0)
		      : (BEST_PAR?0.70:0.0)) * RADIUS_SCALE);
  FERMI_SCALE = ((G4NUCMODEL_FERMI_SCALE ? strtod(G4NUCMODEL_FERMI_SCALE,0)
		  : (BEST_PAR?0.685:(1.932/OLD_RADIUS_UNITS))) * RADIUS_SCALE);
  XSEC_SCALE = (G4NUCMODEL_XSEC_SCALE ? strtod(G4NUCMODEL_XSEC_SCALE,0)
		: (BEST_PAR?0.1:1.0) );
  GAMMAQD_SCALE = (G4NUCMODEL_GAMMAQD?strtod(G4NUCMODEL_GAMMAQD,0):1.);
  DPMAX_DOUBLET = (DPMAX_2CLUSTER ? strtod(DPMAX_2CLUSTER,0) : 0.090);
  DPMAX_TRIPLET = (DPMAX_3CLUSTER ? strtod(DPMAX_3CLUSTER,0) : 0.108);
  DPMAX_ALPHA = (DPMAX_4CLUSTER ? strtod(DPMAX_4CLUSTER,0) : 0.115);
}

G4CascadeParameters::~G4CascadeParameters() {
  delete messenger;
}


// Report any non-default parameters (used by G4CascadeInterface)

void G4CascadeParameters::DumpConfig(std::ostream& os) const {
  if (G4CASCADE_VERBOSE)
    os << "G4CASCADE_VERBOSE = " << G4CASCADE_VERBOSE << endl;
  if (G4CASCADE_USE_PRECOMPOUND)
    os << "G4CASCADE_USE_PRECOMPOUND = " << G4CASCADE_USE_PRECOMPOUND << endl;
  if (G4CASCADE_DO_COALESCENCE)
    os << "G4CASCADE_DO_COALESCENCE = " << G4CASCADE_DO_COALESCENCE << endl;
  if (G4CASCADE_RANDOM_FILE)
    os << "G4CASCADE_RANDOM_FILE = " << G4CASCADE_RANDOM_FILE << endl;
  if (G4NUCMODEL_USE_BEST)
    os << "G4NUCMODEL_USE_BEST = " << G4NUCMODEL_USE_BEST << endl;
  if (G4NUCMODEL_RAD_2PAR)
    os << "G4NUCMODEL_RAD_2PAR = " << G4NUCMODEL_RAD_2PAR << endl;
  if (G4NUCMODEL_RAD_SCALE)
    os << "G4NUCMODEL_RAD_SCALE = " << G4NUCMODEL_RAD_SCALE << endl;
  if (G4NUCMODEL_RAD_SMALL)
    os << "G4NUCMODEL_RAD_SMALL = " << G4NUCMODEL_RAD_SMALL << endl;
  if (G4NUCMODEL_RAD_ALPHA)
    os << "G4NUCMODEL_RAD_ALPHA = " << G4NUCMODEL_RAD_ALPHA << endl;
  if (G4NUCMODEL_RAD_TRAILING)
    os << "G4NUCMODEL_RAD_TRAILING = " << G4NUCMODEL_RAD_TRAILING << endl;
  if (G4NUCMODEL_FERMI_SCALE)
    os << "G4NUCMODEL_FERMI_SCALE = " << G4NUCMODEL_FERMI_SCALE << endl;
  if (G4NUCMODEL_XSEC_SCALE)
    os << "G4NUCMODEL_XSEC_SCALE = " << G4NUCMODEL_XSEC_SCALE << endl;
  if (G4NUCMODEL_GAMMAQD)
    os << "G4NUCMODEL_GAMMAQD = " << G4NUCMODEL_GAMMAQD << endl;
  if (DPMAX_2CLUSTER)
    os << "DPMAX_2CLUSTER = " << DPMAX_2CLUSTER << endl;
  if (DPMAX_3CLUSTER)
    os << "DPMAX_3CLUSTER = " << DPMAX_3CLUSTER << endl;
  if (DPMAX_4CLUSTER)
    os << "DPMAX_4CLUSTER = " << DPMAX_4CLUSTER << endl;
}
