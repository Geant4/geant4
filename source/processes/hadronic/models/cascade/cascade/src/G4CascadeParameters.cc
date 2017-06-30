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
// $Id: G4CascadeParameters.cc 72016 2013-07-03 16:24:15Z mkelsey $
// Encapsulate all user-configurable parameters with associated envvars
//
// 20120912  M. Kelsey -- Add interface to support UI commands
// 20130304  M. Kelsey -- Add flag to collect and display cascade structure
// 20130308  M. Kelsey -- Add flag to use separate 3-body momentum generators
// 20130421  M. Kelsey -- Add flag for CHECK_ECONS, replacing #ifdef's
// 20130702  M. Kelsey -- Add flag to use N-body phase-space generator
// 20140311  G. Cosmo -- Implement standard (non-const) singleton pattern
// 20140929  M. Kelsey -- Enable some parameters as default true (must be set
//		'0' for false): PreCompound, phase-space, clustering,
//		trailing effect
// 20141111  M. Kelsey -- Revert defaults for PreCompound, phase-space,
//		and trailing effect.
// 20141121  Use G4AutoDelete to avoid end-of-thread memory leaks
// 20141211  M. Kelsey -- Change PIN_ABSORPTION flag to double, for energy cut

#include "G4CascadeParameters.hh"
#include "G4CascadeParamMessenger.hh"
#include "G4AutoDelete.hh"
#include <stdlib.h>
#include <iostream>
using std::endl;

#define OLD_RADIUS_UNITS (3.3836/1.2)		// Used with NucModel params

#include "G4HadronicDeveloperParameters.hh"
namespace { 
   G4HadronicDeveloperParameters& HDP = G4HadronicDeveloperParameters::GetInstance();
   class BERTParameters {
      public: 
         BERTParameters(){
            //HDP.SetDefault("NAME",VALUE,LOWER_LIMIT(default=-DBL_MAX),UPPER_LIMIT(default=DBL_MAX));
            HDP.SetDefault( "BERT_RADIUS_SCALE" , OLD_RADIUS_UNITS );
            HDP.SetDefault( "BERT_XSEC_SCALE" , 1.0 , 0. );
         }
   };
   BERTParameters BP;
}


// Singleton accessor

G4CascadeParameters* G4CascadeParameters::fpInstance = 0;

const G4CascadeParameters* G4CascadeParameters::Instance() {
  if (!fpInstance) {
    fpInstance = new G4CascadeParameters;
    G4AutoDelete::Register(fpInstance);
  }

  return fpInstance;
}


// Constructor initializes everything once

//#define OLD_RADIUS_UNITS (3.3836/1.2)		// Used with NucModel params

G4CascadeParameters::G4CascadeParameters()
  : G4CASCADE_VERBOSE(getenv("G4CASCADE_VERBOSE")),
    G4CASCADE_CHECK_ECONS(getenv("G4CASCADE_CHECK_ECONS")),
    G4CASCADE_USE_PRECOMPOUND(getenv("G4CASCADE_USE_PRECOMPOUND")),
    G4CASCADE_DO_COALESCENCE(getenv("G4CASCADE_DO_COALESCENCE")),
    G4CASCADE_SHOW_HISTORY(getenv("G4CASCADE_SHOW_HISTORY")),
    G4CASCADE_USE_3BODYMOM(getenv("G4CASCADE_USE_3BODYMOM")),
    G4CASCADE_USE_PHASESPACE(getenv("G4CASCADE_USE_PHASESPACE")),
    G4CASCADE_PIN_ABSORPTION(getenv("G4CASCADE_PIN_ABSORPTION")),
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
  CHECK_ECONS = (0!=G4CASCADE_CHECK_ECONS);
  USE_PRECOMPOUND = (0!=G4CASCADE_USE_PRECOMPOUND &&
		     G4CASCADE_USE_PRECOMPOUND[0]!='0');
  DO_COALESCENCE = (0==G4CASCADE_DO_COALESCENCE ||
		    G4CASCADE_DO_COALESCENCE[0]!='0');
  SHOW_HISTORY = (0!=G4CASCADE_SHOW_HISTORY);
  USE_3BODYMOM = (0!=G4CASCADE_USE_3BODYMOM);
  USE_PHASESPACE = (0!=G4CASCADE_USE_PHASESPACE &&
		    G4CASCADE_USE_PHASESPACE[0]!='0');
  PIN_ABSORPTION = (G4CASCADE_PIN_ABSORPTION ? strtod(G4CASCADE_PIN_ABSORPTION,0)
		    : 0.);
  RANDOM_FILE = (G4CASCADE_RANDOM_FILE ? G4CASCADE_RANDOM_FILE : "");
  BEST_PAR = (0!=G4NUCMODEL_USE_BEST);
  TWOPARAM_RADIUS = (0!=G4NUCMODEL_RAD_2PAR);
  RADIUS_SCALE = (G4NUCMODEL_RAD_SCALE ? strtod(G4NUCMODEL_RAD_SCALE,0)
  		  : (BEST_PAR?1.0:OLD_RADIUS_UNITS));
  if ( getenv("TEST_HDP") ) {
    HDP.DeveloperGet("BERT_RADIUS_SCALE",RADIUS_SCALE);
  }
  RADIUS_SMALL = ((G4NUCMODEL_RAD_SMALL ? strtod(G4NUCMODEL_RAD_SMALL,0)
		   : (BEST_PAR?1.992:(8.0/OLD_RADIUS_UNITS))) * RADIUS_SCALE);
  RADIUS_ALPHA = (G4NUCMODEL_RAD_ALPHA ? strtod(G4NUCMODEL_RAD_ALPHA,0)
		  : (BEST_PAR?0.84:0.70));
  RADIUS_TRAILING = ((G4NUCMODEL_RAD_TRAILING ? strtod(G4NUCMODEL_RAD_TRAILING,0)
		      : 0.) * RADIUS_SCALE);
  FERMI_SCALE = ((G4NUCMODEL_FERMI_SCALE ? strtod(G4NUCMODEL_FERMI_SCALE,0)
		  : (BEST_PAR?0.685:(1.932/OLD_RADIUS_UNITS))) * RADIUS_SCALE);
  XSEC_SCALE = (G4NUCMODEL_XSEC_SCALE ? strtod(G4NUCMODEL_XSEC_SCALE,0)
  		: (BEST_PAR?0.1:1.0) );
  if ( getenv("TEST_HDP") ) {
    HDP.DeveloperGet("BERT_XSEC_SCALE",XSEC_SCALE);
  }
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
  if (G4CASCADE_CHECK_ECONS)
    os << "G4CASCADE_CHECK_ECONS = " << G4CASCADE_CHECK_ECONS << endl;
  if (G4CASCADE_USE_PRECOMPOUND)
    os << "G4CASCADE_USE_PRECOMPOUND = " << G4CASCADE_USE_PRECOMPOUND << endl;
  if (G4CASCADE_DO_COALESCENCE)
    os << "G4CASCADE_DO_COALESCENCE = " << G4CASCADE_DO_COALESCENCE << endl;
  if (G4CASCADE_PIN_ABSORPTION)
    os << "G4CASCADE_PIN_ABSORPTION = " << G4CASCADE_PIN_ABSORPTION << endl;
  if (G4CASCADE_SHOW_HISTORY)
    os << "G4CASCADE_SHOW_HISTORY = " << G4CASCADE_SHOW_HISTORY << endl;
  if (G4CASCADE_USE_3BODYMOM)
    os << "G4CASCADE_USE_3BODYMOM = " << G4CASCADE_USE_3BODYMOM << endl;
  if (G4CASCADE_USE_PHASESPACE)
    os << "G4CASCADE_USE_PHASESPACE = " << G4CASCADE_USE_PHASESPACE << endl;
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
