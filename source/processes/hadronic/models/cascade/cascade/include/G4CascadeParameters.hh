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
// $Id: G4CascadeParameters.hh 72016 2013-07-03 16:24:15Z mkelsey $
// Encapsulate all user-configurable parameters with associated envvars
//
// 20120912  M. Kelsey -- Add interface to support UI commands
// 20130304  M. Kelsey -- Add flag to collect and display cascade structure
// 20130308  M. Kelsey -- Add flag to use separate 3-body momentum generators
// 20130421  M. Kelsey -- Add flag for CHECK_ECONS, replacing #ifdef's
// 20130702  M. Kelsey -- Add flag to use N-body phase-space generator
// 20140311  G. Cosmo -- Implement standard (non-const) singleton pattern
// 20141030  M. Kelsey -- Add flag to enable direct pi-N absorption
// 20141211  M. Kelsey -- Change PIN_ABSORPTION flag to double, for energy cut

#ifndef G4CascadeParameters_hh
#define G4CascadeParameters_hh 1

#include "globals.hh"
#include <iosfwd>

class G4CascadeParamMessenger;


class G4CascadeParameters {
public:
  static const G4CascadeParameters* Instance();		// Singleton
  ~G4CascadeParameters();

  // Top-level configuration flags
  static G4int verbose()              { return Instance()->VERBOSE_LEVEL; }
  static G4bool checkConservation()   { return Instance()->CHECK_ECONS; }
  static G4bool usePreCompound()      { return Instance()->USE_PRECOMPOUND; }
  static G4bool doCoalescence()       { return Instance()->DO_COALESCENCE; }
  static G4bool showHistory()         { return Instance()->SHOW_HISTORY; }
  static G4bool use3BodyMom()	      { return Instance()->USE_3BODYMOM; }
  static G4bool usePhaseSpace()       { return Instance()->USE_PHASESPACE; }
  static G4double piNAbsorption()     { return Instance()->PIN_ABSORPTION; }
  static const G4String& randomFile() { return Instance()->RANDOM_FILE; }

  // Nuclear structure parameters
  static G4bool useTwoParam()      { return Instance()->TWOPARAM_RADIUS; }
  static G4double radiusScale()    { return Instance()->RADIUS_SCALE; }	
  static G4double radiusSmall()    { return Instance()->RADIUS_SMALL; }
  static G4double radiusAlpha()    { return Instance()->RADIUS_ALPHA; }
  static G4double radiusTrailing() { return Instance()->RADIUS_TRAILING; }
  static G4double fermiScale()     { return Instance()->FERMI_SCALE; }
  static G4double xsecScale()      { return Instance()->XSEC_SCALE; }
  static G4double gammaQDScale()   { return Instance()->GAMMAQD_SCALE; }

  // Final-state clustering cuts
  static G4double dpMaxDoublet() { return Instance()->DPMAX_DOUBLET; }
  static G4double dpMaxTriplet() { return Instance()->DPMAX_TRIPLET; }
  static G4double dpMaxAlpha()   { return Instance()->DPMAX_ALPHA; }

  static void DumpConfiguration(std::ostream& os) { Instance()->DumpConfig(os); }

private:	// Environment variable values, null pointers mean not set
  const char* G4CASCADE_VERBOSE;
  const char* G4CASCADE_CHECK_ECONS;
  const char* G4CASCADE_USE_PRECOMPOUND;
  const char* G4CASCADE_DO_COALESCENCE;
  const char* G4CASCADE_SHOW_HISTORY;
  const char* G4CASCADE_USE_3BODYMOM;
  const char* G4CASCADE_USE_PHASESPACE;
  const char* G4CASCADE_PIN_ABSORPTION;
  const char* G4CASCADE_RANDOM_FILE;
  const char* G4NUCMODEL_USE_BEST;
  const char* G4NUCMODEL_RAD_2PAR;
  const char* G4NUCMODEL_RAD_SCALE;
  const char* G4NUCMODEL_RAD_SMALL;
  const char* G4NUCMODEL_RAD_ALPHA;
  const char* G4NUCMODEL_RAD_TRAILING;
  const char* G4NUCMODEL_FERMI_SCALE;
  const char* G4NUCMODEL_XSEC_SCALE;
  const char* G4NUCMODEL_GAMMAQD;
  const char* DPMAX_2CLUSTER;
  const char* DPMAX_3CLUSTER;
  const char* DPMAX_4CLUSTER;

  void Initialize();		// Fill parameter values from envvar strings

  G4int VERBOSE_LEVEL;		// Top-level configuration flags
  G4bool CHECK_ECONS;
  G4bool USE_PRECOMPOUND;
  G4bool DO_COALESCENCE;
  G4bool SHOW_HISTORY;
  G4bool USE_3BODYMOM;
  G4bool USE_PHASESPACE;
  G4double PIN_ABSORPTION;
  G4String RANDOM_FILE;

  G4bool BEST_PAR;		// Nuclear structure parameters
//BEST_PAR has been used in a project on hold.
//Currently setting BEST_PAR or G4NUCMODEL_USE_BEST does not improve physics performance.
//Developer can get more information about this from cascade/test/README

  G4bool TWOPARAM_RADIUS;
  G4double RADIUS_SCALE;	
  G4double RADIUS_SMALL;
  G4double RADIUS_ALPHA;
  G4double RADIUS_TRAILING;
  G4double FERMI_SCALE;
  G4double XSEC_SCALE;
  G4double GAMMAQD_SCALE;

  G4double DPMAX_DOUBLET;	// Final-state clustering cuts
  G4double DPMAX_TRIPLET;
  G4double DPMAX_ALPHA;

private:	// Singleton -- no public constructor
  G4CascadeParameters();
  void DumpConfig(std::ostream& os) const;

  G4CascadeParamMessenger* messenger;		// For access via UI commands
  friend class G4CascadeParamMessenger;

  static G4CascadeParameters* fpInstance;
};

#endif	/* G4CascadeParameters_hh */
