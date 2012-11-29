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
// Defines an interface to evaporation models of Bertini cascase (BERT)
// based on INUCL code.
//
#ifndef G4ABLAEVAPORATION_h
#define G4ABLAEVAPORATION_h 1

#include "globals.hh"
#include "G4VEvaporation.hh"
#include "G4Fragment.hh"
#include "G4DynamicParticle.hh"

#include "G4Abla.hh"

//#include "G4VCoulombBarrier.hh"

//#define DEBUG

/**
 * Geant4 interface to the ABLA evaporation code.
 */

class G4AblaEvaporation : public G4VEvaporation {
public:
  /**
   * Constructor.
   */
  G4AblaEvaporation();

  /**
   * Destructor.
   */
  ~G4AblaEvaporation();

private:
  G4AblaEvaporation(const G4AblaEvaporation &right);

  const G4AblaEvaporation & operator=(const G4AblaEvaporation &right);
  G4bool operator==(const G4AblaEvaporation &right) const;
  G4bool operator!=(const G4AblaEvaporation &right) const;
  void fillResult( std::vector<G4DynamicParticle *> secondaryParticleVector,
		   G4FragmentVector * aResult );

public:

  /**
   * The method for calling
   */
  G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);
       
  void setVerboseLevel( const G4int verbose );

private:
  /**
   * Seeds for the random number generator.
   */
  G4Hazard *hazard; 

  /**
   * The verbosity of the interface class.
   */
  G4int verboseLevel;

  /**
   * Event number used in the ABLA code.
   */
  G4int eventNumber;
  
  // For Coulomb Barrier calculation
  //  G4VCoulombBarrier * theCoulombBarrierPtr;
  G4double CoulombBarrier;

#ifdef DEBUG

#endif

};

#endif
