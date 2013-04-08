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
//
// $Id:$
//
// Class G4MTHepRandom
// 
// Modified version with MT extensions
// of corresponding CLHEP class HepRandom

// --------------------------------------------------------------------
#ifndef G4MTHepRandom_h
#define G4MTHepRandom_h 1

#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandExponential.h>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Random/RandBit.h>
#include <CLHEP/Random/RandGamma.h>
#include <CLHEP/Random/RandGaussQ.h>
#include <CLHEP/Random/RandGeneral.h>

#include "tls.hh"
#include "G4Types.hh"

class G4MTHepRandom
{

public:

  G4MTHepRandom();
  G4MTHepRandom(G4long seed);
  // Contructors with and without a seed using the default engine
  // (JamesRandom).
 
  G4MTHepRandom(CLHEP::HepRandomEngine & algorithm);
  G4MTHepRandom(CLHEP::HepRandomEngine * algorithm);
  // Constructor taking an alternative engine as argument. If a pointer is
  // given the corresponding object will be deleted by the HepRandom
  // destructor.
  
  virtual ~G4MTHepRandom();
  // Destructor

  G4double flat();
  // Returns the flat value ( interval ]0...1[ ).

  void flatArray(const G4int size, G4double* vect);
  // Fills "vect" array of flat random values, given the size.

  inline G4double flat (CLHEP::HepRandomEngine* theNewEngine);
  // Returns a flat value, given a defined Random Engine.

  inline void flatArray(CLHEP::HepRandomEngine* theNewEngine, 
                        const G4int size, G4double* vect);
  // Fills "vect" array of flat random values, given the size
  // and a defined Random Engine.

  virtual G4double operator()();
  // To get a flat random number using the operator ().

  virtual std::ostream & put ( std::ostream & os ) const;
  virtual std::istream & get ( std::istream & is );
  // Save and restore to/from streams

  // --------------------------------------------------
  // Static member functions using the static generator
  // --------------------------------------------------

  static void setTheSeed(G4long seed, G4int lux=3);
  // (Re)Initializes the generator with a seed.

  static G4long getTheSeed();
  // Gets the current seed of the current generator.

  static void setTheSeeds(const G4long* seeds, G4int aux=-1);
  // (Re)Initializes the generator with a zero terminated list of seeds.

  static const G4long* getTheSeeds();
  // Gets the current array of seeds of the current generator.

  static void getTheTableSeeds (G4long* seeds, G4int index);
  // Gets the array of seeds in the static seedTable at "index" position.

  static G4MTHepRandom * getTheGenerator();
  // Return the current static generator.

  static void setTheEngine (CLHEP::HepRandomEngine* theNewEngine);
  // To set the underlying algorithm object.

  static CLHEP::HepRandomEngine * getTheEngine();
  // Returns a pointer to the underlying algorithm object.

  static void saveEngineStatus( const char filename[] = "Config.conf" );
  // Saves to file the current status of the current engine.

  static void restoreEngineStatus( const char filename[] = "Config.conf" );
  // Restores a saved status (if any) for the current engine.

  static std::ostream& saveFullState ( std::ostream & os );
  // Saves to stream the state of the engine and cached data.

  static std::istream& restoreFullState ( std::istream & is );
  // Restores from stream the state of the engine and cached data.

  static std::ostream& saveDistState ( std::ostream & os ) {return os;}
  // Saves to stream the state of the cached data.

  static std::istream& restoreDistState ( std::istream & is ) {return is;}
  // Restores from stream the state of the cached data.

  static std::ostream& saveStaticRandomStates ( std::ostream & os );
  // Saves to stream the engine and cached data for all distributions.

  static std::istream& restoreStaticRandomStates ( std::istream & is );
  // Restores from stream the engine and cached data for all distributions.

  static void showEngineStatus();
  // Dumps the current engine status on screen.

  static G4int createInstance();
  // used to initialise HepRandom::isActive and instantiate singleton

  static G4int createInstanceOnce();
  // used to initialise HepRandom::isActive and instantiate singleton

private:       // -------- Data members ---------

  static G4ThreadLocal CLHEP::HepRandomEngine * theEngine;
  // The corresponding algorithm.

  static G4ThreadLocal G4MTHepRandom * theGenerator;
  // The common shared static generator
  
  static G4ThreadLocal G4int isActive;
  // Flag notifying singleton instance
  
  G4bool deleteEngine;
  // True if the engine should be deleted on destruction.
};

std::ostream & operator<< (std::ostream & os, const G4MTHepRandom & dist);
std::istream & operator>> (std::istream & is, G4MTHepRandom & dist);

#include "G4MTHepRandom.icc"

#endif
