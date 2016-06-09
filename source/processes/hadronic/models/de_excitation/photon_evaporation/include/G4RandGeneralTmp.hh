//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- G4RandGeneralTmp ---
//                          class header file
// -----------------------------------------------------------------------

// Class defining methods for shooting generally distributed random values,
// given a user-defined probability distribution function.

// =======================================================================
// S.Magni & G.Pieri - Created: 29 April 1998 
// G.Cosmo           - Added constructor using default engine from the
//                     static generator: 20 Aug 1998
// =======================================================================

#ifndef G4RandGeneralTmp_h
#define G4RandGeneralTmp_h 1

#include "CLHEP/Random/Random.h"

class G4RandGeneralTmp : public HepRandom {

public:

  G4RandGeneralTmp ( HepDouble* aProbFunc, HepInt theProbSize );
  G4RandGeneralTmp ( HepRandomEngine& anEngine,
                HepDouble* aProbFunc, HepInt theProbSize );
  G4RandGeneralTmp ( HepRandomEngine* anEngine, 
                HepDouble* aProbFunc, HepInt theProbSize );
  // These constructors should be used to instantiate a G4RandGeneralTmp
  // distribution object defining a local engine for it.
  // The static generator will be skeeped using the non-static methods
  // defined below. In case no engine is specified in the constructor, the
  // default engine used by the static generator is applied.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the G4RandGeneralTmp destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the RandGauss destructor.
  // The probability distribution function (Pdf) must be provided by the user
  // as an array of positive real number. The array size must also be
  // provided. The Pdf doesn't need to be normalized to 1. 

  virtual ~G4RandGeneralTmp();
  // Destructor

  // Methods to shoot random values using the static generator
  // N.B.: The methods are NOT static since they use nonstatic members
  // theIntegralPdf & nBins

  inline HepDouble shoot();

  inline void shootArray ( const HepInt size, HepDouble* vect);

  //  Methods to shoot random values using a given engine
  //  by-passing the static generator.

  HepDouble shoot( HepRandomEngine* anEngine );

  void shootArray ( HepRandomEngine* anEngine, const HepInt size,
                    HepDouble* vect );
			    
  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  HepDouble fire();

  void fireArray ( const HepInt size, HepDouble* vect);

  HepDouble operator()();

private:

  // Private copy constructor. Defining it here disallows use.
  G4RandGeneralTmp(const G4RandGeneralTmp&) : HepRandom() {;}

  HepRandomEngine* localEngine;
  HepBoolean deleteEngine;
  HepDouble* theIntegralPdf;
  HepInt nBins;

};

#include "G4RandGeneralTmp.icc"

#endif
