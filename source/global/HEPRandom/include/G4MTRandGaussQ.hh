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
// Class G4MTRandGaussQ
// 
// Modified version with MT extensions
// of corresponding CLHEP class RandGaussQ

// --------------------------------------------------------------------
#ifndef G4MTRandGaussQ_h
#define G4MTRandGaussQ_h 1

#include "G4MTRandGauss.hh"

class G4MTRandGaussQ : public G4MTRandGauss
{

public:

  inline G4MTRandGaussQ ( CLHEP::HepRandomEngine& anEngine, G4double mean=0.0,
                                                G4double stdDev=1.0 );
  inline G4MTRandGaussQ ( CLHEP::HepRandomEngine* anEngine, G4double mean=0.0,
                                                G4double stdDev=1.0 );
  // These constructors should be used to instantiate a G4MTRandGaussQ
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the G4MTRandGaussQ destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the G4MTRandGaussQ destructor.

  // Destructor
  virtual ~G4MTRandGaussQ();

  //
  // Methods to generate Gaussian-distributed random deviates:
  //
  //   If a fast good engine takes 1 usec, G4MTRandGauss::fire()
  //   adds 1 usec while G4MTRandGaussQ::fire() adds only .4 usec.
  //

  // Static methods to shoot random values using the static generator

  static  inline G4double shoot();

  static  inline G4double shoot( G4double mean, G4double stdDev );

  static  void shootArray ( const G4int size, G4double* vect,
                            G4double mean=0.0, G4double stdDev=1.0 );

  //  Static methods to shoot random values using a given engine
  //  by-passing the static generator.

  static  inline G4double shoot( CLHEP::HepRandomEngine* anotherEngine );

  static  inline G4double shoot( CLHEP::HepRandomEngine* anotherEngine, 
                                  G4double mean, G4double stdDev );


  static  void shootArray ( CLHEP::HepRandomEngine* anotherEngine, 
                            const G4int size,
                            G4double* vect, G4double mean=0.0,
                            G4double stdDev=1.0 );

  //  Instance methods using the localEngine to instead of the static 
  //  generator, and the default mean and stdDev established at construction

  inline G4double fire();

  inline G4double fire ( G4double mean, G4double stdDev );
  
  void fireArray  ( const G4int size, G4double* vect);
  void fireArray  ( const G4int size, G4double* vect,
                    G4double mean, G4double stdDev );

  virtual G4double operator()();
  virtual G4double operator()( G4double mean, G4double stdDev );

protected:

  static G4double transformQuick (G4double r);
  static G4double transformSmall (G4double r);

private:

  // Private copy constructor. Defining it here disallows use.
  G4MTRandGaussQ(const G4MTRandGaussQ& d);

  // All the engine info, and the default mean and sigma,
  // are in the G4MTRandGauss base class.

};

#include "G4MTRandGaussQ.icc"

#endif
