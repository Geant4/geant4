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
// Class G4MTRandBit
// 
// Modified version with MT extensions
// of corresponding CLHEP class RandBit

// --------------------------------------------------------------------
#ifndef G4MTRandBit_h
#define G4MTRandBit_h 1

#include "G4MTRandFlat.hh"

class G4MTRandBit : public G4MTRandFlat
{

public:

  inline G4MTRandBit ( CLHEP::HepRandomEngine& anEngine );
  inline G4MTRandBit ( CLHEP::HepRandomEngine& anEngine, G4double width );
  inline G4MTRandBit ( CLHEP::HepRandomEngine& anEngine, G4double a, G4double b );
  inline G4MTRandBit ( CLHEP::HepRandomEngine* anEngine );
  inline G4MTRandBit ( CLHEP::HepRandomEngine* anEngine, G4double width );
  inline G4MTRandBit ( CLHEP::HepRandomEngine* anEngine, G4double a, G4double b );
  // These constructors should be used to instantiate a G4MTRandBit
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the G4MTRandBit destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the G4MTRandBit destructor.

  virtual ~G4MTRandBit();
  // Destructor

  // Other than the Bit routines, constructors, and destructor, everything is
  // simply inherited from RandFlat.

  static G4int shootBit();

  static G4int shootBit( CLHEP::HepRandomEngine* );

  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  inline G4int fireBit();

};

#include "G4MTRandBit.icc"

#endif
