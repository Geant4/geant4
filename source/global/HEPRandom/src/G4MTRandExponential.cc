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
#if __clang__
  #if ((defined(G4MULTITHREADED) && !defined(G4USE_STD11)) || \
      !__has_feature(cxx_thread_local)) || !__has_feature(c_atomic)
    #define CLANG_NOSTDTLS
  #endif
#endif

#if (defined(G4MULTITHREADED) && \
    (!defined(G4USE_STD11) || (defined(CLANG_NOSTDTLS) || defined(__INTEL_COMPILER))))

#include "G4MTRandExponential.hh"

G4MTRandExponential::~G4MTRandExponential()
{
  if ( deleteEngine ) delete localEngine;
}

G4double G4MTRandExponential::operator()()
{
  return fire( defaultMean );
}

G4double G4MTRandExponential::operator()( G4double mean )
{
  return fire( mean );
}

G4double G4MTRandExponential::shoot()
{
  return -std::log(G4MTHepRandom::getTheEngine()->flat());
}

G4double G4MTRandExponential::shoot(G4double mean)
{
  return -std::log(G4MTHepRandom::getTheEngine()->flat())*mean;
}

void G4MTRandExponential::shootArray( const G4int size, G4double* vect,
                                  G4double mean )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = shoot(mean);
}

void G4MTRandExponential::shootArray(CLHEP::HepRandomEngine* anEngine,
                      const G4int size, G4double* vect, G4double mean )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = shoot(anEngine, mean);
}

void G4MTRandExponential::fireArray( const G4int size, G4double* vect)
{
   for (G4int i=0; i<size; ++i)
     vect[i] = fire( defaultMean );
}

void G4MTRandExponential::fireArray( const G4int size, G4double* vect,
                                     G4double mean )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = fire( mean );
}

#endif
