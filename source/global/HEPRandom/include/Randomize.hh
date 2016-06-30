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
// $Id: Randomize.hh 96593 2016-04-25 10:12:49Z gcosmo $
//
#ifndef randomize_h
#define randomize_h 1

#include <CLHEP/Random/Randomize.h>

#if __clang__
  #if ((defined(G4MULTITHREADED) && !defined(G4USE_STD11)) || \
      !__has_feature(cxx_thread_local)) || !__has_feature(c_atomic)
    #define CLANG_NOSTDTLS
  #endif
#endif

#if (defined(G4MULTITHREADED) && \
    (!defined(G4USE_STD11) || (defined(CLANG_NOSTDTLS) || defined(__INTEL_COMPILER))))

// MT needs special Random Number distribution classes
//
#include "G4MTHepRandom.hh"
#include "G4MTRandBit.hh"
#include "G4MTRandExponential.hh"
#include "G4MTRandFlat.hh"
#include "G4MTRandGamma.hh"
#include "G4MTRandGauss.hh"
#include "G4MTRandGaussQ.hh"
#include "G4MTRandGeneral.hh"

// NOTE: G4RandStat MT-version is missing, but actually currently
// never used in the G4 source
//
#define G4RandFlat G4MTRandFlat
#define G4RandBit G4MTRandBit
#define G4RandGamma G4MTRandGamma
#define G4RandGauss G4MTRandGaussQ
#define G4RandExponential G4MTRandExponential
#define G4RandGeneral G4MTRandGeneral
#define G4Random G4MTHepRandom

#define G4UniformRand() G4MTHepRandom::getTheEngine()->flat()
//
//#include "G4UniformRandPool.hh"
//#define G4UniformRand() G4UniformRandPool::flat()
// Currently not be used in G4 source
//
#define G4RandFlatArray G4MTRandFlat::shootArray
#define G4RandFlatInt G4MTRandFlat::shootInt
#define G4RandGeneralTmp G4MTRandGeneral

#else // Sequential mode or supporting C++11 standard

// Distributions used ...
//
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Random/RandBit.h>
#include <CLHEP/Random/RandGamma.h>
#include <CLHEP/Random/RandGaussQ.h>
#include <CLHEP/Random/RandPoissonQ.h>
#include <CLHEP/Random/RandExponential.h>
#include <CLHEP/Random/RandGeneral.h>

#define G4RandStat CLHEP::HepStat
#define G4RandFlat CLHEP::RandFlat
#define G4RandBit CLHEP::RandBit
#define G4RandGamma CLHEP::RandGamma
#define G4RandGauss CLHEP::RandGaussQ
#define G4RandExponential CLHEP::RandExponential
#define G4RandGeneral CLHEP::RandGeneral
#define G4Random CLHEP::HepRandom

#define G4UniformRand() CLHEP::HepRandom::getTheEngine()->flat()

#endif // G4MULTITHREADED
#endif // randomize_h 
