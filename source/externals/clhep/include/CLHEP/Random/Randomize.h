// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).
//
// This file must be included to make use of the HEP Random module
// On some compilers the static instance of the HepRandom generator
// needs to be created explicitly in the client code. The static
// generator is assured to be correctly initialized by including this
// header in the client code.

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
// Gabriele Cosmo - Last change: 13th February 1996
// Ken Smith      - Added Ranshi and DualRand engines: 4th June 1998
//                - Added Ranlux64 and MTwist engines: 14th July 1998
//                - Added Hurd160, Hurd288m and TripleRand 6th Aug 1998
// =======================================================================

#ifndef Rndmze_h
#define Rndmze_h 1

// Including Engines ...

#include "CLHEP/Random/DualRand.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/MixMaxRng.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "CLHEP/Random/RanluxEngine.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "CLHEP/Random/RanshiEngine.h"

// Including distributions ...

#include "CLHEP/Random/RandBinomial.h"
#include "CLHEP/Random/RandBreitWigner.h"
#include "CLHEP/Random/RandChiSquare.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandExpZiggurat.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandBit.h"
#include "CLHEP/Random/RandGamma.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandGaussZiggurat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandLandau.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandStudentT.h"

namespace CLHEP {

#define HepUniformRand() HepRandom::getTheEngine()->flat()

// On some compilers the static instance of the HepRandom generator
// needs to be created explicitly in the client code (i.e. here).

#if __GNUC__
static const int HepRandomGenActive __attribute__((unused)) = HepRandom::createInstance();
#else
static const int HepRandomGenActive = HepRandom::createInstance();
#endif

}  // namespace CLHEP

#endif
