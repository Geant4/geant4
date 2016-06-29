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

#include "G4MTRandGeneral.hh"

//////////////////
// Constructors
//////////////////

G4MTRandGeneral::G4MTRandGeneral( const G4double* aProbFunc, 
                                  G4int theProbSize, 
                                  G4int IntType  )
  : deleteEngine(false), 
    nBins(theProbSize), 
    InterpolationType(IntType)
{
  localEngine = G4MTHepRandom::getTheEngine();
  prepareTable(aProbFunc);
}

G4MTRandGeneral::G4MTRandGeneral(CLHEP::HepRandomEngine& anEngine,
                                 const G4double* aProbFunc, 
                                 G4int theProbSize, 
                                 G4int IntType  )
: localEngine(&anEngine), 
  deleteEngine(false), 
  nBins(theProbSize),
  InterpolationType(IntType)
{
  prepareTable(aProbFunc);
}

G4MTRandGeneral::G4MTRandGeneral(CLHEP::HepRandomEngine* anEngine,
                                 const G4double* aProbFunc, 
                                 G4int theProbSize, 
                                 G4int IntType )
: localEngine(anEngine), 
  deleteEngine(true), 
  nBins(theProbSize),
  InterpolationType(IntType)
{
  prepareTable(aProbFunc);
}

void G4MTRandGeneral::prepareTable(const G4double* aProbFunc)
{
 //
 // Private method called only by constructors.  Prepares theIntegralPdf.
 //
  if (nBins < 1) {
    std::cerr << 
        "G4MTRandGeneral constructed with no bins - will use flat distribution\n";
    useFlatDistribution();
    return;
  }

  theIntegralPdf.resize(nBins+1);
  theIntegralPdf[0] = 0;
  G4int ptn;
  G4double weight;

  for ( ptn = 0; ptn<nBins; ++ptn ) {
    weight = aProbFunc[ptn];
    if ( weight < 0 ) {
    // We can't stomach negative bin contents, they invalidate the 
    // search algorithm when the distribution is fired.
      std::cerr << 
        "G4MTRandGeneral constructed with negative-weight bin " << ptn <<
        " = " << weight << " \n   -- will substitute 0 weight \n";
      weight = 0;
    }
    // std::cout << ptn << "  " << weight << "  " << theIntegralPdf[ptn] << "\n";
    theIntegralPdf[ptn+1] = theIntegralPdf[ptn] + weight;
  } 

  if ( theIntegralPdf[nBins] <= 0 ) {
    std::cerr << 
      "G4MTRandGeneral constructed nothing in bins - will use flat distribution\n";
    useFlatDistribution();
    return;
  }

  for ( ptn = 0; ptn < nBins+1; ++ptn ) {
    theIntegralPdf[ptn] /= theIntegralPdf[nBins];
    // std::cout << ptn << "  " << theIntegralPdf[ptn] << "\n";
  }

  // And another useful variable is ...
  oneOverNbins = 1.0 / nBins;

  // One last chore:

  if ( (InterpolationType != 0) && (InterpolationType != 1) ) {
    std::cerr << 
      "G4MTRandGeneral does not recognize IntType " << InterpolationType 
      << "\n Will use type 0 (continuous linear interpolation \n";
    InterpolationType = 0;
  }

} // prepareTable()

void G4MTRandGeneral::useFlatDistribution()
{
  //
  // Private method called only by prepareTables in case of user error. 
  //
    nBins = 1;
    theIntegralPdf.resize(2);
    theIntegralPdf[0] = 0;
    theIntegralPdf[1] = 1;
    oneOverNbins = 1.0;
    return;

} // UseFlatDistribution()

//////////////////
//  Destructor
//////////////////

G4MTRandGeneral::~G4MTRandGeneral()
{
  if ( deleteEngine ) delete localEngine;
}


///////////////////
//  mapRandom(rand)
///////////////////

G4double G4MTRandGeneral::mapRandom(G4double rand) const
{
 //
 // Private method to take the random (however it is created) and map it
 // according to the distribution.
 //
  G4int nbelow = 0;       // largest k such that I[k] is known to be <= rand
  G4int nabove = nBins;   // largest k such that I[k] is known to be >  rand
  G4int middle;
  
  while (nabove > nbelow+1) {
    middle = (nabove + nbelow+1)>>1;
    if (rand >= theIntegralPdf[middle]) {
      nbelow = middle;
    } else {
      nabove = middle;
    }
  } // after this loop, nabove is always nbelow+1 and they straddle rad:
    // assert ( nabove == nbelow+1 );
    // assert ( theIntegralPdf[nbelow] <= rand );
    // assert ( theIntegralPdf[nabove] >= rand );  
      // If a defective engine produces rand=1, that will 
      // still give sensible results so we relax the > rand assertion

  if ( InterpolationType == 1 ) {

    return nbelow * oneOverNbins;

  } else {

    G4double binMeasure = theIntegralPdf[nabove] - theIntegralPdf[nbelow];
    // binMeasure is always aProbFunc[nbelow], 
    // but we don't have aProbFunc any more so we subtract.

    if ( binMeasure == 0 ) { 
        // rand lies right in a bin of measure 0.  Simply return the center
        // of the range of that bin.  (Any value between k/N and (k+1)/N is 
        // equally good, in this rare case.)
        return (nbelow + .5) * oneOverNbins;
    }

    G4double binFraction = (rand - theIntegralPdf[nbelow]) / binMeasure;

    return (nbelow + binFraction) * oneOverNbins;
  }

} // mapRandom(rand)

 
void G4MTRandGeneral::shootArray( CLHEP::HepRandomEngine* anEngine,
                                  const G4int size, G4double* vect )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = shoot(anEngine);
}

void G4MTRandGeneral::fireArray( const G4int size, G4double* vect )
{
  for (G4int i=0; i<size; ++i)
     vect[i] = fire();
}
#endif
