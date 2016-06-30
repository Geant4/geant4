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

#include <cmath> // for log()
#include <CLHEP/Units/PhysicalConstants.h>

#include "G4MTRandGaussQ.hh"

G4MTRandGaussQ::~G4MTRandGaussQ()
{
}

G4MTRandGaussQ::G4MTRandGaussQ(const G4MTRandGaussQ& right)
  :  G4MTRandGauss(right)
{
}

G4double G4MTRandGaussQ::operator()()
{
  return transformQuick(localEngine->flat()) * defaultStdDev + defaultMean;
}

G4double G4MTRandGaussQ::operator()( G4double mean, G4double stdDev )
{
  return transformQuick(localEngine->flat()) * stdDev + mean;
}

void G4MTRandGaussQ::shootArray( const G4int size, G4double* vect,
                                 G4double mean, G4double stdDev )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = shoot(mean,stdDev);
}

void G4MTRandGaussQ::shootArray( CLHEP::HepRandomEngine* anEngine,
                            const G4int size, G4double* vect,
                            G4double mean, G4double stdDev )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = shoot(anEngine,mean,stdDev);
}

void G4MTRandGaussQ::fireArray( const G4int size, G4double* vect)
{
   for (G4int i=0; i<size; ++i)
     vect[i] = fire( defaultMean, defaultStdDev );
}

void G4MTRandGaussQ::fireArray( const G4int size, G4double* vect,
                           G4double mean, G4double stdDev )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = fire( mean, stdDev );
}

//
// Table of errInts, for use with transform(r) and quickTransform(r)
//

// Since all these are this is static to this compilation unit only, the 
// info is establised a priori and not at each invocation.

// The main data is of course the gaussQTables table; the rest is all 
// bookkeeping to know what the tables mean.

#define Table0size   250
#define Table1size  1000
#define TableSize   (Table0size+Table1size)

#define Table0step  (2.0E-6) 
#define Table1step  (5.0E-4)

#define Table0scale   (1.0/Table1step)

#define Table0offset 0
#define Table1offset (Table0size)

  // Here comes the big (5K bytes) table, kept in a file ---

static const G4float gaussTables [TableSize] = {
#include "gaussQTables.cdat"
};

G4double G4MTRandGaussQ::transformQuick (G4double r)
{
  G4double sign = +1.0; // We always compute a negative number of 
                        // sigmas.  For r > 0 we will multiply by
                        // sign = -1 to return a positive number.
  if ( r > .5 ) {
    r = 1-r;
    sign = -1.0;
  } 

  G4int index;
  G4double  dx;

  if ( r >= Table1step ) { 
    index = G4int((Table1size<<1) * r); // 1 to Table1size
    if (index == Table1size) return 0.0;
    dx = (Table1size<<1) * r - index;  // fraction of way to next bin
    index += Table1offset-1;
  } else if ( r > Table0step ) {
    G4double rr = r * Table0scale;
    index = G4int(Table0size * rr);    // 1 to Table0size
    dx = Table0size * rr - index;      // fraction of way to next bin
    index += Table0offset-1;
  } else {                             // r <= Table0step - not in tables
    return sign*transformSmall(r);
  }

  G4double y0 = gaussTables [index++];
  G4double y1 = gaussTables [index];
  
  return (G4float) (sign * ( y1 * dx + y0 * (1.0-dx) ));

} // transformQuick()

G4double G4MTRandGaussQ::transformSmall (G4double r)
{
  // Solve for -v in the asymtotic formula 
  //
  // errInt (-v) =  exp(-v*v/2)         1     1*3    1*3*5
  //                ------------ * (1 - ---- + ---- - ----- + ... )
  //                v*sqrt(2*pi)        v**2   v**4   v**6

  // The value of r (=errInt(-v)) supplied is going to less than 2.0E-13,
  // which is such that v < -7.25.  Since the value of r is meaningful only
  // to an absolute error of 1E-16 (double precision accuracy for a number 
  // which on the high side could be of the form 1-epsilon), computing
  // v to more than 3-4 digits of accuracy is suspect; however, to ensure 
  // smoothness with the table generator (which uses quite a few terms) we
  // also use terms up to 1*3*5* ... *13/v**14, and insist on accuracy of
  // solution at the level of 1.0e-7.

  // This routine is called less than one time in a million firings, so
  // speed is of no concern.  As a matter of technique, we terminate the
  // iterations in case they would be infinite, but this should not happen.

  G4double eps = 1.0e-7;
  G4double guess = 7.5;
  G4double v;
  
  for ( G4int i = 1; i < 50; i++ ) {
    G4double vn2 = 1.0/(guess*guess);
    G4double s = -13*11*9*7*5*3 * vn2*vn2*vn2*vn2*vn2*vn2*vn2;
             s +=    11*9*7*5*3 * vn2*vn2*vn2*vn2*vn2*vn2;
             s +=      -9*7*5*3 * vn2*vn2*vn2*vn2*vn2;
             s +=         7*5*3 * vn2*vn2*vn2*vn2;
             s +=          -5*3 * vn2*vn2*vn2;
             s +=             3 * vn2*vn2    - vn2  +    1.0;
    v = std::sqrt ( 2.0 * std::log ( s / (r*guess*std::sqrt(CLHEP::twopi)) ) );
    if ( std::fabs(v-guess) < eps ) break;
    guess = v;
  }
  return -v;

} // transformSmall()

#endif
