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
//                          --- RandGeneralTmp ---
//                      class implementation file
// -----------------------------------------------------------------------

// Class defining methods for shooting generally distributed random values,
// given a user-defined probability distribution function.

// =======================================================================
// S.Magni & G.Pieri - Created: 29th April 1998
// G.Cosmo           - Added constructor using default engine from the
//                     static generator. Simplified shoot() and
//                     shootArray() (not needed in principle!): 20 Aug 1998
// =======================================================================

#include "globals.hh"
#include "G4RandGeneralTmp.hh"

//////////////////
// Constructors
//////////////////

G4RandGeneralTmp::G4RandGeneralTmp(HepDouble* aProbFunc, HepInt theProbSize )
: deleteEngine(false), nBins(theProbSize)
{
  localEngine = HepRandom::getTheEngine();
  register HepInt ptn;
  theIntegralPdf = new HepDouble[theProbSize];
  theIntegralPdf[0] = 0;
  for ( ptn = 1; ptn<theProbSize; ++ptn ) 
  {
    theIntegralPdf[ptn] = theIntegralPdf[ptn-1] + aProbFunc[ptn];
  } 
  for ( ptn = 0; ptn < theProbSize; ++ptn )
  { 
    if (theIntegralPdf[nBins-1] != 0.) theIntegralPdf[ptn] /= theIntegralPdf[nBins-1];
  }
}

G4RandGeneralTmp::G4RandGeneralTmp(HepRandomEngine& anEngine,
                         HepDouble* aProbFunc, HepInt theProbSize )
: localEngine(&anEngine), deleteEngine(false), nBins(theProbSize)
{
  register HepInt ptn;
  theIntegralPdf = new HepDouble[theProbSize];
  theIntegralPdf[0] = 0;
  for ( ptn = 1; ptn<theProbSize; ++ptn ) 
  {
    theIntegralPdf[ptn] = theIntegralPdf[ptn-1] + aProbFunc[ptn];
  } 
  for ( ptn = 0; ptn < theProbSize; ++ptn )
  {
    theIntegralPdf[ptn] /= theIntegralPdf[nBins-1];
  }
}

G4RandGeneralTmp::G4RandGeneralTmp(HepRandomEngine* anEngine,
                         HepDouble* aProbFunc, HepInt theProbSize )
: localEngine(anEngine), deleteEngine(true), nBins(theProbSize)
{
  register HepInt ptn;
  theIntegralPdf = new HepDouble[nBins];
  theIntegralPdf[0] = 0.;
  for ( ptn = 1; ptn<nBins; ++ptn ) {
    theIntegralPdf[ptn] = 0;
    theIntegralPdf[ptn] = theIntegralPdf[ptn-1] + aProbFunc[ptn];
  } 
  for ( ptn = 0; ptn < nBins; ++ptn )
  theIntegralPdf[ptn] /=  theIntegralPdf[nBins-1];


}

//////////////////
//  Destructor
//////////////////

G4RandGeneralTmp::~G4RandGeneralTmp() {
  if ( deleteEngine ) delete localEngine;
  delete [] theIntegralPdf;
}

HepDouble G4RandGeneralTmp::operator()() {
  return fire();
}

HepDouble G4RandGeneralTmp::shoot( HepRandomEngine* anEngine )
{
  HepDouble rand;
  HepInt nabove, nbelow = 0, middle;
  
  nabove = nBins+1;  
  rand = anEngine->flat();
  
  while(nabove-nbelow > 1) {
    middle = ( nabove + nbelow ) / 2;
    if (rand == theIntegralPdf[middle-1]) break;
    if (rand < theIntegralPdf[middle-1]) nabove = middle;
    else nbelow = middle;
  }
  
  return ((HepDouble)nbelow - 1) / nBins;
}

void G4RandGeneralTmp::shootArray( HepRandomEngine* anEngine,
                              const HepInt size, HepDouble* vect )
{
   register HepInt i;

   for (i=0; i<size; ++i)
     vect[i] = shoot(anEngine);
}

HepDouble G4RandGeneralTmp::fire()
{
  HepDouble rand;
  HepInt nabove, nbelow = 0, middle;
  
  nabove = nBins+1;  
  rand = localEngine->flat();
  
  while(nabove-nbelow > 1) {
    middle = ( nabove + nbelow ) / 2;
    if (rand == theIntegralPdf[middle-1]) break;
    if (rand < theIntegralPdf[middle-1]) nabove = middle;
    else nbelow = middle;
  }
  
  return ((HepDouble)nbelow - 1) / nBins;

}

void G4RandGeneralTmp::fireArray( const HepInt size, HepDouble* vect )
{
   register HepInt i;

   for (i=0; i<size; ++i)
     vect[i] = fire();
}
