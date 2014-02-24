// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandGeneral ---
//                      class implementation file
// -----------------------------------------------------------------------

// =======================================================================
// S.Magni & G.Pieri - Created: 5th September 1995
// G.Cosmo           - Added constructor using default engine from the
//                     static generator. Simplified shoot() and
//                     shootArray() (not needed in principle!): 20th Aug 1998
// M.G.Pia & G.Cosmo - Fixed bug in computation of theIntegralPdf in
//                     two constructors: 5th Jan 1999
// S.Magni & G.Pieri - Added linear interpolation: 24th Mar 1999
// M. Fischler	     - General cleanup: 14th May 1999
//			+ Eliminated constructor code replication by factoring 
//			  common code into prepareTable.
//			+ Eliminated fire/shoot code replication by factoring 
//			  out common code into mapRandom.  
//			+ A couple of methods are moved inline to avoid a 
//			  speed cost for factoring out mapRandom:  fire()
//			  and shoot(anEngine).
//			+ Inserted checks for negative weight and zero total 
//			  weight in the bins.
//			+ Modified the binary search loop to avoid incorrect
//			  behavior when rand is example on a boundary.
//			+ Moved the check of InterpolationType up into 
//			  the constructor.  A type other than 0 or 1
//			  will give the interpolated distribution (instead of
//			  a distribution that always returns 0).
//			+ Modified the computation of the returned value
//			  to use algeraic simplification to improve speed.
//			  Eliminated two of the three divisionns, made
//			  use of the fact that nabove-nbelow is always 1, etc.
//			+ Inserted a check for rand hitting the boundary of a
//			  zero-width bin, to avoid dividing 0/0.  
// M. Fischler	      - Minor correction in assert 31 July 2001
//			+ changed from assert (above = below+1) to ==
// M Fischler         - put and get to/from streams 12/15/04
//			+ Modifications to use a vector as theIntegraPdf
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
//
// =======================================================================

#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/DoubConv.h"
#include <cassert>

namespace CLHEP {

std::string RandGeneral::name() const {return "RandGeneral";}
HepRandomEngine & RandGeneral::engine() {return *localEngine;}


//////////////////
// Constructors
//////////////////

RandGeneral::RandGeneral( const double* aProbFunc, 
			  int theProbSize, 
			  int IntType  )
  : HepRandom(),
    localEngine(HepRandom::getTheEngine(), do_nothing_deleter()),
    nBins(theProbSize), 
    InterpolationType(IntType)
{
  prepareTable(aProbFunc);
}

RandGeneral::RandGeneral(HepRandomEngine& anEngine,
                         const double* aProbFunc, 
			 int theProbSize, 
			 int IntType  )
: HepRandom(),
  localEngine(&anEngine, do_nothing_deleter()), 
  nBins(theProbSize),
  InterpolationType(IntType)
{
  prepareTable(aProbFunc);
}

RandGeneral::RandGeneral(HepRandomEngine* anEngine,
                         const double* aProbFunc, 
			 int theProbSize, 
			 int IntType )
: HepRandom(),
  localEngine(anEngine), 
  nBins(theProbSize),
  InterpolationType(IntType)
{
  prepareTable(aProbFunc);
}

void RandGeneral::prepareTable(const double* aProbFunc) {
//
// Private method called only by constructors.  Prepares theIntegralPdf.
//
  if (nBins < 1) {
    std::cerr << 
	"RandGeneral constructed with no bins - will use flat distribution\n";
    useFlatDistribution();
    return;
  }

  theIntegralPdf.resize(nBins+1);
  theIntegralPdf[0] = 0;
  int ptn;
  double weight;

  for ( ptn = 0; ptn<nBins; ++ptn ) {
    weight = aProbFunc[ptn];
    if ( weight < 0 ) {
    // We can't stomach negative bin contents, they invalidate the 
    // search algorithm when the distribution is fired.
      std::cerr << 
	"RandGeneral constructed with negative-weight bin " << ptn <<
	" = " << weight << " \n   -- will substitute 0 weight \n";
      weight = 0;
    }
    // std::cout << ptn << "  " << weight << "  " << theIntegralPdf[ptn] << "\n";
    theIntegralPdf[ptn+1] = theIntegralPdf[ptn] + weight;
  } 

  if ( theIntegralPdf[nBins] <= 0 ) {
    std::cerr << 
      "RandGeneral constructed nothing in bins - will use flat distribution\n";
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
      "RandGeneral does not recognize IntType " << InterpolationType 
      << "\n Will use type 0 (continuous linear interpolation \n";
    InterpolationType = 0;
  }

} // prepareTable()

void RandGeneral::useFlatDistribution() {
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

RandGeneral::~RandGeneral() {
}


///////////////////
//  mapRandom(rand)
///////////////////

double RandGeneral::mapRandom(double rand) const {
//
// Private method to take the random (however it is created) and map it
// according to the distribution.
//

  int nbelow = 0;	  // largest k such that I[k] is known to be <= rand
  int nabove = nBins;     // largest k such that I[k] is known to be >  rand
  int middle;
  
  while (nabove > nbelow+1) {
    middle = (nabove + nbelow+1)>>1;
    if (rand >= theIntegralPdf[middle]) {
      nbelow = middle;
    } else {
      nabove = middle;
    }
  } // after this loop, nabove is always nbelow+1 and they straddle rad:
    assert ( nabove == nbelow+1 );
    assert ( theIntegralPdf[nbelow] <= rand );
    assert ( theIntegralPdf[nabove] >= rand );  
		// If a defective engine produces rand=1, that will 
		// still give sensible results so we relax the > rand assertion

  if ( InterpolationType == 1 ) {

    return nbelow * oneOverNbins;

  } else {

    double binMeasure = theIntegralPdf[nabove] - theIntegralPdf[nbelow];
    // binMeasure is always aProbFunc[nbelow], 
    // but we don't have aProbFunc any more so we subtract.

    if ( binMeasure == 0 ) { 
	// rand lies right in a bin of measure 0.  Simply return the center
	// of the range of that bin.  (Any value between k/N and (k+1)/N is 
	// equally good, in this rare case.)
        return (nbelow + .5) * oneOverNbins;
    }

    double binFraction = (rand - theIntegralPdf[nbelow]) / binMeasure;

    return (nbelow + binFraction) * oneOverNbins;
  }

} // mapRandom(rand)
 
void RandGeneral::shootArray( HepRandomEngine* anEngine,
                            const int size, double* vect )
{
   int i;

   for (i=0; i<size; ++i) {
     vect[i] = shoot(anEngine);
   }
}

void RandGeneral::fireArray( const int size, double* vect )
{
   int i;

   for (i=0; i<size; ++i) {
      vect[i] = fire();
   }
}

std::ostream & RandGeneral::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  os << nBins << " " << oneOverNbins << " " << InterpolationType << "\n";    
  t = DoubConv::dto2longs(oneOverNbins);
  os << t[0] << " " << t[1] << "\n";
  assert (static_cast<int>(theIntegralPdf.size())==nBins+1);
  for (unsigned int i=0; i<theIntegralPdf.size(); ++i) {
    t = DoubConv::dto2longs(theIntegralPdf[i]);
    os << theIntegralPdf[i] << " " << t[0] << " " << t[1] << "\n";
  }
  os.precision(pr);
  return os;
}

std::istream & RandGeneral::get ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != name()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read state of a "
    	      << name() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  if (possibleKeywordInput(is, "Uvec", nBins)) {
    std::vector<unsigned long> t(2);
    is >> nBins >> oneOverNbins >> InterpolationType;
    is >> t[0] >> t[1]; oneOverNbins = DoubConv::longs2double(t); 
    theIntegralPdf.resize(nBins+1);
    for (unsigned int i=0; i<theIntegralPdf.size(); ++i) {
      is >> theIntegralPdf[i] >> t[0] >> t[1];
      theIntegralPdf[i] = DoubConv::longs2double(t); 
    }
    return is;
  }
  // is >> nBins encompassed by possibleKeywordInput
  is >> oneOverNbins >> InterpolationType;
  theIntegralPdf.resize(nBins+1);
  for (unsigned int i=0; i<theIntegralPdf.size(); ++i) is >> theIntegralPdf[i];
  return is;
}

}  // namespace CLHEP
